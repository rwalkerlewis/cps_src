        program sacamp
        implicit none
c-----
c     give a velocity model
c     compute the site response using
c     The Boore and Brown quarter wavelength technique. 
c     the coding follows the siate_amp_batch.for of Boore
c-----
        integer NLAY
        parameter (NLAY=200)
        common/isomod/d(NLAY),va(NLAY),vb(NLAY),rho(NLAY),
     1      qa(NLAY),qb(NLAY),etap(NLAY),etas(NLAY),
     2      frefp(NLAY), frefs(NLAY)
        real d,va,vb,rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep
        common/modtit/title
        character*80 title

c-----
c       internal variables
c-----
        integer*4 mmax, iunit, iiso, iflsph, idimen, icnvel
        integer*4 nerr
        logical listmd
        character mname*40
        logical dop, outtxt, dorayp
        integer i
        real angdeg, rayp
        real fmin, fmax, freq, df
        integer nfreq
        real tmark

        real v(NLAY), z(NLAY), q(NLAY), dens(NLAY)
        real t(NLAY), tq(NLAY)
        real xx, yy

c-----
c       arrays for conversion to Sac
c       x() is the amplitude response
c       y() is amplitude * exp( - pi freq kappa) 
c-----
        integer NFQ
        parameter (NFQ=1000)
        real x(NFQ), y(NFQ)
c-----
c       get the command line arguments
c-----
        call gcmdln(mname,dop,angdeg,rayp,outtxt,dorayp)
        if(mname.eq.' ')then
           call usage('No model file given')
           stop
        endif
c-----
c       get the model, convert Q (defined by values > 1) to q inverse
c-----
        call getmod(1,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,nerr,listmd)
        do i=1,mmax
           if(qa(i).gt.1.0)then
              qa(i) = 1.0/qa(i)
           endif
           if(qb(i).gt.1.0)then
              qb(i) = 1.0/qb(i)
           endif
        enddo
c-----
c       do the processing
c-----

c-----
c       create velocity versus depth 
c       z    - depth
c       vel  - velocity
c       dens -density
c       t    - one-way travel time
c       tq   - one way t*
c----
         
        do i=1,mmax
            if(dop)then
                  v(i)   = va(i)
                  q(i)   = qa(i)
            else
                  v(i)   = vb(i)
                  q(i)   = qb(i)
            endif
            if(i.eq.1)then
               z(i)  = 0.0
               t(i)  = 0.0
               tq(i) = 0.0
            else
               t(i)  = t(i-1)  + d(i-1)/v(i-1)
               tq(i) = tq(i-1) + q(i-1)*d(i-1)/v(i-1)
               z(i)  = z(i-1)  + d(i-1)
            endif
            dens(i) = rho (i)
        enddo
c-----
c       if ray parameter given, then define angle of incidence in
c       bottom layer
c-----
        if(dorayp)then
             angdeg = asin(rayp*vb(mmax))*180.0/3.1415927
        endif
c-----
c       loop over frequencies
c-----
        if(outtxt)then
            open(8,file='sacampl.txt',access='sequential',
     1         status='unknown', form='formatted')
            rewind 8
            write(8,1)
    1   format( ' Freq(Hz)    Amp     Amp(kappa)')
        endif
        fmin = 0.1
        fmax = 20.0
        nfreq = 200
        df = (fmax-fmin)/(nfreq-1)
 
        do i=1,nfreq,1
           freq = fmin +(i-1)*df
c-----
c          according to Boore and Brown SRL 1998, f = 1/[4*OneWaytravelTime]
c-----
           tmark = 1.0/(4.* freq)
c-----
c          now search the arrays
c          t tq to get average velocity,average angle, kappa and 
c-----
           call dosearch(freq,mmax,t,tq,z,dens,tmark,vb(mmax),
     1         rho(mmax),angdeg,xx,yy)
           x(i) = xx
           y(i) = yy
           if(outtxt)then
              write(8,2)freq,xx,yy
           endif
        enddo
    2   format(3(1pe10.3))
        if(outtxt)then
           close(8)
        endif
c-----
c       output sac files for the amplitude and
c       amplitude with kappa response
c-----
        call outsac(fmin,fmax,df,nfreq,x,'AMP.sac.am' )
        call outsac(fmin,fmax,df,nfreq,y,'KAMP.sac.am')
        end

        subroutine outsac(fmin,fmax,df,nfreq,x,fname)
        implicit none
        real fmin, fmax, df
        integer nfreq
        real x(*)
        character fname*(*)

        real depmin,depmax, depmen
        integer indmax, indmin
        integer nerr
        real e, btime


        call newhdr()
        call scmxmn(x,nfreq,depmax,depmin,depmen,indmax,indmin)
        call setfhv('DEPMAX', depmax, nerr)
        call setfhv('DEPMIN', depmin, nerr)
        call setfhv('DEPMEN', depmen, nerr)
        call setnhv('NPTS    ',nfreq	,nerr)
        call setfhv('DELTA   ',df  ,nerr)
        call setfhv('B       ',fmin  ,nerr)
        call setfhv('TIMMAX  ',fmin + indmax*df  ,nerr)
        call setfhv('TIMMIN  ',fmin + indmin*df  ,nerr)
c----
c       this is a hack since Sac does not have a specific
c       value to represents frequency
c       ITIME which happens to be 1   denotes a time series
c       call setihv('IFTYPE  ','ITIME   ',nerr)
c----
        call setihv('IFTYPE  ','IAMPH   ',nerr)
        e = btime + (nfreq -1 )*df
        call setfhv('E       ',e     ,nerr)
        call setlhv('LEVEN   ',.true.,nerr)
        call setlhv('LOVROK  ',.true.,nerr)
        call setlhv('LCALDA  ',.true.,nerr)
        call bwsac(1,nfreq,fname,x)
        return
        end


        subroutine dosearch(freq,mmax,t,tq,z,dens,tmark,vbase,
     1         rhobase,angdeg,amp,ka)
        implicit none
        real freq
        integer mmax
        real t(mmax), tq(mmax), z(mmax), dens(mmax)
        real tmark
        real vbase, rhobase, angdeg
        real amp, ka

        integer i,k
        real p
        real tt, ttq, zz, rr
        real vhat, coshat
        real fac
        real deg, cg,sg
c-----
c       assume that the t, tq and z are continuous
c       rather than discrete layers
c-----
        deg = 3.1415927*angdeg/180.0
        cg = cos(deg)
        sg = sin(deg)
            
        amp = 0.0
        ka = 0.0
c-----
c       if the desired travel time is greater then that 
c       to the base of layer mmax-1, then a special case is made
c       to consider that layer. However the Q in that lower
c       halfspace is infinite
c-----
        if(tmark.gt.t(mmax))then
c-----
c           special care
c-----
            tt = tmark
            zz = z(mmax) + vbase*(tmark-t(mmax))
            rr = dens(mmax)
            vhat = zz/tt
            coshat=sqrt(1. - (vhat*sg/vbase)**2 )
            fac = rhobase*vbase*cg/(rr*vhat*coshat)
            amp = sqrt(fac)
            ka = amp*exp(-3.1415927*freq*tq(mmax))
            ttq = tq(mmax)
        else
c-----
c           the 1/4 wavelength depth is in the layer stack
c-----

            do i=1,mmax-1
               if(tmark.ge.t(i) .and. tmark.le.t(i+1))then
                     p = (tmark-t(i))/(t(i+1)-t(i))
                     tt  = (1-p)*t(i)  + p*t(i+1)
                     ttq = (1-p)*tq(i) + p*tq(i+1)
                     zz  = (1-p)*z(i)  + p*z(i+1)
                     rr  = (1-p)*dens(i)  + p*dens(i+1)
                     vhat = zz/tt
                     coshat=sqrt(1. - (vhat*sg/vbase)**2 )
                     fac = rhobase*vbase*cg/(rr*vhat*coshat)
                     amp = sqrt(fac)
c-----
c      use the kappa from the deepest layer
c-----
                     ka = amp*exp(-3.1415927*freq*tq(mmax))
               endif
            enddo
        endif
               WRITE(6,*)freq,amp,ttq,ka
        return
        end
        subroutine gcmdln(mname,dop,angdeg,rayp,outtxt,dorayp)
        implicit none
c-----
c       parse the command line arguments
c-----
c       mname   C*80    - name of model file
c       dop     L       - .true.  get P
c                         .false. get S
c       angdeg  R       - angle in incident S wave at base degrees
c-----
        character mname*(*)
        logical dop , outtxt, dorayp
        real angdeg, rayp
c-----
c       internal variables
c-----
        character name*40
        integer mnmarg
        integer i, nmarg
c-----
c       set defaults
c-----
        mname = ' '
        dop  = .false.
        outtxt = .false.
        dorayp = .false.
        rayp = -1.0
        angdeg = -1.0
        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
            call mgtarg(i,name)
            if(name(1:2).eq.'-M')then
                 i = i + 1
                 call mgtarg(i,mname)
            else if(name(1:2).eq.'-P')then
                 dop  = .true.
            else if(name(1:2).eq.'-S')then
                 dop  = .false.
            else if(name(1:4).eq.'-TXT')then
                 outtxt  = .true.
            else if(name(1:2).eq.'-A')then
                 i = i + 1
                 call mgtarg(i,name)
                 read(name,'(bn,f20.0)')angdeg
            else if(name(1:5).eq.'-RAYP')then
                 i = i + 1
                 call mgtarg(i,name)
                 read(name,'(bn,f20.0)')rayp
                 dorayp = .true.
            else if(name(1:2).eq.'-?')then
                 call usage(' ')
            else if(name(1:2).eq.'-h')then
                 call usage(' ')
            endif
         go to 1000
 2000    continue
c-----
c        check
c----- 
         if(angdeg.lt.0.0 .and. rayp.lt.0.0)then
             call usage(' -A angdeg or -P ray_parameter required')
         endif
         return
         end

        subroutine usage(str)
        implicit none
     
        integer LER
        parameter (LER=0)

        character str*(*)

        if(str.ne.' ')then
            WRITE(LER,*)str
        endif
        write(LER,*)'Usage: sacamp -M model  ',
     1      ' [-P | -S ] [-A angdeg | -AYPay_parameter]',
     2       ' [-TXT] [-?] [-h]'
        write(LER,*)
     1  '-M model   (default none )  Earth model file'
        write(LER,*)
     1  '-P         (default false )  P response'
        write(LER,*)
     1  '-S         (default true  )  S response'
        write(LER,*)
     1  '-A  angdeg (default none )  S angle incidence at base'
        write(LER,*)
     1  '-RAYP  ray_parameter (default none) (s/km) if model in km'
        write(LER,*)
     1  '       NOTE: Either the angle or the ray parameter is required'
        write(LER,*)
     1  '-TXT       (default false) Output in file sacamp.txt'
        write(LER,*)
     1  '-?         (default none )  this help message '
        write(LER,*)
     1  '-h         (default none )  this help message '
        stop
        end
          
