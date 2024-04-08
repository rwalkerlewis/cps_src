        program slegn96
c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME III                                                      c
c                                                                      c
c      PROGRAM: SLEGN96                                                c
c                                                                      c
c      COPYRIGHT 1996                                                  c
c      R. B. Herrmann                                                  c
c      Department of Earth and Atmospheric Sciences                    c
c      Saint Louis University                                          c
c      221 North Grand Boulevard                                       c
c      St. Louis, Missouri 63103                                       c
c      U. S. A.                                                        c
c                                                                      c
c----------------------------------------------------------------------c
c
c
c     This program calculates the group velocity and partial
c     derivatives of Love waves for any plane multi-layered
c     model.  The propagator-matrix, instead of numerical-
c     integration method is used, in which the Haskell rather
c     than Harkrider formalisms are concerned.
c
c     Developed by C. Y. Wang and R. B. Herrmann, St. Louis
c     University, Oct. 10, 1981.  Modified for use in surface
c     wave inversion, with addition of spherical earth flattening
c     transformation and numerical calculation of group velocity
c     partial derivatives by David R. Russell, St. Louis
c     University, Jan. 1984.
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c Revision history:
c       07 AUG 2002 - make string lengths 120 characters from 80 
c       13 OCT 2006 - verbose output of energy integrals
c       26 SEP 2008 - fixed undefined LOT in subroutine up
c       14 JUN 2009 - reformulate in general terms as per book
c       29 JUL 2017 - modify insert if depth is in halfspace
c       22 FEB 2018 - verbose also gives LAGR/(omega^2 I0)
c       19 MAR 2022 - corrected some numebrical problems
c           This arose when using a model in cm, cm/s and gm/cm^3 
c           instead of km, km,s and km/cm^3.  Numerical precision was
c           lost. Statements such as
c           if(dabs(pr) .lt. 1.0e-5 .and. cdabs(rp).lt.1.0e-5)then
c           were changes to use a much lower limit since double precision
c           permitted this.
c           Note that was only a problem with CGS units and not MKS or
c           the km, km/s gm/cm^3 mixed units
c       09 DEC 2023 - change format statement 2 from AR= to AL=
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        parameter(LER=0,LIN=5,LOT=6)
        parameter(NL=200)
c-----
c       LIN - unit for FORTRAN read from terminal
c       LOT - unit for FORTRAN write to terminal
c       LER - unit for FORTRAN error output to terminal
c       NL  - number of layers in model
c-----
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d,a,b,rho,qa,qb,etap,etas,frefp,frefs

        common/modlly/mmax
        integer mmax

        common/eigfun/ uu(NL),tt(NL),dcdb(NL),dcdh(NL),dcdr(NL),uu0(4)
        real*8 uu, tt, dcdb, dcdh, dcdr, uu0
        real*4 sdcda(NL), sdcdb(NL), sdcdh(NL), sdcdr(NL) 
        real*4 spur(NL), sptr(NL), spuz(NL), sptz(NL)
        common/wateri/iwat(NL)
        common/sumi/   sumi0,sumi1,sumi2,flagr,ale,ugr
        real*8 sumi0, sumi1, sumi2, flagr, ale, ugr
        common/depref/refdep
        real refdep

        real*4 vtp,dtp,rtp
        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        real*4 s1
        parameter(MAXMOD=2000)
        real*8 cp(MAXMOD), t
        real*8 c, freq, omega, wvno, gammal, csph, usph
        logical nwlyrs, nwlyrr
        character*12 fname(2)
c-----
c       wired in file names     - output of this program 
c                           is always in the
c               binary file slegn96.egn or slegn96.der
c                               - input from sdisp96 is always
c               binary file sdisp96.lov
c-----
        logical dolove, dorayl
        character hsfile*120, hrfile*120, title*120 
        logical dotmp
        logical ext
        logical dogam
        logical dderiv

        character mname*120
        integer ipar(10)
        real*4 fpar(10)

        integer nipar(10)

        common/verby/verbose, dout, doauto
        logical verbose, dout, doauto
        real*8 ccausal, gam


        common/control/jref
        integer jref
        integer lmin

c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line information
c-----  
        call gcmdln(hsfile,hrfile,hs,hr,dotmp,dogam,dderiv,
     1      nipar,verbose,doauto)
        if(dotmp)then
            fname(1) = 'tsdisp96.lov'
        else
            fname(1) = 'sdisp96.lov'
        endif
        if(dderiv)then
            fname(2) = 'slegn96.der'
        else
            fname(2) = 'slegn96.egn'
        endif
c-----
c       get control parameters from sdisp96.dat
c-----
        inquire(file='sdisp96.dat',exist=ext)
        if(.not.ext)then
            call usage('Control file sdisp96.dat'//
     1          ' does not exist')
        endif
        open(1,file='sdisp96.dat',form='formatted',status='unknown',
     1      access='sequential')
        rewind 1
        read(1,*,end=9999,err=9999)dt
        read(1,*,end=9999,err=9999)npts,n1,n2
        read(1,'(a)',end=9999,err=9999)mname
        read(1,*)dolove, dorayl
        read(1,*)ohs, ohr
        read(1,*)nmodes
        read(1,*)faclov,facray
        close(1)
c-----
c       get the earth model
c-----
        inquire(file=mname,exist=ext)
        if(.not. ext)call usage('Model file does not exist')
        l = lgstr(mname)

        write(LOT,*)'Model name: ',mname(1:l)
        call getmod(1,mname,mmax,title,iunit,iiso,iflsph,
     1      idimen,icnvel,ierr,.true.)
        nsph = iflsph
        if(nsph.gt.0)then
                call bldsph(mmax)
        endif
c-----
c       get source depth, and sphericity control
c-----
        if(hs.lt.-1.0E+20)then
            depths = ohs
        else
            depths = hs
        endif
        if(hr.lt.-1.0E+20)then
            depthr = ohr
        else
            depthr = hr
        endif
        depthr = depthr + refdep
        depths = depths + refdep
c-----
c       if there is a spherical model, map the source 
c            and receiver depth 
c       in the spherical model to the equivalent depth in the flat model
c-----
        if(nsph.gt.0)then
                dephs = 6371.0 *alog( 6371.0/(6371.0 - (depths-refdep)))
                dephr = 6371.0 *alog( 6371.0/(6371.0 - (depthr-refdep)))
        else
                dephs = depths
                dephr = depthr
        endif
c----
c       see if the file sdisp96.lov exists
c-----
        inquire(file=fname(1),exist=ext)
        if(.not.ext)then
            call usage('Dispersion file: '//fname(1)//
     1          ' does not exist')
        endif
c-----
c       open output file of sdisp96 and of slegn96
c-----
        open(1,file=fname(1),form='unformatted',status='unknown',
     1          access='sequential')
        open(2,file=fname(2),form='unformatted',status='unknown',
     1          access='sequential')
        rewind 1
        rewind 2
c-----
c       obtain the earth model: note if the original was spherical,
c       this will be the transformed model
c-----
        call gtsmdl(1,mmax,d,a,b,rho,qa,qb,nper,
     2              mname,ipar,fpar)
c-----
c       get the Q information into the program
c-----
        do i=1,mmax
            if(.not.dogam)then
                qa(i) = 0.0
                qb(i) = 0.0
            endif
            if(qa(i).gt.1.0)qa(i)=1.0/qa(i)
            if(qb(i).gt.1.0)qb(i)=1.0/qb(i)
            if(frefp(i).le.0.0)frefp(i) = 1.0
            if(frefs(i).le.0.0)frefs(i) = 1.0
        enddo
c-----
c       define modulus of rigidity, also get current transformed model
c       parameters
c       set fpar(1) = refdep
c       set ipar(1) = 1 if medium is spherical
c       set ipar(2) = 1 if source is in fluid
c       set ipar(3) = 1 if receiver is in fluid
c       set ipar(4) = 1 if eigenfunctions are output with -DER flag
c       set ipar(5) = 1 if dc/dh are output with -DER flag
c       set ipar(6) = 1 if dc/da are output with -DER flag
c       set ipar(7) = 1 if dc/db are output with -DER flag
c       set ipar(8) = 1 if dc/dr are output with -DER flag
c-----
        deplw = 0.0
        depup = 0.0
        ipar(2) = 0
        ipar(3) = 0
        do 185 i=1,mmax
            depup = deplw + d(i)
            if(b(i) .lt. 0.0001*a(i))then
                iwat(i) = 1
                if(depths.ge.deplw .and. depths.lt.depup)then
                    ipar(2) = 1
                endif
                if(depthr.ge.deplw .and. depthr.lt.depup)then
                    ipar(3) = 1
                endif
            else
                iwat(i) = 0
            endif
            deplw = depup
  185   continue
        do 186 i=4,8
            ipar(i) = nipar(i)
  186   continue
c-----
c               Insert source and receiver depths into the model used
c               for the dispersion. Note that if the model is spherical
c               we have already mapped the spehrical depth into
c               the flattened model.
c-----
        call putmdl(2,mmax,d,a,b,rho,qa,qb,nper,depths-refdep,
     2      depthr-refdep,mname,ipar,fpar)
c-----
c       split a layer at the source depth
c       the integer ls gives the interface at which the eigenfunctions
c       are to be defined
c-----
C        do  ii=1,mmax
C            write(6,*)'before:',iwat(ii),d(ii),a(ii),b(ii),rho(ii)
C        enddo
        call insert(nsph,dephs,nwlyrs,lsso)
        call insert(nsph,dephr,nwlyrr,lrro)
        call srclyr(dephs,lss)
        call srclyr(dephr,lrr)
c-----DEBUG
c       output the new model with source layer
c-----
C       write(6,*)'lss,lrr:',lss,lrr
C        do  ii=1,mmax
C            write(6,*)'after :',iwat(ii),d(ii),a(ii),b(ii),rho(ii)
C        enddo
c-----
c       EnD DEBUG
c-----
c-----
c       get the layer with the smallest velocity that is closest to the surface
c-----
        call setjref(lmin)
        jref = lmin
        twopi=2.*3.141592654
        ifirst = 0
c-----
c       as a safety artifice increase the halfspace vewlocity slightly
c       to that everything will work for a mdoel with all layer velocities
c       the same
c-----
        a(mmax) = 1.0001 * a(mmax)
        b(mmax) = 1.0001 * b(mmax)
  400   continue
        call gtshed(1,ifunc,mode,t,ierr)
        if(ierr.ne.0)go to 700
        s1= sngl(t)
c-----
c       DEBUG OUTPUT
c-----
c        write(6,*) ifunc,mode,s1
c-----
c       END DEBUG OUTPUT
c-----
        call puthed(2,ifunc,mode,s1)
        omega=twopi/t
        freq = 1.0/t
        if(ifunc.lt.0) go to 700
        if(mode.le.0) go to 400
        read(1) (cp(k),k=1,mode)
        do 600 k=1,mode
                c=cp(k)
c-----
c       main part.
c-----
      
               wvno = omega/c
               call legn(freq,c,uu,tt,dcdb,dcdr,dcdh,ALE,gam,
     1              ccausal,ugr,flagr,sumi0,sumi1,sumi2)
c-----
c       the gamma routine will use the spherical model, but the
c       frequency dependence and Q of the original model
c-----
            if(dogam)then
                gammal = gam
                c = ccausal
                wvno = omega/c
            else
                gammal = 0.0d+00
            endif
c-----
c       output necessary eigenfunction values for
c       source excitation
c-----
            mmaxot = mmax
C        WRITE(6,*)'c(flat)=',c,' u(flat)=',ugr
C        WRITE(6,*)'mmaxot:',mmaxot
C        WRITE(6,*)'dcdb:',dcdb(3)
                if(nsph.gt.0)then
c-----
c               sphericity correction for partial derivatives
c               of the original model
c-----

                 call splove(omega,c,mmaxot,csph,usph,ugr)
                 wvno = omega/csph
C        WRITE(6,*)'c(sph )=',csph,' u(sph )=',usph
C        WRITE(6,*)'dcdb:',dcdb(3)
                 wvno = omega/csph
                endif
c-----
c           check for underflow
c-----
            if(dabs(ugr).lt.1.0d-36)ugr=0.0d+00
            if(dabs(sumi0).lt.1.0d-36)sumi0=0.0d+00
            if(dabs(sumi1).lt.1.0d-36)sumi1=0.0d+00
            if(dabs(sumi2).lt.1.0d-36)sumi2=0.0d+00
            if(dabs(ALE).lt.1.0d-36)ale=0.0d+00
            if(dabs(flagr).lt.1.0d-36)flagr=0.0d+00
            if(dabs(gammal).lt.1.0d-36)gammal=0.0d+00

c-----
c      get the derivatives of the eigenfunctions required
c      for source excitation from the definition of stress. For
c      completeness get the second derivative from the first
c      derivatives and the equation of motion for the medium
c-----
        call getdrv(lss,lrr,omega,wvno,
     1          sut,sdut,sd2ut,suz,sduz,sd2uz,sur0,
     2          rut,rtt,ruz,rtz,wvnrec,rur0)
            sale   = sngl(ALE)
            wvnsrc = sngl(wvno)
            rale   = sngl(ALE)
            wvnrec = wvnsrc

            if( verbose ) then
                 WRITE(LOT,2)t,c,ugr,gammal,
     1                sumi0,sumi1,sumi2,
     2                flagr,ale,flagr/(sumi0*omega**2)
    2    format(' T=',e15.7,'  C=',e15.7,'  U=',e15.7,' G=',e15.7/
     1          'I0=',e15.7,' I1=',e15.7,' I2=',e15.7/
     2          ' L=',e15.7,' AL=',e15.7,' L/om2 I0',e15.7)
C                 WRITE(LOT,'(3e26.17/3e26.17/3e26.17)')t,c,ugr,gammal,
C     1                sumi0,sumi1,sumi2,
C     2                flagr,ale
            endif


            sumkr = 0.0
            sumgr = 0.0
            sumgv = 0.0
            if(nsph.gt.0)then
                wvno = omega/ csph
                ugr = usph
            endif
        if(dderiv)then
c-----
c               if a layer was inserted, get the partial
c               derivative for the original layer.
c           The sequence is cannot be changed.
c           Originally the model grows because the source layer
c           is added and then the receiver layer
c           This we must first strip the receiver and lastly the
c           source
c----
c       initialize
c-----
            if(nwlyrr)then
                call collap(lrro+1,mmaxot)
            endif
            if(nwlyrs)then
                call collap(lsso+1,mmaxot)
            endif
C            WRITE(6,*)'444'
C            do i = 1,mmax
C            WRITE(6,'(i5,6f10.6)')
C    1            I,d(i),uu(i),tt(i),dcdh(i),dcdb(i),dcdr(i)
C            enddo
            call chksiz(uu,spur,mmaxot)
            call chksiz(tt,sptr,mmaxot)
            call chksiz(dcdh,sdcdh,mmaxot)
            call chksiz(dcdb,sdcdb,mmaxot)
            call chksiz(dcdr,sdcdr,mmaxot)

            call putder(2,5,sngl(wvno),sngl(ugr), 
     1          sngl(gammal), 
     1          sut,sdut,sd2ut,suz,sduz,sd2uz,sale,wvnsrc,sur0,
     2          rut,rtt,ruz,rtz,rale,wvnrec,rur0,
     3          sumkr,sumgr,sumgv,mmaxot,
     4          sdcdh,sdcda,sdcdb,sdcdr,
     5          spur,sptr,spuz,sptz,ipar)
        else
            call putegn(2,1,1,sngl(wvno),sngl(ugr),
     1          sngl(gammal),
     1          sut,sdut,suz,sduz,sale,wvnsrc,sur0,
     2          rut,rtt,ruz,rtz,rale,wvnrec,rur0,
     3          sumkr,sumgr,sumgv)
        endif
C           WRITE(6,*)sngl(wvno),sngl(ugr)
C           WRITE(6,*)sngl(gammal)
C           WRITE(6,*)sut,sdut,sd2ut,suz,sduz,sd2uz,sale,wvnsrc,sur0
C           WRITE(6,*)rut,rtt,ruz,rtz,rale,wvnrec,rur0
C           WRITE(6,*)sumkr,sumgr,sumgv

c-----
c       DEBUG OUTPUT
c-----
c           write(6,*) wvno,c,ugr,are
c           write(6,*) uu(ls),dut
cc-----
c       END DEBUG OUTPUT
c-----
  600   continue
        go to 400
  700   continue
c-----
c       close input file from sdisp96 and output file of this program
c-----
        do 900 i=1,2
                close(i,status='keep')
  900   continue

 9999   continue
        end

        subroutine insert(nsph,dph,newlyr,ls)
        implicit none
c-----
c       procdure arguments
c-----
        integer nsph
        real*4 dph
        logical newlyr
        integer ls
c-----
c       common blocks
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qai(NL),qbi(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d,a,b,rho,qai,qbi,etap,etas,frefp,frefs
        common/modlly/mmax
        integer mmax

        integer iwat
        common/wateri/iwat(NL)

        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        real*4 vtp,dtp,rtp
c-----
c       local variables
c------
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        real dep, dp, dphh, hsave
        integer m, m1,i
        real stackh
c-----
c       stackh is todal depth to bottom of the model. If the inserted
c       depth is greater than stackh, adjust the haflspace layer
c       thickness
c-----
C       WRITE(6,*)'insert:',dph
        stackh = 0.0
        do i=1,mmax
           stackh = stackh + d(i)
        enddo
        if(dph.gt.stackh)then
            d(mmax) = (dph-stackh) + 1.0
        endif
c-----
c       Insert a depth point into the model by splitting a layer
c       so that the point appears at the top boundary of the layer
c       dph = depth of interest
c       newlyr  - L .true. layer added
c               .false. no layer added to get source eigenfunction
c-----
c-----
c       determine the layer in which the depth dph lies.
c       if necessary, adjust  layer thickness at the base
c-----
c       Here determine layer corresponding to specific depth dph
c       If the bottom layer is not thick enough, extend it
c
c       dep     - depth to bottom of layer
c       dphh    - height of specific depth above bottom of the layer
c-----
        dep = 0.0
        dp = 0.0
        dphh=-1.0
        ls = 1
        do 100 m =1, mmax
                dp = dp + d(m)
                dphh = dp - dph
                if(m.eq.mmax)then
                        if(d(mmax).le.0.0 .or. dphh.lt.0.0)then
                                d(mmax) = (dph - dp)
                        endif
                endif
                dep = dep + d(m)
                dphh = dep - dph
                ls = m
                if(dphh.ge.0.0) go to 101
  100   continue
  101   continue
c-----
c       In the current model, the depth point is in the ls layer
c       with a distance dphh to the bottom of the layer
c
c       Do not create unnecessary layers, e.g., 
c            at surface and internally
c       However do put in a zero thickness layer 
c            at the base if necessary
c-----
        if(dph .eq. 0.0)then
            newlyr = .false.
                return
        else if(dphh .eq. 0.0 .and. ls.ne.mmax)then
            ls = ls + 1
            newlyr = .false.
                return
        else
            newlyr = .true.
c-----
c               adjust layering
c-----
                 do 102 m = mmax,ls,-1
                       m1=m+1
                        d(m1)     = d(m)
                        a(m1)     = a(m)
                        b(m1)     = b(m)
                        rho(m1)   = rho(m)
                        qai(m1)   = qai(m)
                        qbi(m1)   = qbi(m)
                        frefp(m1) = frefp(m)
                        frefs(m1) = frefs(m)
                        etap(m1)  = etap(m)
                        etas(m1)  = etas(m)
                        if(nsph.gt.0)then
                            vtp(m1)  = vtp(m)
                            dtp(m1)  = dtp(m)
                            rtp(m1)  = rtp(m)
                        endif
                iwat(m1) = iwat(m)
  102           continue
                hsave=d(ls)
                d(ls) = hsave - dphh
                d(ls+1) = dphh
                mmax = mmax + 1
        endif
        return
        end

        subroutine srclyr(depth,lmax)
        implicit none
        real depth
        integer lmax
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qai(NL),qbi(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d,a,b,rho,qai,qbi,etap,etas,frefp,frefs
        common/modlly/mmax
        integer mmax

        integer i
        real dep
c-----
c       Find source/receiver boundary. It is assumed that
c       it will lie upon a boundary
c
c       lmax = source layer 
c       depth = source depth 
c-----
        dep = 0.0
        do 100 i=1,mmax
            if(abs(depth - dep).le.0.001*d(i))then
                lmax = i
                return
            endif
            dep = dep + d(i)
  100   continue
        return
        end 


        subroutine splove(om,c,mmax,csph,usph,ugr)
c-----
c       Transform spherical earth to flat earth
c       and relate the corresponding flat earth dispersion to spherical
c
c       Schwab, F. A., and L. Knopoff (1972). Fast surface wave 
c            and free
c       mode computations, in  
c            Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c            B. A. Bolt (ed),
c       Academic Press, New York
c
c       Rayleigh Wave Equations 111, 114 p 144
c
c       Partial with respect to parameter uses the relation
c       For phase velocity, for example,
c
c       dcs/dps = dcs/dpf * dpf/dps, c is phase velocity, p is
c       parameter, and f is flat model
c
c       om      R*8     angular frequency
c       c       R*8     phase velocity
c       mmax    I*4     number of layers 
c-----
        real*8 om, tm
        real*8 c, usph, csph, ugr
        parameter(NL=200)
        common/eigfun/ ut(NL),tt(NL),dcdb(NL),dcdh(NL),dcdr(NL),uu0(4)
        real*8 ut, tt, dcdb, dcdh, dcdr, uu0
        real*4 vtp,dtp,rtp
        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        a = 6371.0d0
        tm=sqrt(1.+(3.0*c/(2.*a*om))**2)
            do 20 i=1,mmax
C       WRITE(6,*)'i,dcdb,vtp,tm:',i,dcdb(i),vtp(i),tm
                dcdb(i)=dcdb(i)*  vtp(i)/(tm**3)
                dcdh(i)=dcdh(i)*  dtp(i)/(tm**3)
                dcdr(i)=dcdr(i)*  rtp(i)/(tm**3)
   20       continue
        csph = c / tm
        usph = ugr * tm
c       write(6,*)'c flat=',c,' csph=',c/tm
        return
        end

        subroutine bldsph(mmax)
c-----
c       Transform spherical earth to flat earth
c
c       Schwab, F. A., and L. Knopoff (1972). 
c            Fast surface wave and free
c       mode computations, 
c            in  Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c            B. A. Bolt (ed),
c       Academic Press, New York
c
c       Love Wave Equations  44, 45 , 41 pp 112-113
c       Rayleigh Wave Equations 102, 108, 109 pp 142, 144
c
c     Revised 28 DEC 2007 to use mid-point, assume linear variation in
c     slowness instead of using average velocity for the layer
c     Use the Biswas (1972:PAGEOPH 96, 61-74, 1972) density mapping
c
c       This is for Love Waves
c
c-----
        parameter (NL=200)
        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        real*4 vtp,dtp,rtp
        common/isomod/d(NL),ai(NL),bi(NL),rhoi(NL),
     1      qai(NL),qbi(NL),etapi(NL),etasi(NL), 
     2      frefpi(NL), frefsi(NL)
        real d,ai,bi,rhoi,qai,qbi,etapi,etasi,frefpi,frefsi
        real*4 ar, dr, r0, r1, z0, z1
c-----
c       vtp is the factor used to convert 
c            spherical velocity to flat velocity
c       rtp is the factor used to convert 
c            spherical density  to flat density 
c       dtp is the factor used to convert 
c            spherical boundary to flat boundary
c-----
c-----
c       duplicate computations of srwvds(IV)
c-----
        ar=6371.0
        dr=0.0
        r0=ar
        z0 = 0.0
        d(mmax)=1.0
        do 10 i=1,mmax
C          WRITE(6,*)'i,d,r1,r0:',i,d(i),r1,r0
            r1 = r0 * exp((-d(i))/ar)
            if(i.lt.mmax)then
                z1 = z0 + (d(i))
            else
                z1 = z0 + 1.0
            endif
            TMP=(ar+ar)/(r0+r1)
             vtp(i) = TMP
             rtp(i) = TMP**(-5)
            dtp(i) = (ar/r0)
            r0 = r1
            z0 = z1
   10   continue
C       write(6,*)'vtp:',(vtp(i),i=1,mmax)
C       write(6,*)'rtp:',(rtp(i),i=1,mmax)
C       write(6,*)'dtp:',(dtp(i),i=1,mmax)
c-----
c       at this point the model information is no longer used and
c       will be overwritten
c-----
        return
        end

        subroutine gcmdln(hsfile,hrfile,hs,hr,dotmp,dogam,dderiv,
     1      nipar,verbose,doauto)
c-----
c       parse the command line arguments
c-----
c       hsfile  C*120   - name of source depth file
c       hrfile  C*120   - name of receiver depth file
c       hs  R*4 source depth (single one specified
c       hr  R*4 receiver depth (single one specified
c       dotmp   L   - .true. use file tsdisp96.lov
c       dogam   L   - .true. incorporate Q
c       dderiv  L   - .true. output depth dependent values
c       nipar   I*4 - array o integer controls
c           set nipar(4) = 1 if eigenfunctions are output with -DER flag
c           set nipar(5) = 1 if dc/dh are output with -DER flag
c           set nipar(6) = 1 if dc/da are output with -DER flag
c           set nipar(7) = 1 if dc/db are output with -DER flag
c           set nipar(8) = 1 if dc/dr are output with -DER flag
c       verbose L   - .true. output information on energy integrals
c       doauto  L   - .true. reference layer is one with lowest velocity near surface
c                     .false. reference is the top surface (default)
c-----
        character hrfile*120, hsfile*120
        logical dotmp, dogam, dderiv
        integer nipar(10)
        logical verbose,doauto


        character name*40
        integer mnmarg

        hrfile = ' '
        hsfile = ' '
        hs = -1.0E+21
        hr = -1.0E+21
        dotmp = .false.
        dogam = .true.
        dderiv = .false.
        nipar(4) = 0
        nipar(5) = 0
        nipar(6) = 0
        nipar(7) = 0
        nipar(8) = 0
        verbose = .false.
        doauto  = .false.
        nmarg = mnmarg()
        i = 0
 1000   continue
            i = i + 1
            if(i.gt.nmarg)go to 2000
                call mgtarg(i,name)
                if(name(1:3).eq.'-HR')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')hr
                else if(name(1:3).eq.'-HS')then
                    i = i + 1
                    call mgtarg(i,name)
                    read(name,'(bn,f20.0)')hs
C               else if(name(1:4).eq.'-FHR')then
C                   i = i + 1
C                   call mgtarg(i,hrfile)
C               else if(name(1:4).eq.'-FHS')then
C                   i = i + 1
C                   call mgtarg(i,hsfile)
                else if(name(1:2) .eq. '-T')then
                    dotmp = .true.
                else if(name(1:4) .eq. '-NOQ')then
                    dogam = .false.
                else if(name(1:4) .eq. '-DER')then
                    nipar(4) = 1
                    nipar(5) = 1
                    nipar(6) = 0
                    nipar(7) = 1
                    nipar(8) = 1
                    dderiv = .true.
                else if(name(1:3).eq.'-DE')then
                    dderiv = .true.
                    nipar(4) = 1
                else if(name(1:3).eq.'-DH')then
                    dderiv = .true.
                    nipar(5) = 1
                else if(name(1:3).eq.'-DB')then
                    dderiv = .true.
                    nipar(7) = 1
                else if(name(1:3).eq.'-DR')then
                    dderiv = .true.
                    nipar(8) = 1
                else if(name(1:2).eq.'-V')then
                    verbose = .true.
               else if(name(1:5).eq.'-AUTO')then
                     doauto = .true.
               else if(name(1:4).eq.'-TOP')then
                     doauto = .false.
                else if(name(1:2) .eq. '-?')then
                    call usage(' ')
                else if(name(1:2) .eq. '-h')then
                    call usage(' ')
                endif
        go to 1000
 2000   continue
        return
        end

        subroutine usage(ostr)
        integer LER
        parameter (LER=0)
        character ostr*(*)
        if(ostr .ne. ' ')then
            lostr = lgstr(ostr)
            write(LER,*)ostr(1:lostr)
        endif
        write(LER,*)'Usage: slegn96 ',
C    1      ' -FHR recdep -FHS srcdep -HS hs -HR hr ' ,
     1      ' -HS hs -HR hr ' ,
     2       ' [-NOQ] [-T] [-DER -DE -DH -DB -DR ]',
     3       ' [-AUTO | -TOP ] [-V] [-?] [-h]'
C       write(LER,*)
C     1 '-FHS srcdep (overrides -HS )  Name of source depth  file'
C       write(LER,*)
C     1 '-FHR recdep (overrides -HR )  Name of receiver depth  file'
        write(LER,*)
     1  '-HS hs      (default 0.0 )  Source depth '
        write(LER,*)
     1  '-HR hr      (default 0.0 )  Receiver depth'
        write(LER,*)
     1  '-NOQ        (default Q used) perfectly elastic'
        write(LER,*)
     1  '-T          (default false) out tdisp96.lov not disp96.lov',
     2   ' to test TI code'
        write(LER,*)
     1  '-DER        (default false) output depth dependent values'
        write(LER,*)
     1  '-DE         (default false) output eigenfunctions(depth)'
        write(LER,*)
     1  '-DH         (default false) output DC/DH(depth)'
        write(LER,*)
     1  '-DB         (default false) output DC/DB(depth)'
        write(LER,*)
     1  '-DR         (default false) output DC/DR(depth)'
        write(LER,*)
     1  '-V          (default false) list energy integrals'
        write(LER,*)
     1  '-TOP        (default true) Top boundary is reference'
        write(LER,*)
     1  '-AUTO       (default false) Lowest velocity near surface',
     2             ' is reference'
        write(LER,*)
     1  '-?          (default none )  this help message '
        write(LER,*)
     1  '-h          (default none )  this help message '
        stop
        end

        subroutine collap(ls,mmaxot)
        parameter(NL=200)
        common/eigfun/ uu(NL),tt(NL),dcdb(NL),dcdh(NL),dcdr(NL),uu0(4)
        real*8 uu, tt, dcdb, dcdh, dcdr, uu0
        do 501 i = ls-1,mmaxot
            if(i .eq. ls -1)then
                dcdb(i) = dcdb(i) + dcdb(i+1)
                dcdr(i) = dcdr(i) + dcdr(i+1)
            endif
            if(i.gt.ls)then
                dcdb(i-1) = dcdb(i)
                dcdh(i-1) = dcdh(i)
                dcdr(i-1) = dcdr(i)
                uu(i-1) = uu(i)
                tt(i-1) = tt(i)
            endif
  501   continue
        mmaxot = mmaxot - 1
        return
        end

        subroutine chksiz(dp,sp,mmaxot)
c-----
c       correctly convert double precision to single precision
c-----
        real*8 dp(mmaxot)
        real*4 sp(mmaxot)
            do 610 i=1,mmaxot
                if(dabs(dp(i)).lt.1.0d-36)then
                    sp(i) = 0.0
                else
                    sp(i) = sngl(dp(i))
                endif
  610       continue
        return
        end

        subroutine getdrv(lss,lrr,omega,wvno,
     1          sut,sdut,sd2ut,suz,sduz,sd2uz,sur0,
     2          rut,rtt,ruz,rtz,wvnrec,rur0)
        implicit none
c-----
c       get the eigenfunctions and derivatives at the source depth
c       get the eigenfunctions at the receiver depth
c-----
        integer lss, lrr
        real*8 omega, wvno
        real sut,sdut,sd2ut,suz,sduz,sd2uz,wvnsrc,sur0
        real rut,rtt,ruz,rtz,wvnrec,rur0

        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real*4 d, a, b, rho, qa, qb, etap, etas, frefp, frefs

        common/modlly/mmax
        integer mmax

        common/wateri/iwat(NL)
        integer iwat

        common/eigfun/ uu(NL),tt(NL),dcdb(NL),dcdh(NL),dcdr(NL),uu0(4)
        real*8 uu, tt, dcdb, dcdh, dcdr, uu0

c-----
c       local variables
c-----
        real xmu
        real*8 wvno2
        real*8 Eut, Edut, Ed2ut, Eut0, Ett0

        wvno2 = wvno*wvno
c-----  
c       compute the eigenfuntions and depth derivatives
c       at the source depth
c-----  
        if(iwat(lss).eq.0)then
            Eut = uu(lss)
            xmu = rho(lss)*b(lss)*b(lss)
            Edut = tt(lss)/xmu
            Ed2ut = ( - (omega/b(lss))**2 + wvno2)*Eut
        else 
            Eut = 0.0d+00
            Edut = 0.0d+00
            Ed2ut = 0.0d+00
        endif
        if(iwat(lrr).eq.0)then
            Eut0 = uu(lrr)
            Ett0 = tt(lrr) 
        else
            Eut0 = 0.0d+00
            Ett0 = 0.0d+00
        endif

            if(dabs(Eut).lt. 1.0d-36)Eut = 0.0d+00
            if(dabs(Edut).lt. 1.0d-36)Edut = 0.0d+00
            if(dabs(Ed2ut).lt. 1.0d-36)Ed2ut = 0.0d+00
            if(dabs(Eut0).lt. 1.0d-36)Eut0 = 0.0d+00
            if(dabs(Ett0).lt. 1.0d-36)Ett0 = 0.0d+00

c-----
c      get the derivatives of the eigenfunctions required
c      for source excitation from the definition of stress. For
c      completeness get the second derivative from the first
c      derivatives and the equation of motion for the medium
c-----
            sut = sngl(Eut)
            sdut = sngl(Edut)
            sd2ut = sngl(Ed2ut)
            suz = 0.0
            sduz = 0.0
            sd2uz = 0.0
            wvnsrc = sngl(wvno)
            sur0 = 0.0

            rut = sngl(Eut0)
            rtt = sngl(Ett0)
            ruz = 0.0
            rtz = 0.0
            wvnrec = sngl(wvno)
            rur0 = 0.0
        return
        end
 
        subroutine setjref(lmin)
c-----
c     determine the layer with lowest velocity
c     nearest the surface. then set the jref value for use 
c     with rsh.f or rpsv.f to perform the recursion
c-----
      implicit none
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etapi(NL),etasi(NL),
     2      frefpi(NL), frefsi(NL)
        real d,a,b,rho,qa,qb,etapi,etasi,frefpi,frefsi

        common/modlly/mmax
        integer mmax
        integer lmin

        integer  i
        real vmin

        if(b(mmax).gt.0.0)then
            vmin = b(mmax)
        else
            vmin = a(mmax)
        endif
        lmin = mmax

        do i=mmax-1,1,-1
           if(b(i).eq.0.0)then
                if(a(i).le.vmin)then
                     vmin=a(i)
                     lmin = i
                endif
           else
                if(b(i).le.vmin)then
                     vmin=b(i)
                     lmin = i
                endif
           endif
        enddo
        return
        end
       
