        program sregn96
c----------------------------------------------------------------------c
c                                                                      c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                 c
c      VOLUME III                                                      c
c                                                                      c
c      PROGRAM: SREGN96                                                c
c                                                                      c
c      COPYRIGHT 1996, 2010                                            c
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
c     Rewrite of theory to agree more with hspec96 wavenumber
c     integration code
c
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c Revision history:
c       07 AUG 2002 - make string lengths 120 characters from 80 
c       13 OCT 2006 - verbose output of energy integrals
c       26 SEP 2008 - fixed undefined LOT in subroutine up
c       14 JUN 2009 - reformulate in general terms as per book
c       01 AUG 2010 - change code to agree with book  for
c                     migration to TI. Also clean up code using
c                     implicit none
c       20 SEP 2012 - cleaned up some compiler warnings
c       25 FEB 2017 - fixed error in subroutine varsv. The
c           computation of sin (nu z)/nu yielded NaN when nu = 0
c           The corrected code gives lim  sin(nu z) / nu = z
c                                         nu->0
c       29 JUL 2017 - modify insert if depth is in halfspace
c       22 FEB 2018 - verbose also gives LAGR/(omega^2 I0)
c     19 MAR 2022 - corrected some numebrical problems
c           This arose when using a model in cm, cm/s and gm/cm^3 
c           instead of km, km,s and km/cm^3.  Numericla precision was
c           lost. Statements such as
c           if(dabs(pr) .lt. 1.0e-5 .and. cdabs(rp).lt.1.0e-5)then
c           were changes to use a much lower limit since double precision
c           permitted this.
c           Note that was only a problem with CGS units and not MKS or
c           the km, km/s gm/cm^3 mixed units
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        implicit none
c-----
c       common blocks
c-----
        integer NL
        parameter(NL=200)
c-----
c       LIN - unit for FORTRAN read from terminal
c       LOT - unit for FORTRAN write to terminal
c       LER - unit for FORTRAN error output to terminal
c       NL  - number of layers in model
c-----
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qai(NL),qbi(NL),etapi(NL),etasi(NL),
     2      frefpi(NL), frefsi(NL)
        real d,a,b,rho,qai,qbi,etapi,etasi,frefpi,frefsi

        common/modlly/mmax
        integer mmax

        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0, dcda, dcdb, dcdr, dcdh
        common/depref/refdep
        real refdep



        real*4 sdcda(NL), sdcdb(NL), sdcdh(NL), sdcdr(NL) 
        real*4 spur(NL), sptr(NL), spuz(NL), sptz(NL)

        common/wateri/iwat(NL)
        integer iwat
c-----
c       model characterization
c-----
        common/modspec/allfluid
        logical allfluid

        common/sumi/ sumi0,sumi1,sumi2,sumi3,flagr,are,ugr
        real*8 sumi0, sumi1, sumi2, sumi3, flagr, are, ugr

        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        real*4 vtp,dtp,rtp

c-----
c       function prototype
c-----
        integer lgstr

c-----
c       local variables
c-----
        integer LER, LIN, LOT
        parameter(LER=0,LIN=5,LOT=6)
        integer ipar(10)
        integer nipar(10)
        integer i,k,l
        integer iunit,iiso,iflsph,idimen,icnvel,ierr
        integer MAXMOD
        parameter(MAXMOD=2000)
        integer lss, lsso, lrr, lrro, n1, n2, npts, ifirst
        integer ifunc, mode, nsph, nper, nmodes, mmaxot

        logical dotmp
        logical ext
        logical dogam
        logical dderiv

        common/verby/verbose, dout, doauto
        logical verbose, dout, doauto
        real*8 ccausal, gam

        common/control/jref
        integer jref
        integer lmin

        logical nwlyrs, nwlyrr
        logical dolove, dorayl

        character mname*120
        character hsfile*120, hrfile*120, title*120 
        character*12 fname(2)

        real*4 fpar(10)
        real*4 s1, dephs, dephr, faclov, facray, hr, hs, dt
        real*4 deplw, depthr, depths, depup
        real*4 ohs, ohr
        real*4 rare, rtr, rtz, rur0, rur, ruz
        real*4 sare, sd2ur, sd2uz, sduz, sur, sur0, suz, sdur
        real*4 sumgr, sumgv, sumkr
        real*4 wvnsrc, wvnrec
        real*4 twopi

        real*8 c, freq, omega, wvno, gammar, csph, usph
        real*8 cp(MAXMOD), t


c-----
c       machine dependent initialization
c-----
        call mchdep()
c-----
c       parse command line information
c-----  
        call gcmdln(hsfile,hrfile,hs,hr,dotmp,dogam,dderiv,
     1      nipar,verbose, doauto)
c-----
c       wired in file names     - output of this program 
c                           is always in the
c               binary file sregn96.egn or sregn96.der
c                               - input from sdisp96 is always
c               binary file sdisp96.ray
c-----
        if(dotmp)then
            fname(1) = 'tsdisp96.ray'
        else
            fname(1) = 'sdisp96.ray'
        endif
        if(dderiv)then
            fname(2) = 'sregn96.der'
        else
            fname(2) = 'sregn96.egn'
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
c-----
c       check for water only
c-----
        allfluid = .true.
        do 1200 i=1,mmax
            if(b(i).gt.0.0)then
                allfluid = .false.
            endif
 1200   continue
        nsph = iflsph
        if(nsph.gt.0)then
                call bldsph()
        endif
c-----
c       get the Q information into the program
c       since this is not carried through in the sdisp96 output
c-----
        do 1234 i=1,mmax
            if(.not.dogam)then
                qai(i) = 0.0
                qbi(i) = 0.0
            endif
            if(qai(i).gt.1.0)qai(i)=1.0/qai(i)
            if(qbi(i).gt.1.0)qbi(i)=1.0/qbi(i)
            if(frefpi(i).le.0.0)frefpi(i) = 1.0
            if(frefsi(i).le.0.0)frefsi(i) = 1.0
 1234   continue
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
c       see if the file sdisp96.ray exists
c-----
        inquire(file=fname(1),exist=ext)
        if(.not.ext)then
            call usage('Dispersion file: '//fname(1)//
     1          ' does not exist')
        endif
c-----
c       open output file of sdisp96 and of sregn96
c-----
        open(1,file=fname(1),form='unformatted',status='unknown',
     1          access='sequential')
        open(2,file=fname(2),form='unformatted',status='unknown',
     1          access='sequential')
        rewind 1
        rewind 2
c-----
c       obtain the earth model: note if the original was spherical,
c       this will be the transformed model that must be used to
c       get the eigenfunctions
c-----
        call gtsmdl(1,mmax,d,a,b,rho,qai,qbi,nper,
     2              mname,ipar,fpar)
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
        call putmdl(2,mmax,d,a,b,rho,qai,qbi,nper,depths-refdep,
     2      depthr-refdep,mname,ipar,fpar)
c-----
c               Insert source and receiver depths into the model used
c               for the dispersion. Note that if the model is spherical
c               we have already mapped the spehrical depth into
c               the flattened model.
c-----
c-----
c       split a layer at the source depth
c       the integer ls gives the interface at which the eigenfunctions
c       are to be defined
c-----
        call insert(nsph,dephs,nwlyrs,lsso)
        call insert(nsph,dephr,nwlyrr,lrro)
        call srclyr(dephs,lss)
        call srclyr(dephr,lrr)
c-----DEBUG
c       output the new model with source layer
c-----
C       write(6,*)'lss,lrr:',lss,lrr
C        do 9902 i=1,mmax
C                write(6,*)iwat(i),d(i),a(i),b(i),rho(i)
C 9902   continue
c-----
c       END DEBUG
c-----
c-----
c       get the layer with the smallest velocity that is closest to the surface
c-----
        call setjref(lmin)
c-----
c       set the jref
c-----
        if(lmin.eq.1)then
              jref = 0
        else
              jref = lmin
        endif
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
        s1=sngl(t)
c-----
c       DEBUG OUTPUT
c-----
c        write(6,*) ifunc,mode,s1
c-----
c       END DEBUG OUTPUT
c-----
        call puthed(2,ifunc,mode,s1)
        omega=twopi/t
        freq = 1./t
        if(ifunc.lt.0) go to 700
        if(mode.le.0) go to 400
        read(1) (cp(k),k=1,mode)
        do 600 k=1,mode
                c=cp(k)
c-----
c       main part.
c-----
            wvno=omega/c
            call regn(freq,c,UZ,UR,TZ,TR,dcda,dcdb,dcdr,dcdh,
     1            ARE,gam,ccausal,ugr,
     1            flagr,sumi0,sumi1,sumi2,sumi3)
c-----  
c       the gamma routine will use the spherical model, but the
c       frequency dependence and Q of the original model
c-----
            if(dogam)then
                gammar = gam
                c = ccausal
                wvno = omega/c
            else
                gammar = 0.0d+00
            endif
            rur0 = sngl(UR(1))
            sur0 = sngl(UR(1))
c------
c     also check for possible conversion errors in IEEE
c     conversion from double precision to single precision
c-----
            if(dabs(uu0(1)).lt.1.0d-36)uu0(1)=0.0d+00
            if(dabs(uu0(3)).lt.1.0d-36)uu0(3)=0.0d+00
            if(dabs(c).lt.1.0d-36)c=0.0d+00
            if(dabs(ugr).lt.1.0d-36)ugr=0.0d+00
            if(dabs(are).lt.1.0d-36)are=0.0d+00
            if(dabs(flagr).lt.1.0d-36)flagr=0.0d+00
            if(dabs(gammar).lt.1.0d-36)gammar=0.0d+00

c-----
c      get the derivatives of the eigenfunctions required
c      for source excitation from the definition of stress. For
c      completeness get the second derivative from the first
c      derivatives and the equation of motion for the medium
c-----
        call  getdrv(lss,lrr,omega,wvno,
     1           sur,sdur,sd2ur,suz,sduz,sd2uz,sur0,
     2            rur,rtr,ruz,rtz,wvnrec,rur0)
            sare   = sngl(ARE)
            wvnsrc = sngl(wvno)
            rare   = sngl(ARE)
            wvnrec = wvnsrc


c-----
c       output necessary eigenfunction values for
c       source excitation
c-----
            mmaxot = mmax
            if(nsph.gt.0)then
c-----
c           sphericity correction for partial derivatives 
c           of the original model
c-----
                call sprayl(omega,c,mmaxot,csph,usph,ugr)
                wvno = omega / csph

            endif

            if( verbose ) then
                 WRITE(LOT,2)t,c,ugr,gammar,
     1                sumi0,sumi1,sumi2,sumi3,
     2                flagr,are,sur0,flagr/(sumi0*omega**2)
    2    format(' T=',e15.7,'  C=',e15.7,'  U=',e15.7,'  G=',e15.7/
     1          'I0=',e15.7,' I1=',e15.7,' I2=',e15.7,' I3=',e15.7/
     2          ' L=',e15.7,' AR=',e15.7,'  E=',f15.7,' L/om2 I0',e15.7)
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
c           If a layer was inserted, get the partial
c           derivative for the original layer.
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
            call chksiz(ur,spur,mmaxot)
            call chksiz(tr,sptr,mmaxot)
            call chksiz(uz,spuz,mmaxot)
            call chksiz(tz,sptz,mmaxot)
            call chksiz(dcdh,sdcdh,mmaxot)
            call chksiz(dcdb,sdcdb,mmaxot)
            call chksiz(dcda,sdcda,mmaxot)
            call chksiz(dcdr,sdcdr,mmaxot)

            call putder(2,6,sngl(wvno),sngl(ugr), 
     1          sngl(gammar), 
     1          sur,sdur,sd2ur,suz,sduz,sd2uz,sare,wvnsrc,sur0,
     2          rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     3          sumkr,sumgr,sumgv,mmaxot,
     4          sdcdh,sdcda,sdcdb,sdcdr,
     5          spur,sptr,spuz,sptz,ipar)
            if(verbose)then
                WRITE(LOT,3)sur,sdur,sd2ur,suz,sduz,sd2uz,sare
                WRITE(LOT,4)rur,rtr,ruz,rtz,rare
    3   format(' SUR=',e15.7,' SDUR=',e15.7,' SD2UR=',e15.7,
     1    ' SUZ=',e15.7,' SDUZ=',e15.7,' SD2DZ=',e15.7,' SARE=',e15.7)
    4   format(' RUR=',e15.7,' RTR=',e15.7,' RUZ=',e15.7,
     1    ' RTZ=',e15.7,' RARE=',e15.7)

 
            endif
        else
            call putegn(2,2,1,sngl(wvno),sngl(ugr),
     1          sngl(gammar),
     1          sur,sdur,suz,sduz,sare,wvnsrc,sur0,
     2          rur,rtr,ruz,rtz,rare,wvnrec,rur0,
     3          sumkr,sumgr,sumgv)
        endif
C           WRITE(6,*)sngl(wvno),sngl(ugr)
C           WRITE(6,*)sngl(gammar)
C           WRITE(6,*)sur,sdur,sd2ur,suz,sduz,sd2uz,sare,wvnsrc,sur0
C           WRITE(6,*)rur,rtr,ruz,rtz,rare,wvnrec,rur0
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
                        if(d(mmax).le.0.0d+00 .or. dphh.lt.0.0)then
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

        subroutine sprayl(om,c,mmax,csph,usph,ugr)
c-----
c       Transform spherical earth to flat earth
c       and relate the corresponding flat earth dispersion to spherical
c
c       Schwab, F. A., and L. Knopoff (1972). 
c          Fast surface wave and free
c       mode computations, in  
c          Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c          B. A. Bolt (ed),
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
        implicit none
c-----
c       procedure arguments
c-----
        real*8 om, c, csph,usph,ugr
        integer mmax
c-----
c       common blocks
c-----
        integer NL
        parameter(NL=200)
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0
        real*8 dcda, dcdb, dcdr, dcdh
        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        real*4 vtp,dtp,rtp
c-----
c       local arguments
c-----
        real*8 ar, tm
        integer i

        ar = 6370.0d0
        tm=sqrt(1.+(c/(2.*ar*om))**2)
            do 20 i=1,mmax
C        WRITE(6,*)'i,dcdb,vtp,tm:',i,dcdb(i),vtp(i),tm
                dcda(i)=dcda(i)*  vtp(i)/(tm**3)
                dcdb(i)=dcdb(i)*  vtp(i)/(tm**3)
                dcdh(i)=dcdh(i)*  dtp(i)/(tm**3)
                dcdr(i)=dcdr(i)*  rtp(i)/(tm**3)
   20       continue
c       write(6,*)'c flat=',c,' csph=',c/tm
        usph = ugr*tm
        csph = c/tm
        return
        end

        subroutine bldsph()
c-----
c       Transform spherical earth to flat earth
c
c       Schwab, F. A., and L. Knopoff (1972). 
c          Fast surface wave and free
c       mode computations, in  
c          Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c          B. A. Bolt (ed),
c       Academic Press, New York
c
c       Love Wave Equations  44, 45 , 41 pp 112-113
c       Rayleigh Wave Equations 102, 108, 109 pp 142, 144
c
c     Revised 28 DEC 2007 to use mid-point, assume linear variation in
c     slowness instead of using average velocity for the layer
c     Use the Biswas (1972:PAGEOPH 96, 61-74, 1972) density mapping
c
c       This is for Rayleigh Waves
c
c-----
        implicit none
        integer NL
        parameter (NL=200)
c-----
c       common blocks
c-----
        common/sphere/vtp(NL),dtp(NL),rtp(NL)
        real*4 vtp,dtp,rtp

        common/isomod/d(NL),a(NL),b(NL),rho(NL),zqai(NL),zqbi(NL),
     1      zetap(NL),zetas(NL),zfrefp(NL),zfrefs(NL)
        real d, a, b, rho, zqai, zqbi, zetap, zetas, 
     1      zfrefp, zfrefs
        common/modlly/mmax
        integer mmax
c-----
c       local arguments
c-----
        real ar, dr, r0, r1, z0, z1
        real tmp
        integer i
c-----
c       vtp is the factor used to convert spherical 
c          velocity to flat velocity
c       rtp is the factor used to convert spherical 
c          density  to flat density 
c       dtp is the factor used to convert spherical 
c          boundary to flat boundary
c-----
c-----
c       duplicate computations of srwvds(IV)
c-----
        ar=6370.0d0
        dr=0.0d0
        r0=ar
        z0 = 0.0d+00
        d(mmax)=1.0
        do 10 i=1,mmax
            r1 = r0 * exp(-d(i)/ar)
            if(i.lt.mmax)then
                z1 = z0 + d(i)
            else
                z1 = z0 + 1.0
            endif
            TMP=(ar+ar)/(r0+r1)
            vtp(i) = TMP
            rtp(i) = TMP**(-2.275)
            dtp(i) = (ar/r0)
            r0 = r1
            z0 = z1
   10   continue
C        write(6,*)'vtp:',(vtp(i),i=1,mmax)
C        write(6,*)'rtp:',(rtp(i),i=1,mmax)
C        write(6,*)'dtp:',(dtp(i),i=1,mmax)
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
        implicit none
c-----
c       procedure arguments
c-----
        character hrfile*120, hsfile*120
        real hs, hr
        logical dotmp, dogam, dderiv
        integer nipar(10)
        logical verbose,doauto

c-----
c       procedure prototypes
c-----
        integer mnmarg
c-----
c       local arguments
c-----
        character name*40
        integer i, nmarg

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
        doauto = .false.
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
                    nipar(6) = 1
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
        implicit none
c-----
c       procedure arguments
c-----
        character ostr*(*)
c-----
c       local variables
c-----
        integer LER
        parameter (LER=0)
        integer lostr
c-----
c       function prototype
c-----
        integer lgstr

        if(ostr .ne. ' ')then
            lostr = lgstr(ostr)
            write(LER,*)ostr(1:lostr)
        endif
        write(LER,*)'Usage: sregn96 ',
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
        implicit none
c-----
c       routine arguments
c-----
        integer ls, mmaxot
c-----
c       common blocks
c-----
        integer NL
        parameter(NL=200)
        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),uu0(4),
     *                 dcda(NL),dcdb(NL),dcdr(NL),dcdh(NL)
        real*8 ur, uz, tz, tr, uu0, dcda, dcdb, dcdr, dcdh
c-----
c       local arguments
c-----
        integer i

        do 501 i = ls-1,mmaxot
            if(i .eq. ls -1)then
                dcda(i) = dcda(i) + dcda(i+1)
                dcdb(i) = dcdb(i) + dcdb(i+1)
                dcdr(i) = dcdr(i) + dcdr(i+1)
            endif
            if(i.gt.ls)then
                dcda(i-1) = dcda(i)
                dcdb(i-1) = dcdb(i)
                dcdh(i-1) = dcdh(i)
                dcdr(i-1) = dcdr(i)
                ur(i-1) = ur(i)
                tr(i-1) = tr(i)
                uz(i-1) = uz(i)
                tz(i-1) = tz(i)
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
     1          sur,sdur,sd2ur,suz,sduz,sd2uz,sur0,
     2          rur,rtr,ruz,rtz,wvnrec,rur0)
        implicit none
c-----
c       get the eigenfunctions and derivatives at the source depth
c       get the eigenfunctions at the receiver depth
c-----
        integer lss, lrr
        real*8 omega, wvno
        real sur,sdur,sd2ur,suz,sduz,sd2uz,wvnsrc,sur0
        real rur,rtr,ruz,rtz,wvnrec,rur0

        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs

        common/modlly/mmax
        integer mmax

        common/wateri/iwat(NL)
        integer iwat

        common/eigfun/ ur(NL),uz(NL),tz(NL),tr(NL),dcda(NL),
     1     dcdb(NL),dcdh(NL),dcdr(NL),uu0(4)
        real*8 ur,uz,tz,tr,dcda, dcdb, dcdh, dcdr, uu0

        common/verby/verbose, dout, doauto
        logical verbose, dout, doauto

c-----
c       local variables
c-----
        real xmu,xl2m,xlam
        real*8 wvno2
        real*8 duzdz, durdz, d2urdz, d2uzdz
        real*8 ruzdz, rurdz

        wvno2 = wvno*wvno

c-----
c      get the derivatives of the eitenfunctions required
c      for source excitation from the definition of stress. For
c      completeness get the second derivative from the first
c      derivatives and the equation of motion for the medium
c-----
c-----
c      parameters for source
c----
            xl2m = rho(lss)*a(lss)*a(lss)
            xmu  = rho(lss)*b(lss)*b(lss)
            xlam = xl2m - 2* xmu
            duzdz = ( tz(lss) + wvno*xlam*ur(lss))/xl2m
            if(iwat(lss).eq.1)then
                durdz = wvno*uz(lss)
                d2urdz = wvno*duzdz
            else
                durdz = -wvno*uz(lss) + tr(lss)/xmu
                d2urdz = (-wvno*duzdz*(xlam+xmu) 
     1              - ur(lss)*(rho(lss)*omega*omega -
     2              wvno*wvno*xl2m))/xmu
            endif
            d2uzdz = ( - uz(lss) *( rho(lss)*omega*omega - 
     1          xmu*wvno*wvno) + 
     2          wvno*durdz*(xlam+xmu))/ xl2m

c-----
c     parameters for receiver
c-----
            xl2m = rho(lrr)*a(lrr)*a(lrr)
            xmu  = rho(lrr)*b(lrr)*b(lrr)
            xlam = xl2m - 2* xmu
            ruzdz = ( tz(lrr) + wvno*xlam*ur(lrr))/xl2m
            if(iwat(lrr).eq.1)then
                rurdz = wvno*uz(lrr)
            else
                rurdz = -wvno*uz(lrr) + tr(lrr)/xmu
            endif
            if(dabs(duzdz).lt.1.0d-36)duzdz=0.0d+00
            if(dabs(durdz).lt.1.0d-36)durdz=0.0d+00
            if(dabs(d2uzdz).lt.1.0d-36)d2uzdz=0.0d+00
            if(dabs(d2urdz).lt.1.0d-36)d2urdz=0.0d+00
            if(dabs(ruzdz).lt.1.0d-36)ruzdz=0.0d+00
            if(dabs(rurdz).lt.1.0d-36)rurdz=0.0d+00


            sur = sngl(ur(lss))
            sdur = sngl(durdz)
            sd2ur = sngl(d2urdz)
            suz = sngl(uz(lss))
            sduz = sngl(duzdz)
            sd2uz = sngl(d2uzdz)
            wvnsrc = sngl(wvno)
            sur0 = sngl(ur(1))

            rur = sngl(ur(lrr))
            rtr = sngl(tr(lrr))
            ruz = sngl(uz(lrr))
            rtz = sngl(tz(lrr))
            wvnrec = sngl(wvno)
            rur0 = sngl(ur(1))

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
       
