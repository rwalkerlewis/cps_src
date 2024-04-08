      subroutine legn(f,c,UT,TT,dcdb,dcdr,dcdh,AL,gam,ccausal,u,
     1   flagr,sumi0,sumi1,sumi2)
      implicit none
      real*8 f,c
      
      integer NL
      parameter (NL=200)
      real*8 UT(NL), TT(NL)
      real*8 dcdb(NL), dcdr(NL), dcdh(NL)
      real*8 AL
      real*8  gam
      real*8  ccausal
      real*8 u
      real*8 flagr,sumi0, sumi1,sumi2



      common/downsh/RD,TD
      complex*16 RD(NL), TD(NL)

      common/upsh/RU,TU
      complex*16 RU(NL), TU(NL)

      common/coeffinalsh/CDN,CUP
      complex*16 CUP(NL), CDN(NL)

      common /elmat/el,elinv,es,nub
      complex*16 el(NL,2,2)
      complex*16 elinv(NL,2,2)
      complex*16 es(NL),nub(NL)

      common/modlly/mmax
      integer mmax
 
        common/shcontrol/ltop,lbot
        integer ltop, lbot

        common/control/jref
        integer jref

        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep

c-----
c     el = E layer
c     elinv = E^-1 layer
c     es = exp(-nu d) layer i
c-----


c-----
c     internal variables
c-----
      common/bcsh/alpsh,betsh
      complex*16 alpsh(2,2), betsh(2,2)
      integer jtop, jbot
      
      complex*16 om, wvno
      integer j
      integer i, jlyr, nbc

      complex*16 TUT(NL),TTT(NL)
      complex*16 TUT0
      complex*16 Iij,IUU, IdUdU

      complex*16 cdcdb(NL),cdcdr(NL)
      complex*16 csumi0, csumi1,csumi2
      real*8 mu
      complex*16 ash(2,2)
      real x
      real*8  dc, pi

c-----
c      test
c-----
CRBH     WRITE(0,*)'f=',f,' c=',c,' mode=',mode
CRBH     WRITE(6,*)'mmax=',mmax,' a=',a(mmax),' b=',b(mmax),'r=',rho(mmax)


c-----

      ltop = 1
      lbot = mmax

      om = dcmplx(6.2831853d+00*f,0.0d+00)
      wvno = om/c
c-----
c     for SH neither the source nor the receiver can be
c     in a fluid. All layers between the source and the
c     receiver must be solid
c-----
c-----
c         evaluate the E Einv matrices, exponentials
c-----
       call doel(om,wvno)
c-----
c         define the boundary matrices
c         jbot = 3 halfspace
c         jtop = 1 free surface
c-----
       jbot = 3
       call dobotsh(alpsh, jbot)
       jtop = 1
       call dotopsh(betsh, jtop)
c-----
c         evaluate the reflection coeffcients
c-----
       call botupsh()
       call topdownsh()
c-----
c         now do the recursion to get the Cu Cd and then the solution
c          initialize
c
c-----
       if(jref.eq.0)then
            jlyr = 1
       else if(jref.ge.mmax)then
            jlyr =mmax -1
       else
            jlyr = jref
       endif

       call cinitsh(CUP(jlyr))
       CDN(jlyr) = RU(jlyr)*CUP(jlyr)
       do j=jlyr-1,1,-1
          CUP(j) = TU(j+1)*CUP(j+1)
          CDN(j) = RU(j)  *CUP(j)
       enddo
       do j=jlyr+1,mmax-1
          CDN(j) = TD(j-1)*CDN(j-1)
          CUP(j) = RD(j)  *CDN(j)
       enddo
           CDN(mmax) = TD(mmax-1)*CDN(mmax-1)
           CUP(mmax) = dcmplx(0.0d+00,0.0d+00)
c-----
c         now get the eigenfunctions
c-----
C   2     FORMAT(i5,2(1x,2e10.3,1x))
C   3     FORMAT('*',i4,2(1x,2e10.3,1x))
c-----
c      EIGEN FUNCTIONS
c
c      top of layer
c      B = E11 enu CU + E12 CD
c          E21 enu CU + E22 CD
c      bottom of layer
c      B = E11 Cu + E12 enu CD
c          E21 Cu + E22 enu CD
c      here we compute at the top of each layer.
c      for the halfspace CUP = 0 and es(mamx) = =
c      so use the same formula rewrite so that es(mmax)
c      is not used
       do j=1,mmax
               TUT(j) = el(j,1,1)*es(j)*CUP(j)+el(j,1,2)*CDN(j)
               TTT(j) = el(j,2,1)*es(j)*CUP(j)+el(j,2,2)*CDN(j)
            if(j.eq.1)then
                 TUT0 = tut(j)
            endif
        enddo
        do j=1,mmax
            UT(j) = sngl(dreal(TUT(j)/TUT0))
            TT(j) = sngl(dreal(TTT(j)/TUT0))
        enddo
c-----
c         get the energy integrals
c-----
          csumi0 = dcmplx(0.0d+00, 0.0d+00)
          csumi1 = dcmplx(0.0d+00, 0.0d+00)
          csumi2 = dcmplx(0.0d+00, 0.0d+00)
c-----
c         layers
c-----
       do i=1,mmax
           if(i.lt.mmax)then
              nbc = 1
           else
              nbc = 2
           endif
            
           call amatsh(om,wvno,ash,i)
           mu = rho(i)*b(i)*b(i)
c-----
c          IUU = int U^2 dz
           call shIijlay(i,1,1,CUP(i),CDN(i),Iij,nbc)
c-----
c          IUU = int U^2 dz
c-----
           IUU = Iij
c-----
c          IdUdU = int (dU/dz)^2 dz
c-----
           call shIijlay(i,2,2,CUP(i),CDN(i),Iij,nbc)
           IdUdU = Iij*ash(1,2)*ash(1,2)
           csumi0 = csumi0 + rho(i)*IUU
           csumi1 = csumi1 + mu    *IUU
           csumi2 = csumi2 + mu * IdUdU
           cdcdb(i) = rho(i)*b(i)*(IUU + IdUdU/(wvno**2) )
           cdcdr(i) = (0.5*b(i)/rho(i))*cdcdb(i) - 0.5*c*c*IUU
       enddo
c-----
c         get group velocity
c-----
        sumi0 = dreal(csumi0/TUT0**2)
        sumi1 = dreal(csumi1/TUT0**2)
        sumi2 = dreal(csumi2/TUT0**2)
CRBH        WRITE(6,*)'I0:',sumi0
CRBH        WRITE(6,*)'I1:',sumi1
CRBH        WRITE(6,*)'I2:',sumi2
        u =  sumi1/(c*sumi0)
        AL=1.0/(2*c*u*sumi0)
        flagr=dreal(om*om)*sumi0
     1      - dreal(wvno*wvno)*sumi1
     2      - sumi2
c-----
c       get the dcdh
c-----
        call ddh(c,f,sumi1,UT,TT,dcdh)
c-----
c       get final partials and gamma
c-----
        gam = 0.0
        dc = 0.0
        pi = 3.141592653589493d+00
        do i=1,mmax
           dcdb(i) = dreal(cdcdb(i)/csumi0)/u
           dcdr(i) = dreal(cdcdr(i)/csumi0)/u
           
           if(qb(i).gt.1.0)then
              x = sngl(dcdb(i)) * b(i)/qb(i)
           else 
              x = sngl(dcdb(i)) * b(i)*qb(i)
           endif
           gam = gam + x
           dc = dc + dlog(f/frefs(i))*x/pi
        enddo
        gam=0.5*6.2831853*f*gam/(c*c)
        ccausal = c + sngl(dc)
CRBH        write(6,21)mode,1./f,ccausal,u,AL,gam, refdep
CRBH   21   format(//17x,'  LOVE WAVE        MODE #',i3/
CRBH     1  '        T =',e11.4,' C =   ',e11.4,' U   =',e11.4/
CRBH     2  '        AL=',e11.4,' GAMMA=',e11.4,' ZREF=',e11.4/
CRBH     3  '    M       UT         TT       ',
CRBH     4  'DC/DH      DC/DB      DC/DR')
CRBH       do i=1,mmax
CRBH         write(6,1)i,UT(i),TT(i),dcdh(i),dcdb(i),dcdr(i)
CRBH       enddo
CRBH    1   format(i5,1x,5e11.3)
       
       return
       end
      
   
       subroutine cinitsh(cval)
c-----
c         initialize the C value for the recutsion
c         the SH case is very simple
c-----
       implicit none
       complex*16 cval
       cval = dcmplx(1.0d+00, 0.0d+00)
       return 
       end

       subroutine amatsh(omega,wvno,ash,lyr)
       implicit none
c-----
c         get the SH A matrix of
C         dB/dz = A B
c-----

       complex*16 omega, wvno
       complex*16 ash(2,2)
       integer lyr


      integer NL
      parameter (NL=200)
       common/modlly/mmax
       integer mmax
 
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep

c-----
c       internal variables
c-----
        complex*16 cdzero
        real mu

        mu = rho(lyr)*b(lyr)*b(lyr)

        cdzero = dcmplx(0.0d+00,0.0d+00)
        ash(1,1) = cdzero
        ash(2,2) = cdzero
        ash(1,2) = dcmplx(1.0d+00,0.0d+00)/mu
        ash(2,1) = mu*wvno*wvno - rho(lyr)*omega*omega
        return
        end

      subroutine shIijlay(lyr,i,j,CUP,CDN,Iij,ind)
      implicit none

      integer lyr,i,j,ind
      complex*16 CUP,CDN,Iij
c-----
c     lyr - layer to be evaluated
c     i j   To get the Iij integral
c     omega - angular frequency
c     wavo  - wavenumber
c     CUP, CDN  - coefficients for the layer
c     Iij   - integral
c     ind   - 1 evaluate for layer
c             2 evaluate for lower halfspace
c-----

      integer NL
      parameter (NL=200)

      common /elmat/el,elinv,es,nub
      complex*16 el(NL,2,2)
      complex*16 elinv(NL,2,2)
      complex*16 es(NL)
      complex*16 nub(NL)

        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep

c-----
c     internal values
c-----
c        (9.7.2)
c-----
      complex*16 gv, fv

      if(ind.eq.1)then
         call getfv(fv,nub(lyr),d(lyr))
         call getgv(gv,nub(lyr),d(lyr))
      Iij  = el(lyr,i,1)*el(lyr,j,1)*CUP*CUP*fv 
     1  + ( el(lyr,i,1)*el(lyr,j,2)+el(lyr,i,2)*el(lyr,j,1))*CUP*CDN*gv
     2  + el(lyr,i,2)*el(lyr,j,2)*CDN*CDN*fv
      else if (ind.eq.2)then
           Iij = el(lyr,i,2)*el(lyr,j,2)*CDN*CDN/(2.*nub(lyr))
      endif

      return
      end

        subroutine ddh(c,f,sumi1,UT,TT,dcdh)
        implicit none

      integer NL
      parameter (NL=200)
c-----
c       get the dcdh partials
c-----
        real*8 f,c
        real*8 sumi1
        real*8 UT(NL), TT(NL)
        real*8 dcdh(NL)
      
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep


      common/modlly/mmax
      integer mmax

c-----
c       internal variables
c-----
        real omega, omega2
        real wvno, wvno2
        real drho, dmu, dvdz, fac, ale
        real dfac
        integer k
        real mu1, mu2
        integer i, j, llflag
        real sumd

        omega = sngl(6.2831853*f)
        omega2 = omega*omega
        wvno = sngl(omega/c)
        wvno2 = wvno*wvno

c-----
c       define partial with respect to layer thickness
c-----
c       fac = 0.5d+00*c**3/(omega2*sumi1)
       
        ale=sngl(0.5d+00/sumi1)
 
        fac = sngl(ale*c/wvno2)
        llflag = 0
        do  k=1,mmax
             if(b(k).ne.0.0)then
                 if(llflag.eq.0)then
                     drho = rho(k)
                     dmu  = rho(k)*b(k)*b(k)
                     dvdz = 0.0
                 else 
                     mu2  = rho(k)*b(k)*b(k)
                     mu1  = rho(k-1)*b(k-1)*b(k-1)
                     drho = rho(k) - rho(k-1)
                     dmu  = mu2 - mu1
                     dvdz = sngl(TT(k)*TT(k)*(1.0/mu2 - 1.0/mu1))
                 endif
                 dfac = fac * sngl(( UT(k)*UT(k))*
     1               (omega2*drho - wvno2*dmu) + dvdz)
                 if(abs(dfac).lt.1.0e-38)then
                     dcdh(k) = 0.0
                 else
                     dcdh(k) = (dfac)
                 endif
                 llflag = llflag + 1
             else
                 dcdh(k) = 0.0
             endif
        enddo
c-----  
c           up to this point the dcdh are changes to phase velocity if
c           if the layer boundary changes. Here we change this to mean
c           the dc/dh for a change in layer thickness
c    
c           A layer becomes thicker if the base increases and the top
c           decreases its position. The dcdh to this point indicates
c           the effect of moving a boundary down. Now we convert to
c           the effect of changing a layer thickness.
c-----  
            do i=1,mmax-1
                sumd = 0.0
                do j=i+1,mmax
                    sumd = sumd + sngl(dcdh(j))
                enddo
                dcdh(i) = sumd
CDBH                WRITE(6,*)i,dcdh(i)
            enddo
            dcdh(mmax) = 0.0
        return
        end
