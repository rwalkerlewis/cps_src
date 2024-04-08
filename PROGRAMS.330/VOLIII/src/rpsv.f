      subroutine regn(f,c,UZ,UR,TZ,TR,dcda,dcdb,dcdr,dcdh,
     1    AR,gam,ccausal,u,
     1   flagr,sumi0,sumi1,sumi2,sumi3)
      implicit none

      integer NL
      parameter (NL=200)
      real*8 f,c
      real*8 UZ(NL), TZ(NL)
      real*8 UR(NL), TR(NL)
      real*8 dcda(NL),dcdb(NL), dcdr(NL), dcdh(NL)
      real*8 AR
      real*8  gam
      real*8  ccausal
      real*8 u
      real*8 flagr,sumi0, sumi1,sumi2,sumi3



      common/downpsv/RD, TD
      complex*16 RD(NL,2,2), TD(NL,2,2)

      common/uppsv/RU, TU
      complex*16 RU(NL,2,2), TU(NL,2,2)

      common/coeffinalpsv/CDN,CUP
      complex*16 CDN(NL,2), CUP(NL,2)

      common /ermat/er,erinv,ep,es,elog,nua,nub
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL),es(NL)
      complex*16 nua(NL),nub(NL)
      real elog(NL)


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

      common/bcpsv/alppsv,betpsv
      complex*16 alppsv(4,4), betpsv(4,4)

      integer jtop, jbot
      
      complex*16 om, wvno
      integer j, jlyr

      complex*16 TUR(NL),TUZ(NL),TTZ(NL),TTR(NL)
      complex*16 TUZ0
      complex*16 cdval(2)
      complex*16  cdcdb(NL), cdcda(NL), cdcdr(NL)
      complex*16 csumi0, csumi1, csumi2, csumi3

      integer i
      integer nbc
      complex*16 apsv(4,4)
      real mu, lam, l2mu
      real*8 x, dc, pi
      complex*16 IUZUZ, IURUR, IUZDUR, IURDUZ
      complex*16 IDURDUR, IDUZDUZ
      complex*16 I11,I22,I33,I44,I24,I14,I13
      real kr, or
      real z
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
          call doer(om,wvno)
c-----
c         define the boundary matrices
c         jbot = 3 halfspace
c         jtop = 1 free surface
c-----
          jbot = 3
          call dobotpsv(alppsv, jbot)
          jtop = 1
          call dotoppsv(betpsv, jtop)
c-----
c         evaluate the reflection coeffcients
c-----
          call botuppsv()
          call topdownpsv()
C         WRITE(6,*)'TU'
C         DO I=1,MMAX-1
C            WRITE(6,11)I,TU(i,1,1),TU(i,1,2),TU(i,2,1),TU(i,2,2)
C         ENDDO
C         WRITE(6,*)'RU'
C         DO I=1,MMAX-1
C            WRITE(6,11)I,RU(i,1,1),RU(i,1,2),RU(i,2,1),RU(i,2,2)
C         ENDDO
C         WRITE(6,*)'TD'
C         DO I=1,MMAX-1
C            WRITE(6,11)I,TD(i,1,1),TD(i,1,2),TD(i,2,1),TD(i,2,2)
C         ENDDO
C         WRITE(6,*)'RD'
C         DO I=1,MMAX-1
C            WRITE(6,11)I,RD(i,1,1),RD(i,1,2),RD(i,2,1),RD(i,2,2)
C         ENDDO

c-----
c         now do the recursion to get the Cu Cd and then the solution
c      jcond
c        1    top                    use jref = 1
c        2    within the layer stack use jref = 1,...,mmax-1
c        3    bottom                 use jref = mmax-1
c      jref   reference layer
c            
c-----
          if(jref.eq.0)then
               jlyr = 1
               call  cdinit(cdval,1,jlyr)
          else if(jref.ge.mmax)then
               jlyr = mmax-2 
               call cdinit(cdval,3,jlyr)
          else
               jlyr =jref
               call cdinit(cdval,2,jlyr)
          endif
            



C  10 format(a,2i5,4(1x,2e10.3,1x))
c-----
c-----
c          initialize and then apply recursion (9.8.6)(
c-----
          CDN(jlyr,1) = cdval(1)
          CDN(jlyr,2) = cdval(2)
          CUP(jlyr,1) = RD(jlyr,1,1)*CDN(jlyr,1) 
     1                + RD(jlyr,1,2)*CDN(jlyr,2)
          CUP(jlyr,2) = RD(jlyr,2,1)*CDN(jlyr,1) 
     1                + RD(jlyr,2,2)*CDN(jlyr,2)
          do j=jlyr-1,1,-1
             CUP(j,1) = TU(j+1,1,1)*CUP(j+1,1)
     1                + TU(j+1,1,2)*CUP(j+1,2)
             CUP(j,2) = TU(j+1,2,1)*CUP(j+1,1)
     1                + TU(j+1,2,2)*CUP(j+1,2)
             CDN(j,1) = RU(j  ,1,1)*CUP(j  ,1)
     1                + RU(j  ,1,2)*CUP(j  ,2)
             CDN(j,2) = RU(j  ,2,1)*CUP(j  ,1)
     1                + RU(j  ,2,2)*CUP(j  ,2)
          enddo
          do j=jlyr+1,mmax-1
             CDN(j,1) = TD(j-1,1,1)*CDN(j-1,1)
     1                + TD(j-1,1,2)*CDN(j-1,2)
             CDN(j,2) = TD(j-1,2,1)*CDN(j-1,1)
     1                + TD(j-1,2,2)*CDN(j-1,2)
             CUP(j,1) = RD(j  ,1,1)*CDN(j  ,1)
     1                + RD(j  ,1,2)*CDN(j  ,2)
             CUP(j,2) = RD(j  ,2,1)*CDN(j  ,1)
     1                + RD(j  ,2,2)*CDN(j  ,2)
          enddo
c-----
c         Get the coeffcients for the halfspace
c         CUP = 0 CDN = Td CDN
c----
          CDN(mmax,1) = TD(mmax-1,1,1)*CDN(mmax-1,1) 
     1                        + TD(mmax-1,1,2)*CDN(mmax-1,2)
          CDN(mmax,2) = TD(mmax-1,2,1)*CDN(mmax-1,1) 
     1                        + TD(mmax-1,2,2)*CDN(mmax-1,2)
          CUP(mmax,1) = dcmplx(0.0d+00,0.0d+00)
          CUP(mmax,2) = dcmplx(0.0d+00,0.0d+00)
c-----
c     safety
c-----   
          TUZ0=dcmplx(1.0d+00,0.0d+00)
c-----
c         now get the eigenfunctions (9.8.1)
c         B = E L C
c-----
C  12 format(a,i5,4(1x,2e10.3,1x))
C  13 format('*',a,i4,4(1x,2e10.3,1x))
C         WRITE(6,*)'CU1,Cu2,CD1,CD2'
          do j=1,mmax
C       WRITE(6,21)j,CUP(j,1),CUP(j,2),CDN(j,1),CDN(j,2)
C  21 format(i5,8e12.4)

               TUR(j) = er(j,1,1)*ep(j)*CUP(j,1)
     1                + er(j,1,2)*es(j)*CUP(j,2)
     2                + er(j,1,3)      *CDN(j,1)
     3                + er(j,1,4)      *CDN(j,2)
               TUZ(j) = er(j,2,1)*ep(j)*CUP(j,1)
     1                + er(j,2,2)*es(j)*CUP(j,2)
     2                + er(j,2,3)      *CDN(j,1)
     3                + er(j,2,4)      *CDN(j,2)
               TTZ(j) = er(j,3,1)*ep(j)*CUP(j,1)
     1                + er(j,3,2)*es(j)*CUP(j,2)
     2                + er(j,3,3)      *CDN(j,1)
     3                + er(j,3,4)      *CDN(j,2)
               TTR(j) = er(j,4,1)*ep(j)*CUP(j,1)
     1                + er(j,4,2)*es(j)*CUP(j,2)
     2                + er(j,4,3)      *CDN(j,1)
     3                + er(j,4,4)      *CDN(j,2)
               if(j.eq.1)then
                    TUZ0 = TUZ(j)
               endif
           enddo
           do j=1,mmax
               UR(j) = dreal(TUR(j)/TUZ0)
               UZ(j) = dreal(TUZ(j)/TUZ0)
               TZ(j) = dreal(TTZ(j)/TUZ0)
               TR(j) = dreal(TTR(j)/TUZ0)
           enddo
C  11  format(i5,4(1x,2e10.3,1x))
C   1     format(i5,1pd14.6,1pd14.6,1pd14.6,1pd14.6)
c-----
c         get the integrals
c-----
          csumi0 = 0.0d+00
          csumi1 = 0.0d+00
          csumi2 = 0.0d+00
          csumi3 = 0.0d+00
c-----
c     layers
c-----
       do i=1,mmax
           if(i.lt.mmax)then
              nbc = 1
           else
              nbc = 2
           endif
           call  amatpsv(om,wvno,apsv,i)
           mu = rho(i) *b(i) *b(i)
           l2mu  = rho(i)*a(i)*a(i) 
           lam  = l2mu - 2.*mu
           call psvIijlay(i,1,1,CUP,CDN,I11,nbc)
           call psvIijlay(i,2,2,CUP,CDN,I22,nbc)
           call psvIijlay(i,3,3,CUP,CDN,I33,nbc)
           call psvIijlay(i,4,4,CUP,CDN,I44,nbc)
           call psvIijlay(i,1,4,CUP,CDN,I14,nbc)
           call psvIijlay(i,2,4,CUP,CDN,I24,nbc)
           call psvIijlay(i,1,3,CUP,CDN,I13,nbc)
c   INT U_z^2
           IUZUZ = I22
c   INT U_r^2
           IURUR = I11
c   INT [ dU_z/dz]^2
           IDUZDUZ = apsv(2,1)*apsv(2,1)*I11
     1             + 2.*apsv(2,1)*apsv(2,3)*I13
     2             + apsv(2,3)*apsv(2,3)*I33
c   INT [ dU_r/dz]^2
           IDURDUR = apsv(1,2)*apsv(1,2)*I22
     1             + 2.*apsv(1,2)*apsv(1,4)*I24
     2             + apsv(1,4)*apsv(1,4)*I44
c   INT U_r dU_z/dz
           IURDUZ  = apsv(2,1)*I11 + apsv(2,3)*I13
c   INT U_z dU_r/dz
           IUZDUR  = apsv(1,2)*I22 + apsv(1,4)*I24
           csumi0 = csumi0 + rho(i)*(IURUR + IUZUZ)
           csumi1 = csumi1 + mu*IUZUZ + l2mu*IURUR
           csumi2 = csumi2 + mu*IUZDUR - LAM*IURDUZ
           csumi3 = csumi3 + mu*IDURDUR + l2mu*IDUZDUZ
           cdcdb(i) = rho(i)*b(i)*(4*IURDUZ/wvno +
     1         IUZUZ + 2*IUZDUR/wvno + IDURDUR/wvno**2)
           cdcda(i) = rho(i)*a(i)*(IURUR - 2.*IURDUZ/wvno +
     1                IDUZDUZ/wvno**2)
           cdcdr(i) = 0.5*a(i)*cdcda(i)/rho(i)
     1              + 0.5*b(i)*cdcdb(i)/rho(i)
     2              - 0.5*c*c*(IUZUZ + IURUR)
        enddo
c-----
c         get group velocity
c-----
        sumi0 = dreal(csumi0/TUZ0**2)
        sumi1 = dreal(csumi1/TUZ0**2)
        sumi2 = dreal(csumi2/TUZ0**2)
        sumi3 = dreal(csumi3/TUZ0**2)
        kr = sngl(dreal(wvno))
        or = sngl(dreal(om))
        u =  (kr*sumi1+sumi2)/(or*sumi0)
        AR=1.0/(2*c*u*dreal(csumi0/TUZ0**2))
c-----
c       get the dcdh
c-----
        call ddh(c,f,AR,UR,UZ,TZ,TR,dcdh)
c-----
c       get final partials and gamma
c-----
        gam = 0.0
        dc = 0.0 
        pi = 3.141592653589493d+00
        do i=1,mmax
           dcda(i) = dreal(cdcda(i)/csumi0)/u
           dcdb(i) = dreal(cdcdb(i)/csumi0)/u
           dcdr(i) = dreal(cdcdr(i)/csumi0)/u
      
           x = 0.0
           if(qb(i).gt.1.0)then
              x = x + dcdb(i) * b(i)/qb(i) 
           else
              x = x + dcdb(i) * b(i)*qb(i)
           endif
           if(qa(i).gt.1.0)then
              x = x + dcda(i) * a(i)/qa(i) 
           else
              x = x + dcda(i) * a(i)*qa(i)
           endif
           gam = gam + x
           dc = dc + dlog(f/frefs(i))*x/pi 
        enddo
        gam=0.5*6.2831853*f*gam/(c*c)
        ccausal = c + dc
C  31    format(//22x,'RAYLEIGH WAVE      MODE #',i3/
C    1  '        T =',e11.4,' C =   ',e11.4,' U   =',e11.4/
C    2  '        AR=',e11.4,' GAMMA=',e11.4,' ZREF=',e11.4/
C    3  '    M       UR         TR        UZ'
C    4  ,'         TZ        DC/DH'
C    5  ,'      DC/DA      DC/DB      DC/DR')
C       write(6,331)mode,1./f,c,u,AR,gam,refdep
C 331    format(22x,'RAYLEIGH WAVE      MODE #',i3/
C    1  '        T =',e11.4,' C =   ',e11.4,' U   =',e11.4/
C    2  '        AR=',e11.4,' GAMMA=',e11.4,' ZREF=',e11.4/
C    3  '    M          Z      UR         TR        UZ'
C    4  ,'         TZ        DC/DH'
C    5  ,'      DC/DA      DC/DB      DC/DR')
        flagr = dreal(om)*dreal(om)*sumi0 -
     1     dreal(wvno)*dreal(wvno)*sumi1 -
     2     2*dreal(wvno)*sumi2           -
     3     sumi3
C       write(6,*)'L=',flagr, 'L/om2=',flagr/dreal(om)**2
C                WRITE(6,212)1/f,c,u,gam,
C    1                sumi0,sumi1,sumi2,sumi3,
C    2                flagr,ar,ur(1)/uz(1),dreal(flagr/(sumi0*om**2))
C       write(6,331)mode,1./f,c,u,AR,gam,refdep
C   2   format(i5,1x,8e11.3)
C  22   format(i5,1x,f10.3,8e11.3)

        z = 0.0
        do i=1,mmax
C       write(6,2)i,ur(i),tr(i),uz(i),tz(i),dcdh(i),
C    1      dcda(i),dcdb(i),dcdr(i)
C       write(6,22)i,z,ur(i),tr(i),uz(i),tz(i),dcdh(i),
C    1      dcda(i),dcdb(i),dcdr(i)
        z = z + d(i)
        enddo
C 212    format(' T=',e15.7,'  C=',e15.7,'  U=',e15.7,'  G=',e15.7/    
C    1          'I0=',e15.7,' I1=',e15.7,' I2=',e15.7,' I3=',e15.7/
C    2          ' L=',e15.7,' AR=',e15.7,'  E=',f15.7,' L/om2 I0',e15.7)
      
          return
          end


       subroutine cdinit(cdval,jcond,jref)
c-----
c      initialize the P-SV CD value depending on where the
c      initialization is done
c
c      cdval   CD
c      jcond
c        1    top
c        2    within the layer stack
c        3    bottom
c      jref   reference layer
c-----
c      To obtain a starting value for CD it is necessary to call upon 
c      the fact that the secular funciton used is zero.
c      The secular function can be evaluated at the top, bottom or within the layer stack.
c
c      1 Top - this is written n the expanded form not as matrix multiplication Equation 9.8.2
c            for a free surface condition, e.g., H = I
c        |  z11 z12 | | cd1 |   | 0 |
c        |          | |     | = |   |
c        |  z21 z22 | | cd2 |   | 0 |
c        where
c           z11 = e31 ep RD11 + e32 es RD21 + e33
c           z21 = e41 ep RD11 + e42 es RD21 + e43
c           z12 = e31 ep RD12 + e32 es RD22 + e34
c           z22 = e41 ep RD12 + e42 es RD22 + e44
c        ep = e^{- \nu_alpha_1 d_1 }  es = e^{- \nu \beta_1 d_1}
c
c      2 Within
c         Equation 9.8.5
c           
C         [ I - R_U^j R_D^j)]C_D^j = 0
c         Note this will not work if the complete model is a uniform halfspace
c         since the R_D == 0
c      3 bottom written an expanded form (9.8.3) for an elastic halfspace,. e.g., G = E_mmax^-1
c        and here we let g = E_{mmax}^-1 E_{mmax-1}
c
c
c        or
c           z11 = g11 + g13 ep RU11  + g14 es Ru21
c           z12 = g12 + g13 ep RU12  + g14 es Ru22
c           z21 = g21 + g23 ep RU11  + g24 es Ru21
c           z22 = g22 + g23 ep RU12  + g24 es Ru22
c        ep = e^{- \nu_alpha_{N-1} d_{N-1} }  es = e^{- \nu \beta_{N-1} d_{N-1}}
c
c      In these secualr equations,
c      the  C_D is a 2x1 and the matrix within the square brackets is 2x2. From Exercise 9.4, if the system is
c
c      | a   b | | x |    | 0 |
c      |       | |   | =  |   |
c      | c   d | | y |    | 0 |
c     
c      this has meaning if the determinaant = 0 AND 
c
c      | x |             1          |  b |
c      |   |  =  _____________      |    |
c      | y |     [a^2 + b^2]^{1/2}  | -a |
c
       implicit none
c
       complex*16 cdval(2)
       integer jcond, jref

        integer NL
        parameter (NL=200)
      common/downpsv/RD, TD
      complex*16 RD(NL,2,2), TD(NL,2,2)
        
      common/uppsv/RU, TU
      complex*16 RU(NL,2,2), TU(NL,2,2)

      common /ermat/er,erinv,ep,es,elog,nua,nub
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL),es(NL)
      complex*16 nua(NL),nub(NL)
      real elog(NL)

      common/modlly/mmax
      integer mmax
c-----
c     internal variables
c-----
       complex*16 c3(3,3)
       real*8 cnorm
       complex*16 g(2,4)
       integer i,j,k
       complex*16 z1,z2

     
c----
c      initialize - note the use of the matrix algebra routines would be interbnally of
c           this were programmed in FORTRAN90 and later
c-----


       if(jcond.eq.1)then
            c3(1,1) = er(1,3,1)*ep(1)*RD(1,1,1) 
     1              + er(1,3,2)*es(1)*RD(1,2,1) + er(1,3,3)
            c3(2,1) = er(1,4,1)*ep(1)*RD(1,1,1)
     1              + er(1,4,2)*es(1)*RD(1,2,1) + er(1,4,3)
            c3(1,2) = er(1,3,1)*ep(1)*RD(1,1,2) 
     1              + er(1,3,2)*es(1)*RD(1,2,2) + er(1,3,4)
            c3(2,2) = er(1,4,1)*ep(1)*RD(1,1,2) 
     1              + er(1,4,2)*es(1)*RD(1,2,2) + er(1,4,4)
       else if(jcond.eq.2)then
c-----
c           Form I = RU RD
c-----
            do i=1,2
                 do j=1,2
                     if(i.eq.j)then
                        c3(i,j) = dcmplx(1.0d+00,0.0d+00)
                     else
                        c3(i,j) = dcmplx(0.0d+00,0.0d+00)
                     endif
                     do k=1,2
                       c3(i,j) = c3(i,j) -RU(jref,i,k)*RD(jref,k,j)
                     enddo
                 enddo
             enddo
       else if(jcond.eq.3)then
c            GE
             do i=1,2
                do j=1,4
                   g(i,j) = dcmplx(0.0d+00,0.0d+00)
                   do k=1,4
                       g(i,j) = g(i,j) + erinv(mmax,i,k)*er(mmax-1,k,j)
                   enddo
                enddo
             enddo
c-----
c           This will give Cu(N-1)
            z1 = g(1,1) + g(1,3)*ep(mmax-1)*RU(mmax-1,1,1)  
     1          + g(1,4)*es(mmax-1)*Ru(mmax-1,2,1)
            z2 = g(1,2) + g(1,3)*ep(mmax-1)*RU(mmax-1,1,2)  
     1          + g(1,4)*es(mmax-1)*Ru(mmax-1,2,2)
c           To get CD use CD = Ru CU
            c3(1,1) = RU(mmax-1,1,1)*z1 + RU(mmax-1,1,2)*z2
            c3(1,2) = RU(mmax-1,2,1)*z1 + RU(mmax-1,2,2)*z2
       endif
c-----
c           This part is common to all cases. 
c           now define the normalized vector - trhis is a bt odd since
c           The CD vectory is
c             (x1 + i y1) e_1 + ( x2 + i y2) e_2
c           Do extend the idea of a dot product of z1 e_1 + z2 e_2 to be
c                 z1 z1* + z2 z2*
c----
                z1 = c3(1,1)
                z2 = c3(1,2)
                cnorm = dsqrt(dreal(z1*dconjg(z1) + z2*dconjg(z2)))
                CDVAL(1) =   z2/cnorm
                CDVAL(2) = - z1/cnorm
       return
       end

       subroutine cjcopy(R,jref,C)
       implicit none
       integer NL
       parameter (NL=200)
       complex*16 R(NL,2,2)
       integer jref
       complex*16 C(2,2)
c-----
c      internal variables
c-----
       integer i,j
       do i=1,2  
          do j=1,2
              c(i,j) = R(jref,i,j)
          enddo
       enddo
       return
       end

       subroutine csub(a,b,c)
       implicit none
c-----
c      c = a - b
c-----
       complex*16 a(2,2), b(2,2), c(2,2)
c-----
c      internal variables
c-----
       integer i,j
       do i=1,2
          do j=1,2
              c(i,j) = a(i,j) - b(i,j)
          enddo
       enddo
       return
       end

       subroutine amatpsv(omega,wvno,apsv,lyr)
       implicit none
c-----
c         get the SH A matrix of
C         dB/dz = A B
c-----

       complex*16 omega, wvno
       complex*16 apsv(4,4)
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
        complex*16 cdzero, cdone
        real mu,lam,l2mu
        integer i,j
c-----
c       initialize
c-----
        cdzero = dcmplx(0.0d+00,0.0d+00)
        cdone  = dcmplx(1.0d+00,0.0d+00)
        do i=1,4
           do j=1,4
              apsv(i,j) = cdzero
           enddo
        enddo
          

        mu = rho(lyr)*b(lyr)*b(lyr)
        l2mu = rho(lyr)*a(lyr)*a(lyr)
        lam = l2mu - 2*mu

        apsv(1,2) = -wvno
        apsv(1,4) = cdone/mu
        apsv(2,1) = wvno*lam/l2mu
        apsv(2,3) = cdone/l2mu
        apsv(3,2) = - rho(lyr)*omega*omega
        apsv(3,4) = - apsv(1,2)
        apsv(4,1) = - rho(lyr)*omega*omega 
     1              + 4*wvno*wvno*mu*(lam+mu)/l2mu
        apsv(4,3) = - apsv(2,1)
        return
        end

      subroutine psvIijlay(lyr,i,j,CUP,CDN,Iij,ind)
      implicit none

      integer NL
      parameter (NL=200)

      integer lyr,i,j,ind
      complex*16 CUP(NL,2),CDN(NL,2),Iij
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
      common /ermat/er,erinv,ep,es,elog,nua,nub
      complex*16 er(NL,4,4)
      complex*16 erinv(NL,4,4)
      complex*16 ep(NL),es(NL)
      real elog(NL)
      complex*16 nua(NL),nub(NL)

        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep

c-----
c     internal values
c-----
      complex*16 C(4)
      complex*16 EE(4,4)
      complex*16 II(4,4)

      complex*16 gva, fva
      complex*16 gvb, fvb
      complex*16 fvapb, fvamb
      complex*16 cdtmp
 
c-----
c   (9.7.1)
c   Note K_{m-1}^{D,P}  = CDN(1)
c   Note K_{m-1}^{D,SV} = CDN(2)
c   Note K_{m  }^{U,P}  = CUP(1)
c   Note K_{m  }^{U,SV} = CUP(2)
c   this is message since (9.7.1) has  8 terms
c
c   We want to evaluate the layer integrals
c   Iij = bi bj where the  bi and bj are the eigenfunctions
c   Recall that b_i = E_im L_mm Cm
c   where L is a diagonal matrix depending on z and C = [Cup, Cdn]
c   So  INT b_i b_j = SUM SUM E_im E_jn INT (L_mm L_nn) C_mm C_nn
c   Rather than writing out the expressions in long hand, use some
c   simpler terms to make it easier to write the final expression
c-----

      C(1) = CUP(lyr,1)
      C(2) = CUP(lyr,2)
      C(3) = CDN(lyr,1)
      C(4) = CDN(lyr,2)
      EE(1,1) = er(lyr,1,1)
      EE(1,2) = er(lyr,1,2)
      EE(1,3) = er(lyr,1,3)
      EE(1,4) = er(lyr,1,4)
      EE(2,1) = er(lyr,2,1)
      EE(2,2) = er(lyr,2,2)
      EE(2,3) = er(lyr,2,3)
      EE(2,4) = er(lyr,2,4)
      EE(3,1) = er(lyr,3,1)
      EE(3,2) = er(lyr,3,2)
      EE(3,3) = er(lyr,3,3)
      EE(3,4) = er(lyr,3,4)
      EE(4,1) = er(lyr,4,1)
      EE(4,2) = er(lyr,4,2)
      EE(4,3) = er(lyr,4,3)
      EE(4,4) = er(lyr,4,4)

c-----
c     The diagonal elements of L are
c    [e^{- \nu_\alpha (z_j -z ) ,
c    [e^{- \nu_\beta  (z_j -z ) ,
c    [e^{- \nu_\alpha (z -z_{j-1} ) ,
c    [e^{- \nu_\beta  (z -z_{j-1} ) ]
c
c    where z_j - z_{j-1} = d_m, the thikness of layer m. The integrals requirted will be
c    Imn = INT L_m L_n We will need 10 such integrals
c-----

      if(ind.eq.1)then
         call getfv(fva,nua(lyr),d(lyr))
         call getgv(gva,nua(lyr),d(lyr))
         call getfv(fvb,nub(lyr),d(lyr))
         call getgv(gvb,nub(lyr),d(lyr))
         cdtmp = nua(lyr) + nub(lyr)
         call getfvab(fvapb,cdtmp,d(lyr))
         fvamb = (es(lyr) -ep(lyr))/(nua(lyr) - nub(lyr))

         II(1,1) = fva
         II(1,2) = fvapb
         II(1,3) = gva
         II(1,4) = fvamb

         II(2,1) = fvapb
         II(2,2) = fvb
         II(2,3) = fvamb
         II(2,4) = gvb

         II(3,1) = gva
         II(3,2) = fvamb
         II(3,3) = fva
         II(3,4) = fvapb

         II(4,1) =  fvamb
         II(4,2) = gvb
         II(4,3) =  fvapb
         II(4,4) = fvb

         Iij = 
     1                           EE(i,1)*EE(j,1)*II(1,1)*C(1)*C(1)
     1       + (EE(i,1)*EE(j,2) + EE(i,2)*EE(j,1)) *II(1,2)*C(1)*C(2)
     1       + (EE(i,1)*EE(j,3) + EE(i,3)*EE(j,1)) *II(1,3)*C(1)*C(3)
     1       + (EE(i,1)*EE(j,4) + EE(i,4)*EE(j,1)) *II(1,4)*C(1)*C(4)
     1       +                      EE(i,2)*EE(j,2)*II(2,2)*C(2)*C(2)
     1       +                      EE(i,3)*EE(j,3)*II(3,3)*C(3)*C(3)
     1       +                      EE(i,4)*EE(j,4)*II(4,4)*C(4)*C(4)
     1       + (EE(i,2)*EE(j,3) + EE(i,3)*EE(j,2)) *II(2,3)*C(2)*C(3)
     1       + (EE(i,2)*EE(j,4) + EE(i,4)*EE(j,2)) *II(2,4)*C(2)*C(4)
     1       + (EE(i,3)*EE(j,4) + EE(i,4)*EE(j,3)) *II(3,4)*C(3)*C(4)
      else if (ind.eq.2)then
c-----
c     only include the CD terms or C(3), C(4) and assume that d_m -> \infty
c-----
           II(3,3) = dcmplx(1.0d+00,0.0d+00)/(2.*nua(lyr) )
           II(4,4) = dcmplx(1.0d+00,0.0d+00)/(2.*nub(lyr) )
           II(3,4) = dcmplx(1.0d+00,0.0d+00)/(nua(lyr)+nub(lyr))
         Iij = 
     1       +  EE(i,3)*EE(j,3)*II(3,3)*C(3)*C(3)
     1       +  EE(i,4)*EE(j,4)*II(4,4)*C(4)*C(4)
     1       + (EE(i,3)*EE(j,4) + EE(i,4)*EE(j,3)) *II(3,4)*C(3)*C(4)
      endif

      return
      end

        subroutine ddh(c,f,AR,UR,UZ,TZ,TR,dcdh)
        implicit none

      integer NL
      parameter (NL=200)
c-----
c       get the dcdh partials
c-----
        real*8 f,c
        real*8 AR
        real*8 UR(NL), UZ(NL), TZ(NL), TR(NL), dcdh(NL)
      
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
        real omega, om2
        real wvno, wvno2
        real fac
        integer i, j
        real*8  sumd
        real*8 tur, tuz, ttz, ttr
        real*8 dfac, gfac1, gfac2, gfac3, gfac4, gfac5, gfac6
        
        real*8 drho, dmu, dlm
        real*8 duzdzp, dlur2
        real*8 dl2mu, duzdzm, drur2
        real*8 URB, DURDZM, DURDZP
        real xmu, xmum1
        real xlam, xlamm1
        real xl2mm, xl2mp
        integer m

        omega = sngl(6.2831853*f)
        om2 = omega*omega
        wvno = sngl(omega/c)
        wvno2 = wvno*wvno

        fac = sngl(AR*c/wvno2 )

c-----

        do  m=1,mmax
            
            tuz = uz(m)
            ttz = tz(m)
            ttr = tr(m)
            if(b(m).eq.0.0)then
                tur = -wvno*ttz/(rho(m)*om2)
            else
                tur = ur(m)
            endif
c-----
c       this assumes that the top is a halfspace
c-----
                xmu = rho(m)*b(m)*b(m)
                xlam = rho(m)*a(m)*A(m) - 2.* xmu
            if(m.eq.1)then
                xmu = rho(m)*b(m)*b(m)
                xlam = rho(m)*a(m)*A(m) - 2.* xmu
                drho = rho(1) - 0.0
                dmu  = xmu - 0.0
                dlm = xlam - 0.0
                dl2mu = dlm + dmu + dmu
                xl2mp = xlam + xmu + xmu
                duzdzp = (ttz + wvno*xlam*tur)/xl2mp
                if(b(m) .eq.0.0)then
                    durdzp = wvno*tuz
                else
                    durdzp = (ttr/xmu) - wvno*tuz
                endif
                drur2 = tur*tur*drho
                dlur2 = tur*tur*dl2mu

                gfac1 =  om2*drho*tuz**2
                gfac2 =  om2*drur2
                gfac3 = -wvno2*dmu*tuz**2
                gfac4 = -wvno2*dlur2
                gfac5 =  (xl2mp*duzdzp**2)
                gfac6 =  (xmu*durdzp**2 )
            else
                xmum1 = rho(m-1)*b(m-1)*b(m-1)
                xlamm1 = rho(m-1)*a(m-1)*a(m-1) -2.*xmum1
                drho = rho(m) - rho(m-1)
                dmu = xmu - xmum1
                dlm = xlam - xlamm1
                dl2mu = dlm + dmu + dmu
                xl2mp = xlam   + xmu   + xmu
                xl2mm = xlamm1 + xmum1 + xmum1
                duzdzp = (ttz + wvno*xlam*tur)/xl2mp
                if(xmu.eq.0.0)then
                    durdzp = wvno*tuz
                else
                    durdzp = (ttr/xmu) - wvno*tuz
                endif
                if(xmum1.eq.0.0 )then
                    durdzm = wvno*tuz
                else
                    durdzm = (ttr/xmum1) - wvno*tuz
                endif
c-----
c       attempt to fix for water layer, since Ur is not continuous
c       across fluid - solid interface or fluid-fluid interface
c-----
                if(b(m-1).eq.0.0 .and. b(m).ne.0 )then
                    URB = -wvno*tz(m)/(rho(m-1)*om2)
                    drur2 = tur*tur*rho(m)-URB*URB*rho(m-1)
                    dlur2 = tur*tur*xl2mp -URB*URB*xl2mm
                    duzdzm = (ttz + wvno*xlamm1*URB)/xlamm1
                else if(b(m-1).eq.0.0 .and. b(m).eq.0.0 )then
                    URB = -wvno*tz(m)/(rho(m-1)*om2)
                    drur2 = tur*tur*rho(m)- URB*URB*rho(m-1)
                    dlur2 = tur*tur*xl2mp  - URB*URB*xl2mm
                    duzdzm = (ttz + wvno*xlamm1*URB)/xl2mm
                else
                    drur2 = tur*tur*drho
                    dlur2 = tur*tur*dl2mu
                    duzdzm = (ttz + wvno*xlamm1*tur)/xl2mm
                endif


                gfac1 =  om2*drho*tuz**2
                gfac2 =  om2*drur2
                gfac3 = -wvno2*dmu*tuz**2
                gfac4 = -wvno2*dlur2
                gfac5 =  (xl2mp*duzdzp**2-xl2mm*duzdzm**2)
                gfac6 =  (xmu*durdzp**2 - xmum1*durdzm**2)
            endif
            dfac = fac * (
     1          gfac1 + gfac2 + gfac3 + gfac4
     2          + gfac5 + gfac6 )
            if(dabs(dfac).lt.1.0d-80)then
                dfac = 0.0d+00
            endif
            dcdh(m) = dfac
        enddo

 
C       fac = ale*c/wvno2
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
                    sumd = sumd + dcdh(j)
                enddo
                dcdh(i) = sumd
            enddo
            dcdh(mmax) = 0.0
        return
        end
