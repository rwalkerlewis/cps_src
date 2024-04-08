        program hudson96
c---------------------------------------------------------------------c
c                                                                     c
c      COMPUTER PROGRAMS IN SEISMOLOGY                                c
c      VOLUME VI                                                      c
c                                                                     c
c      PROGRAM: HUDSON96                                              c
c                                                                     c
c      COPYRIGHT 2008                                                 c
c      R. B. Herrmann                                                 c
c      Department of Earth and Atmospheric Sciences                   c
c      Saint Louis University                                         c
c      221 North Grand Boulevard                                      c
c      St. Louis, Missouri 63013                                      c
c      U. S. A.                                                       c
c                                                                     c
c---------------------------------------------------------------------c
c       CHANGES
c       17 MAR 2008 - created this program for the purpose of 
c             synthesizing a teleseismic P-wave signal
c             for a high sample rate. To save time, the
c             approach of Hudson (1969) is used that
c             separates the problem into a source crust,
c             receiver crust and teleseism problems, with the source
c             and receiver crust response fcomputed using propagator
c             matrices.
c
c             This program may seem  complicated, primarily because
c             large sections of the programs hspec96p, hrftn96 and
c             time96 are used as the building blocks of this program.
c             However it does make synthetics fast
c          I recommend using the P wave synthetics from first arrival to
c          the fime of PP or PcP, which ever is earliest.
c
c          I do not recommend the use of Z and R (e.g., SV) S synthetics
c          since too much is missing.
c
c          The SH synthetics can be used up to the arrival of SS, or ScS,
c          whichever is earliest. However do not use part 80 degrees
c          this the timing and geometrical spreading are that of SKS and
c          not S
c
c       11 APR 2008 - corrected an error that led to incorrect synthetics
c          if the source was at the default 60 km of the source model
c       25 MAY 2008 - SH times are written into the header so that they 
c          appear in the SAC file
c       28 MAY 2008 - to properly set the arrival time, the start of the
c          teleseismic attenuation pulse is determined empirically using
c          two FFT's. This must be done using the same DT as for the 
c          synthetic because of the frequency dispersion 
c          inherent in the causal T* 
c          operator
c       29 MAY 2008 - corrected indexing when merging source and 
c          receiver structure with teleseism structure. The halfspace 
c          was not inserted correctly
c       18 FEB 2009 -  took care to add common/depref/refdep in frstar
c       17 MAR 2009 -  the default offset of P is 10 sec, and
c                      is NOW 20 sec not 10 sec for S
c                      corrected numerical underflow for S-wave incident
c                      at the receiver by using a compound matrix
c       19 MAR 2009 -  added -TSTAR tstar to override the computed T*
c                      to permit the use of the same code, if T* = 0
c                      we use 0.0001 internally
c       07 MAY 2009 -  Add -TTONLY flag to get travel times only
c       23 JUL 2009 - ensured consistent call to FFT routine with
c               complex arguments. We now call zfour
c       11 OCT 2010  - forced extension of source/receiver models to
c                    a predefined maximum depth
c       01 JUN 2013 -  Modified subroutine to prevent NaN for direct 
c                      ray by replacing
c                      pnew = 0.999d+00 * pupper to
c                      pnew = 0.99d+00 * pupper
c                      also modified subroutine frstar to have the 
c                      dogeom argument if we do not want to compute 
c                      teleseismic geometrical spreading
c       13 MAY 2015 - fixed error in subroutine modmerge in which 
c                     the number of layers for the source and receiver 
c                     models were not specified correctly changed
c                        call getindex(layerdepth(i),layerdepth(i+1),Sd,Tmmax,k)
c                     to
c                        call getindex(layerdepth(i),layerdepth(i+1),Sd,Smmax,k)
c                     and changed
c                        call getindex(layerdepth(i),layerdepth(i+1),Rd,Tmmax,k)
c                     to
c                        call getindex(layerdepth(i),layerdepth(i+1),Rd,Rmmax,k)
c                     Also reworded logic on modmerge for k < 0 return
c       22 APR 2018 - correct defaults from DT=1024 to NPTS=1 to DT=1, 
c                     NPTS=1024
c       02 MAY 2018 - set the tstaradj from telop to zero since this 
c                     will make
c                     everything agree with wavenumber integration 
c                     synthetics.
c       28 SEP 2019 - This now also creates the Green's functions for a 
c                     point force
c        
c          models stored in the common blocks srcmod recmod
c       04 JUL 2023 - complete rewrite of code to follow the
c                     derivation in the book. 
c  TODO
c       
c  2. For the S first arrival do not permit SKS
c  6. redo frstar so that the model read is in the fstarr which
c      then returns the medium properties of the original model
c      This will make geometrical simpler
c  7. See effect of Earth flattening on p-tau.  Currently for
c     short distance and deep event, the ray parameters of P pP and sP
c     are different and hence the whole validity of things is questioned
c     the flattening is a kludge, but actually works for hspec96?
c     THIS IS REQUIRED
c  8. 
c----- 
        implicit none

c-----
c       input/output default LUN
c-----
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
c-----
c       velocity model file
c-----
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep

        integer mmax, iunit, iiso, iflsph, idimen, icnvel
        integer ierr
        character title*80
c-----
c       command line information
c-----
        character modeltel*100
        character modelsrc*100
        character modelrec*100

        real dt, gcarc, hs, offset, utstar
        integer npts

        logical dosrc, dorec, dotel, dop, dokjar, dottonly,verbose
        logical zdown

        common/control/verbose
c-----
c       teleseismic parameters
c-----
        real hstel, hrtel
        real zsrc, zrec
        real rayp,  geom,  tstar
        real raypp
        real rayps
        real vpsrc, vssrc, densrc
        real vprec, vsrec, denrec

        real tstaradj

        real TP, TSV
        real TPmantle, TSVmantle

        real srcdelay, recdelay
        real srcdist, recdist
        real srctime, rectime
c-----
c       Earth parameters
c-----
        common/earth/radius
        real radius
        real kmdeg
c-----
c       receiver region arrays
c-----
        integer MAXPTS
        parameter (MAXPTS=16384)
        complex recz(MAXPTS), recr(MAXPTS), rect(MAXPTS)
        complex srcp(21,MAXPTS)
c-----
c       teleseismic tstar operator
c-----
        complex zts(MAXPTS)
c-----
c       time series parameters
c-----
        real alpha
c-----
c       internal parameters
c-----
        integer lnblnk
        real ttime

c-----
c       ksrc = temporary array for jsrc for output
c       here jsrc != 0 reflects source information
c       if receiver is not in a fluid DO NOT output pressure field
c-----
        integer ksrc(21)
c-----
c       lsrc maps jsrc to output Green s functions. e.g., if
c       jsrc(8) = radial explosion term, but in final output it
c       occupies position 10, or jsrc(lsrc(10)) = computed
c-----
        integer*4 lsrc(21)
        data lsrc/1,2,3,4,13,5,6,14,7,8,9,10,11,12,15,16,17,18,19,20,21/

        

c-----
c       initialize
c-----
        radius = 6371.
        kmdeg=radius*3.1415927/180.0
c-----
c       parse command line arguments
c-----
        call gcmdln(dt, gcarc, hs, npts,
     1    modeltel,modelsrc,modelrec,offset,
     2    dosrc,dorec,dotel,dop,dokjar,utstar,dottonly,
     3    zsrc,zrec,verbose,zdown)
c-----
c       define alpha for time domain damping
c-----
        alpha = 2.5/(npts*dt)
C       ALPHA = 0.0
c-----
c       make npts a power of 2, and check the size
c-----
        call npow2(npts)
        if(npts.gt.MAXPTS)npts = MAXPTS
c-----
c       define moment/force system for wave type
c-----
        call dosetupsrc(dop,ksrc)
c-----
c       output command line 
c-----
        if(dop)then
             WRITE(LOT,*)'Synthetic seismogram parameters for P'
        else
             WRITE(LOT,*)'Synthetic seismogram parameters for S'
        endif
        WRITE(LOT,*)'Offset               : ',offset,' s'
        WRITE(LOT,*)'Source depth (km)    : ',hs    ,' km'
        WRITE(LOT,*)'Arc distance (deg)   : ',gcarc ,' deg'
        WRITE(LOT,*)'npts                 : ',npts
        WRITE(LOT,*)'dt                   : ',dt    ,' s'
        WRITE(LOT,*)'Teleseism model      : ',
     1                                 modeltel(1:lnblnk(modeltel))
        WRITE(LOT,*)'Source region model  : ',
     1                                 modelsrc(1:lnblnk(modelsrc))
        WRITE(LOT,*)'Receiver region model: ',
     1                                 modelrec(1:lnblnk(modelrec))
c-----
c       The order here is important
c       1. From the mixed model get the ray parameter and total
c          travel time from the source at depth hs to the receiver
c          at the surface
c       2. Using this ray parameter get the source and receiver
c          region response.  There are several important parameters
c          srcdelay - recdelay : time delays of the synthetics
c          srcdist  - recdist  : horizontal distance in flattened
c                        model from source/receiver to where ray 
c                        intersects the base of layered structure
c          srctime  - rectime  : ray time delay along this path
c       3. Get geometrical spreading and T* for the mantle
c          path for propagation of dist - srcdist - recdist
c          
c-----
        call defineprop(modeltel,modelsrc,modelrec,hs,
     1      hstel, hrtel, zsrc, zrec,
     2      vpsrc,vssrc,densrc, vprec,vsrec,denrec)
        WRITE(6,*)vpsrc,vssrc,densrc
        WRITE(6,*)vprec,vsrec,denrec

        WRITE(LOT,*)'Teleseism model used from for depths',
     1      ' greater than',hstel,' (km) at source'
        WRITE(LOT,*)'Teleseism model used to   for depths',
     1      ' greater than',hrtel,' (km) at receiver'
c-----
c     Get true travel time and rayp between source and surface receiver
c     The travel time is used to check the computations of
c     source + mantle + crust done in outputgreen
c-----
        WRITE(LOT,*)'Total travel time and ray parameter ',
     1   'from the source to reciever at surface:'

        if(dop)then
        call htrav(gcarc*111.195,hs,0.0,.true. ,TP ,raypp,tstar)
           rayp = raypp
           WRITE(LOT,*)'   TP        :',TP, ' s'
           WRITE(LOT,*)'   RAYPp     :',raypp,' s/km',
     1       111.195*raypp, 's/deg (at surface)'
        else 
        call htrav(gcarc*111.195,hs,0.0,.false.,TSV,rayps,tstar)
           rayp = rayps
           WRITE(LOT,*)'   TS        :',TSV, ' s'
           WRITE(LOT,*)'   RAYPs     :',rayps,' s/km',
     1       111.195*rayps, 's/deg (at surface)'
        endif
        WRITE(6,*)'Total T*',tstar
        WRITE(6,*)'tstar:',tstar
c-----
c       Get source region  input to teleseismic propagation
c-----
        WRITE(LOT,*)'-------------------------------------------'
        WRITE(LOT,*)'Computing source region input (negative layer',
     1    ' index indicates above the source) :'
        call getsrc(dt,alpha,rayp,npts,srcp,'hudsonsrc.mod',hs,hstel,
     1      srcdelay,srcdist,srctime,dosrc,vpsrc,vssrc,densrc,dop)
        WRITE(LOT,*)' P-velocity src base  : ',vpsrc,' (km/s)'
        WRITE(LOT,*)' S-velocity src base  : ',vssrc,' (km/s)'
        WRITE(LOT,*)' Density    src base  : ',densrc,' (gm/cm^3)'
c-----
c       Get receiver region response due to an incidnt P or S wave of
c       unit amplitude
c-----
        WRITE(LOT,*)'-------------------------------------------'
        WRITE(LOT,*)'Computing receiver region surface response',
     1   ' (negative layer index indicates above the source):'
        call getrec(dt,alpha,rayp,npts,recz,recr,rect,'hudsonrec.mod',
     1      hrtel,recdelay,recdist,rectime,
     1      dorec,vprec,vsrec,denrec,dop)
        WRITE(LOT,*)' P-velocity rec base  : ',vprec,' (km/s)'
        WRITE(LOT,*)' S-velocity rec base  : ',vsrec,' (km/s)'
        WRITE(LOT,*)' Density    rec base  : ',denrec,' (gm/cm^3)'

        WRITE(LOT,*)'-------------------------------------------'
        WRITE(LOT,*)'Teleseismic propagation'
c----
c       get the travel time and T* for the teleseismic propagation
c       through the mantle. This is from the base of each layered
c       structure through the mantle. All is OK if the ray parameter
c       is the same
c       Since this give the time from the base of the source and 
c       receiver layers, it give the correct velocity and density
c       of the mantle at this depth for use as halfspace in getrec
c-----
C ORDER IMPORTANT HERE FOR RAYP AND TSTAR
        call firstarr(gcarc*111.195-srcdist-recdist,
     1      hstel,hrtel,dop,ttime,vprec, vsrec, denrec,
     2   vpsrc, vssrc, densrc, rayp, geom, tstar,
     1   .true.)

        if(dop)then
             TPmantle = ttime
        
             WRITE(LOT,*)'Ray parameter        : ',rayp,' (sec/km)',
     1        111.195*rayp,' (sec/deg)'
             WRITE(LOT,*)'P travel time mantle : ',TPmantle,'(s)'
             WRITE(LOT,*)'Geometrical spreading: ',geom
             WRITE(LOT,*)'T*(P) from model     : ',tstar  ,'(s)'
        else 
             TSVmantle = ttime
             WRITE(LOT,*)'Ray parameter        : ',rayp,' (sec/km)',
     1        111.195*rayp,' (sec/deg)'
             WRITE(LOT,*)'S travel time mantle : ',TSVmantle,'(s)'
             WRITE(LOT,*)'Geometrical spreading: ',geom
             WRITE(LOT,*)'T*(S) from model     : ',tstar   ,'(s)'
        endif
        if(utstar.ge.0.0)then
           if(dop)then
             WRITE(LOT,*)'Command line T*(P) used           : ',utstar
           else
             WRITE(LOT,*)'Command line T*(S) used           : ',utstar
           endif
        endif
c-----
c       get the mantle response
c-----

        call tstarpulse(utstar,tstar,tstar,TPmantle,TSVmantle,
     1     geom,geom,dt,alpha,npts,zts,dotel,tstaradj,
     2     dokjar,offset,dop)
        if(dottonly)then
           if(dop)then
             WRITE(LOT,'(a,1x,f5.1,1x,f5.1,1x,f6.2)')
     1          'TT:',gcarc,HS,TPmantle+srctime+rectime
           else
        WRITE(6,*)TSVmantle,srctime,rectime
             WRITE(LOT,'(a,1x,f5.1,1x,f5.1,1x,f6.2)')
     1          'TT:',gcarc,HS,TSVmantle+srctime+rectime
           endif
        else
c-----
c       output synthetics in Green function format for use by hpulse956
c-----
        call outputgreen(npts,dt,alpha,modeltel,modelsrc,
     1  modelrec,srcdelay,srctime,recdelay,rectime,
     1  TPmantle,TSVmantle,ksrc,
     2  tstaradj,recz, recr,rect,zts,srcp,gcarc,offset,dop,hs ,
     4  zdown)
c-----
c       
c-----
        endif
        WRITE(LOT,*)'Program Hudson completed'
        end

        subroutine dosetupsrc(dop,ksrc)
        logical dop
        integer ksrc(21)

             do i=1,21
                 ksrc(i) = 0
             enddo
            if(dop)then
c-----
c                force P for MT sources
c-----
                  ksrc( 1) = 1
                  ksrc( 2) = 1
                  ksrc( 3) = 1
                  ksrc( 4) = 1
                  ksrc( 6) = 1
                  ksrc( 7) = 1
                  ksrc( 9) = 1
                  ksrc(10) = 1
                  ksrc(11) = 1
                  ksrc(12) = 1
                  ksrc(13) = 1
                  ksrc(14) = 1
             else
                  ksrc( 1) = 1
                  ksrc( 2) = 1
                  ksrc( 3) = 1
                  ksrc( 4) = 1
                  ksrc( 5) = 1
                  ksrc( 6) = 1
                  ksrc( 7) = 1
                  ksrc( 8) = 1
                  ksrc( 9) = 1
                  ksrc(10) = 1
                  ksrc(11) = 1
                  ksrc(12) = 1
                  ksrc(13) = 1
                  ksrc(14) = 1
                  ksrc(15) = 1
             endif
           return
           end
        subroutine tstarpulse(utstar,tstarp,tstars,TP,TSV,geomp,geoms,
     1     dt,alpha, npts,zts, dotel,tstaradj, dokjar,offset,dop)
        implicit none
c-----
c       subroutine argumentsa
c-----
        real utstar,tstarp,tstars,TP,TSV,geomp,geoms,dt,alpha
        integer npts
        complex zts(npts)
        logical dotel,dokjar,dop
        real tstaradj, offset
c-----
c       input/output default LUN
c-----
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
c-----
c       utstar  - user provided t* from command line
c                 if == 0  use a very small t* to basically have no t* effect
c                    >  0  use the user provide t* from the command line
c                    <  0  use the t* from the model
c       tstarp  - P wave t* from the teleseismic model
c       tstars  - S wave t* from the teleseismic model
c       TP      - teleseismic P  travel time from base of source model to base of receiver model
c       TSV     - teleseismic SV travel time from base of source model to base of receiver model
c       geomp   - P-wave geometrical spreading from the teleseismic model
c       geoms   - S-wave geometrical spreading from the teleseismic model
c       dt      - sample interval for synthetic 
c       alpha   - time domain damping for avoiding DFT periuodicity
c       npts    - number of samples for final sysnthetic, a power of 2
c       zts     - Fourier transform of the t* pulse
c       dotel   - if true compute t* operator
c       tstaradj- adjustment to get initial pulse. The causal Q oeprator has a group delay
c                 this is an attempt to get the beginning of the t* pulse to have
c                 zero second shift with respect to the travel time. Otherwise the
c                 initial movement will not occur at the ray theorty treavel time
c       dokjar  - .true. use Kjartasson operator, else use  Futterman
c       offset  - source pulse offset - note this is used here only for the telop correction
c                 time offset before the signal. 
c       dop     -   .true. P-save, .false. S-wave
c-----
c-----
c       get the t* attenuation pulse
c-----
        if(dop)then
              if(utstar.eq.0.0)then     
              call telop(0.0001,TP ,geomp,dt,alpha,npts,zts,dotel,
     1             tstaradj,dokjar,offset)
              else if(utstar.gt.0.0)then     
              call telop(utstar,TP ,geomp,dt,alpha,npts,zts,dotel,
     1             tstaradj,dokjar,offset)
              else
              call telop(tstarp,TP ,geomp,dt,alpha,npts,zts,dotel,
     1             tstaradj,dokjar,offset)
              endif
        else
              if(utstar.eq.0.0)then     
              call telop(0.0001,TSV,geoms,dt,alpha,npts,zts,dotel,
     e             tstaradj,dokjar,offset)
C        WRITE(6,*)'nyq2:', nyq2
              else if(utstar.gt.0.0)then     
              call telop(utstar,TSV,geoms,dt,alpha,npts,zts,dotel,
     1             tstaradj,dokjar,offset)
              else
              call telop(tstars,TSV,geoms,dt,alpha,npts,zts,dotel,
     1             tstaradj,dokjar,offset)
              endif
        endif
        WRITE(LOT,*)'T* pulse adjustment  : ',tstaradj,' (s)'
        return
        end
      subroutine outputgreen(npts,dt,alpha,modeltel,modelsrc,
     1  modelrec,srcdelay,srctime,recdelay,rectime,TPmantle,TSVmantle,
     1  ksrc,tstaradj,
     2  recz, recr,rect,zts,srcp,gcarc,offset,dop,hs ,zdown)
       
c-----
c     combine the contributions of the source, receiver and
c     teleseismic propagation to make the complete
c     Green's function
c-----
      implicit none
        integer MAXPTS
        parameter (MAXPTS=16384)
c-----
c     arguments from the subroutine call
c-----
      integer npts
      real dt, alpha,srcdelay,srctime,recdelay,rectime,TPmantle,
     1   TSVmantle,tstaradj
      character modeltel*(*), modelsrc*(*), modelrec*(*)
      integer ksrc(21)
        complex recz(MAXPTS), recr(MAXPTS), rect(MAXPTS)
       complex zts(MAXPTS)
        complex srcp(21,MAXPTS)
       real gcarc, offset
       logical dop, zdown
       real hs
c-----
c       dt      - sample interval for synthetic 
c       alpha   - time domain damping for avoiding DFT periuodicity
c       npts    - number of samples for final sysnthetic, a power of 2
c-----
c-----
c       input/output default LUN
c-----
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
c-----
c     internal variables
c-----
       integer i,iprog, n, n1, n2, nyq2,jk
       real df, fl, fu
       real SR, SL, SN, SA,SC,SF
       real VSA, VSB, VSR
       real VRA, VRB
       real datar, datai, tmp
       real rr, tt0, freq, tshft
       complex ztmp
       real flip
c-----
c       velocity model file
c-----
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep
        real refdepsrc, refdeprec
        real dph

        integer mmax, iunit, iiso, iflsph, idimen, icnvel
        integer lmax
        integer ierr
        character title*80

        real tb, ta, t0
      

c-----
c       output synthetics in Green function format for use by hpulse956
c-----
        open(unit=4,file='hspec96.grn',status='unknown',
     1      form='unformatted',access='sequential')
        rewind 4
c------
c       note I need the A B RHO at source depth and receiver depth
c       not at base of layers as above
c-----
        iprog = 4
        n =npts
        n1 = 1
        n2 = npts/2 + 1
        df = 1./(npts*dt)
        fl = (n1 -1 ) * df
        fu = (n2 -1 ) * df
        nyq2 = 2*(npts/2 + 1)
        write(4)iprog
        write(4) alpha,fl,fu,dt,n,n1,n2,df,nyq2
        write(4)modeltel
c-----
c       safety do an insert
c----
        call getmod(1,modelsrc,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
        if(ierr .lt. 0)stop
        refdepsrc = refdep
        call insert(hs)
            call srclay(hs,lmax,dph)
c-----
c       determine position of source and receive in new layers
c-----

            VSB = b(lmax)
            VSA = a(lmax)
            VSR = rho(lmax)
        if(VSB .le. 0.0001*VSA)then
            do  i=1,21
                if(i.lt.9 .or. i.gt.10)then
                      ksrc(i) = 0
                endif
            enddo
        endif

c-----
c       define TI constants
c-----
            SR = VSR
            SL = VSR*VSB*VSB
            SN = VSR*VSB*VSB
            SA = VSR*VSA*VSA
            SC = VSR*VSA*VSA
            SF = SA - 2.0*SN
        call getmod(1,modelrec,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
        if(ierr .lt. 0)stop
        refdeprec = refdep
            call srclay(0.0,lmax,dph)
            VRB = b(lmax)
            VRA = a(lmax)
c-----
c      if the source in in a fluid, then
c      only output the REX PEX
c      also do not permit S
c-----
c-----
c       TSHFT is the absolute time of the first sample point of the
c       time series. We need this done carefully because of the
c       use of complex frequency (e.g., time damping) for control of
c       singularities. The time of the first sample is
c         -offset + timep + recdelay + srcdelay which accounts for the
c       initial offset, the propagation from the source down to the
c       teleseism model, through the deep earth model, and up from the
c       teleseism model at the receiver to the surface.
c       To do this correctly we must also synthetically delay the
c       teleseism part for propagation. Everything in unraveled in
c       hpulse96
c-----
        
        if(dop)then
            ta = TPmantle + rectime + srctime
            WRITE(LOT,*)'P travel time source layers  :',srctime, '(s)'
            WRITE(LOT,*)'P travel time mantle         :',TPmantle,'(s)'
            WRITE(LOT,*)'P travel time receiver layers:',rectime, '(s)'
            WRITE(LOT,*)'P travel time  total         :',ta,      '(s)'
            TSHFT = -offset + TPmantle + srcdelay + recdelay
            t0 = -12345.
            tb = -offset + ta
        else
            t0 = TSVmantle  + rectime + srctime
            WRITE(LOT,*)'S travel time source layers  :',srctime,  '(s)'
            WRITE(LOT,*)'S travel time mantle         :',TSVmantle,'(s)'
            WRITE(LOT,*)'S travel time receiver layers:',rectime,  '(s)'
            WRITE(LOT,*)'S travel time total          :',t0  ,     '(s)'
            WRITE(LOT,*)'Note P delays in fluid for EX source'
            TSHFT = -offset + TSVmantle + srcdelay + recdelay
            t0 = TSVmantle  + rectime + srctime
            ta = -12345.
            tb = -offset + t0
        endif
C       WRITE(6,*)' First sample (b)  :',tB , ' s'
C       WRITE(6,*)' Origin time  (o)  :',0.0, ' s'
C       WRITE(6,*)' P time       (a)  :',ta , ' s'
C       WRITE(6,*)' S time       (t0) :',t0 , ' s'
            write(4)gcarc*111.195,tb,hs-refdepsrc,
     1          0.0-refdeprec,
     2          ta,t0,t0,
     3          SA, SC, SF, SL, SN, SR

c-----
c       shift the time series to account for T* adjustment
c       and propagation through the source and receiver layers
c-----
        TSHFT = TSHFT + tstaradj 
        write(4)ksrc
        do i=n1,n2
           freq = (i-1)*df
c-----
c          handle the first point as origin time
c-----
           tmp = 6.2831853*freq*(TSHFT )
           ztmp = cmplx(cos(tmp),sin(tmp))
       
           do jk=1,21
c-----
c             Check if this Green is to be output
c-----
              if(ksrc(jk).eq.1)then
c-----
c                  vertical compoent
c-----
                   if(jk.eq.1 .or. jk.eq.3 .or. jk.eq.6
     1                        .or. jk.eq.9 .or. jk.eq.11 
     2                        .or. jk.eq.13)then
c-----
c                     At this point Uz is positive down. Invert
c                     if positive up is desired
c-----
                      if(zdown)then
                        flip = 1.0
                      else 
                        flip = -1
                      endif
                      datar = flip*real (recz(i)*ztmp*zts(i)*srcp(jk,i))
                      datai = flip*aimag(recz(i)*ztmp*zts(i)*srcp(jk,i))
c-----
c                  radial component
c-----
                   else if(jk.eq.2 .or. jk.eq.4 .or. jk.eq.7
     1                        .or. jk.eq.10 .or. jk.eq.12 
     2                        .or. jk.eq.14)then
                      datar = real (recr(i)*ztmp*zts(i)*srcp(jk,i))
                      datai = aimag(recr(i)*ztmp*zts(i)*srcp(jk,i))
c-----
c                  transverse component
c-----
                   else if(jk.eq.5 .or. jk.eq.8 .or. jk.eq.15)then  
                      datar = real (rect(i)*ztmp*zts(i)*srcp(jk,i))
                      datai = aimag(rect(i)*ztmp*zts(i)*srcp(jk,i))
                   endif
                   write(4)datar,datai
              endif
           enddo
        enddo
        rr = -1.0
        tt0 = 0.0
        write(4)rr,tt0

        close(4)
        return
        end
        subroutine zfour(zarr,nn,isign,dt,df) 
c-----
c     THE input is a complex array
c     which has numbers stored in memory as
c     R1, I1, R2, I2, ..., Rnn, Inn
c     where nn must be a power of 2 R and I are the real and imaginary
c     parts of the complex number
c
c     For isign -1 this is a complex time series
c     For isign +1 this is a complex frequency series with
c        index 1 (in fortran corresponding to f=0
c              2                              f=df
c            nn/2 + 1                         f = 1/2dt = Nyquist
c            nn - 1                           f = -2df
c             nn                              f = -df

c-----
c     the cooley-tookey fast fourier transform in usasi basic fortran
c     transform(j) = sum(datc(i)*w**((i-1)(j-1)), where i and j run
c     from 1 to nn and w = exp(isign*2*pi*sqrt(-1)/nn).  datc is a one-
c     dimensional complex array (i.e., the real and imaginary parts of
c     datc are located immediately adjacent in storage, such as fortran
c     places them) whose length nn is a power of two.  isign
c     is +1 or -1, giving the sign of the transform.  transform values
c     are returned in array datc, replacing the input datc.  the time is
c     proportional to n*log2(n), rather than the usual n**2
c     rms resolution error being bounded by 6*sqrt(i)*log2(nn)*2**(-b),
c     b is the number of bits in the floating point fraction.
c
c     the program computes df from dt, dt from df and checks to see
c     if they are consistent. In addition, the transforms are multiplied
c     by dt or df to make the results dimensionally correct
c
c     This is a slightly modified version of the original Brenner routine
c     The changes were to add physical dimensions to the transform
c     and to make it all complex
c-----
        integer nn, isign
        complex zarr(nn) 
        real dt, df

        complex ztemp
c-----
c       ensure that the dt and df are defined and
c       consistent
c-----
        if(dt.eq.0.0) dt = 1./(nn*df) 
        if(df.eq.0.0) df = 1./(nn*dt) 
        if(dt.ne.(nn*df)) df = 1./(nn*dt) 
c-----
c       now begin the transform
c-----
        jj = 1
        do 5 ii=1,nn 
        if(ii .lt. jj) then
              ztemp = zarr(jj)
              zarr(jj) = zarr(ii)
              zarr(ii) = ztemp
        endif
        mm = nn/2
    3   continue
        if(jj.le.mm) then
            go to 55
        else 
              jj = jj - mm
              mm = mm /2
              if(mm.lt.1)then
                  go to 55
              else
                  go to 3
              endif
        endif
   55   continue
        jj = jj + mm
    5   continue
        mmax = 1 
c-----
    6 continue
        if(mmax .lt. nn)then
            go to 7
        else if(mmax .ge. nn)then
            go to 10
        endif
    7   continue
        istep= 2*mmax 
        theta = 6.283185307/(isign*2.0*mmax) 
        sinth=sin(theta/2.) 
        wstpr=-2.*sinth*sinth 
        wstpi=sin(theta) 
        wr=1.0 
        wi=0.0 
        do 9 mm=1,mmax
              do 8 ii=mm,nn,istep
                    jj=ii+mmax
                    ztemp=cmplx(wr,wi)*zarr(jj)
                    zarr(jj) = zarr(ii) - ztemp
                    zarr(ii) = zarr(ii) + ztemp
    8         continue
c-----
c       use trig relations to compute the next sin/cos
c       without actually calling sin() or cos()
c-----
              tempr = wr 
              wr = wr*wstpr-wi*wstpi + wr 
              wi = wi*wstpr+tempr*wstpi + wi 
    9   continue
        mmax = istep 
        go to 6 
c-----
c       transform is done
c-----
   10   continue 
c-----
c     give the arrays the proper physical dimensions
c-----
        if(isign.lt.0)then
c-----
c             time to frequency domain
c-----
              do  ii = 1,nn
                    zarr(ii) = zarr(ii) * dt
              enddo
        else
c-----
c             frequency to time domain
c-----
              do ii = 1,nn
                    zarr(ii) = zarr(ii) * df
              enddo
        endif
        return
        end
        subroutine getindex(dlow,dhgh,d,mmax,k)
c-----
c       convert the thickness array d into depths and
c       then find the bounding index for depths between dlow and dhgh
c-----
        implicit none
        real dlow, dhgh, d(*)
        integer mmax, k

        integer i
        real dphl,dphh
        real dmid
        dmid = 0.5*(dlow+dhgh)
        dphl = 0.0
        do i=1,mmax 
          dphh = dphl + d(i)
          if(dmid.ge.dphl .and. dmid.le.dphh)then
              k = i
              return
          endif
          dphl = dphh
        enddo
        k = -1
        return
        end
       

        subroutine causlq(freq,alpha,qi,tp,ztmp,dokjar)
c-----
c       Kjartansson, E. (1979).
c          Constant Q-wave propagation and attenuation,
c       J. Geophys. Res. 84, 4737-4748.
c-----
            real freq, tp,  qi
            complex Ztmp
            logical dokjar

            real gam, fac, om1, fref
            real pi, pi2
            complex atn, omega
            double complex zatn, zfac
            real zr, zi

            real ts

            ts = tp * qi
            omega = cmplx(6.2831853*freq, -alpha)
            fref = 1.0
            om1 = 6.2831853*fref
            pi = 3.1415927
            pi2 = pi / 2.

            

            if(dokjar)then

                   gam = atan(qi)/pi
                   if(gam.le.0.0)then
                      zatn = cmplx(1.0,0.0)
                   else
                       fac = pi2*gam
                       rfac = sin(fac)/cos(fac)
                       zatn =
     1                     (omega/om1)**dble(1.0-gam)  *
     2                     dcmplx(dble(rfac),1.0d+00)
                   endif
                   atn=cmplx(sngl(real(zatn)),sngl(dimag(zatn)))
c-----
c           form the propagation term
c-----
                   
                   zfac = - om1*tp*atn 
            else
                 zfac =  omega*ts*clog(omega/om1)/pi + 
     1             cmplx(0.0,1.0)*omega*ts/2
                 zfac = zfac * cmplx(0.0,1.0) 
            endif
            zr = sngl(dreal(zfac))
            zi = sngl(dimag(zfac))
            if(zr.gt. -88)then
                 ztmp = exp(zr)*cmplx(cos(zi),sin(zi))
            else
                 ztmp = cmplx(0.0,0.0)
            endif
        return
        end

        subroutine adomod()
c-----
c       just fill the rhosh, bsh and qbsh arrays 
c-----
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        do  i=1,mmax
            bsh(i)=b(i)
            qbsh(i)=qb(i)
            rhosh(i) = rho(i) 
        enddo
        return
        end

        subroutine adosph()
c-----
c       Transform spherical earth to flat earth
c
c       Schwab, F. A., and L. Knopoff (1972). 
c           Fast surface wave and free
c       mode computations, in  
c           Methods in Computational Physics, Volume 11,
c       Seismology: Surface Waves and Earth Oscillations,  
c           B. A. Bolt (ed),
c       Academic Press, New York
c
c       Love Wave Equations  44, 45 , 41 pp 112-113
c       Rayleigh Wave Equations 102, 108, 109 pp 142, 144
c
c
c-----
c       mmax    I*4 number of layers
c       ipsvsh  I*4     1 - get P time
c                       2 - get SV time
c                       3 - get SH time
c-----
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        double precision z0,z1,r0,r1,dr,ar,tmp

        common/earth/radius
        real radius

        ar=radius
        dr=0.0d0
        r0=ar + refdep
        d(mmax)=1.0
        do 10 i=1,mmax
            r1=r0-dble(d(i))
            z0=ar*dlog(ar/r0)
            z1=ar*dlog(ar/r1)
            d(i)=sngl(z1-z0)
c-----
c        attempt 7 15 2007 - use standard rule but at mid layer depth as per DGH
c-----
            TMP=(ar+ar)/(r0+r1)

            a(i)=a(i)*sngl(tmp)
            btmp = b(i)
            b(i)=sngl(btmp*tmp)
            bsh(i)=sngl(btmp*tmp)
            qbtmp = qb(i)
            qbsh(i)=qbtmp
            rhotmp=rho(i)
            rhosh(i) = sngl(rhotmp * tmp **(-5.0))
            rho(i) = sngl(rhotmp * tmp **(-2.275))
            r0 = r1
   10   continue
        d(mmax)=0.0
        return
        end
        subroutine uniq(n,x)
        integer n
        real x(n)
      
        m = 1
        do i=2,n
              if(x(i).ne.x(m))then
                  m = m + 1
                  x(m) = x(i)
              endif
        enddo
        n = m
        return
        end

        subroutine werror(ostr)
c-----
c       output error message and terminate program
c-----
        parameter(LER=0, LIN=5, LOT=6)
        character ostr*(*)
        write(LER,*)'PROGRAM TERMINATION'
        write(LER,*)ostr
        stop
        end

        subroutine dezero()
        parameter (LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
c-----
c       ultimately get rid of zero thickness layers - this
c       will require readjusting the model from top down, and
c       also readjusting the source and receiver indices.
c----
c       Here just guarantee that the halfspace is not of zero thickness
c-----
        dmin = 1.0e+30
        do 100 i=1,mmax-1
            if(d(i) .lt. dmin .and. d(i).gt.0.0)dmin = d(i)
  100   continue
c       if(d(mmax).le.0.0)then
c           d(mmax) = 0.1*dmin
c       endif
        return
        end

        subroutine chkmod()
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/jout/jsrc(21) , jbdrys, jbdryh
        common/lwater/lfluid
        logical lfluid
c-----
c       check model for inconsistencies
c-----
c-----
c       Model cannot consist entirely of water layers
c       Also determine first solid layer from top
c-----
        iw = 0  
        isoldd = 0
        do 100 i=1,mmax
            if(b(i).eq.0.0)then
                iw = iw + 1
            else
                if(isoldd .eq.0)isoldd=i
            endif
  100   continue
        if(iw .eq. mmax)then
            lfluid = .true.
C           call werror('MODEL CONSISTS ONLY OF LIQUID LAYERS')
        else
            lfluid = .false.
        endif
c-----
c       Determine first solid layer from bottom
c-----
        iw = 0  
        isoldu = 0
        do 101 i=mmax,1,-1
            if(b(i).eq.0.0)then
                iw = iw + 1
            else
                if(isoldu .eq.0)isoldu=i
            endif
  101   continue
c-----
c       Check for interior water layer
c-----
        if(iw.gt.0 .and. .not. lfluid)then
            do 102 i=isoldd,isoldu
                if(b(i).eq.0.0)then
                call werror('MODEL HAS INTERIOR  FLUID LAYERS')
                endif
  102       continue
        endif
c-----
c       If boundary condition is rigid, and the adjacent layer is
c       fluid, terminate 
c-----
C       if(b(1).le.1.0e-04 .and. jbdrys.eq.-1 .and. lfluid)then
C           call werror('TOP LAYER IS FLUID AND RIGID BOUNDARY')
C       endif
C       if(b(mmax).le.1.0e-04 .and. jbdryh.eq.-1 .and. lfluid)then
C           call werror('BOTTOM LAYER IS FLUID AND RIGID BOUNDARY')
C       endif
        return
        end

        subroutine insert(dph)
        implicit none
        real dph
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        real dp, dphh, hsave, dep
        integer m, ls
c-----
c       Insert a depth point into the model by splitting a layer
c       so that the point appears at the top boundary of the layer
c       dph = depth of interest
c-----
c       determine the layer in which the depth dph lies.
c       if necessary, adjust  layer thickness at the base
c-----
c       Here determine layer corresponding to specific depth dph
c       If the bottom layer is not thick enough, extend it
c
c       dep - depth to bottom of layer
c       dphh    - height of specific depth above bottom of the layer
c-----
        ls = 1
        dphh =  0
        if(dph.le.0)then
            d(1) = d(1) - dph
            return
        else if(dph.ge.0)then
            dep = 0.0 
            dp = 0.0 
            dphh = -1.0
            do 100 m = 1,mmax 
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
  100       continue 
  101       continue 
        endif
c-----
c       In the current model, the depth point is in the ls layer
c       with a distance dphh to the bottom of the layer
c
c       Do not create unnecessary layers, e.g., 
c           at surface and internally
c       However do put in a zero thickness layer 
c           at the base if necessary
c-----
        if(dphh .eq. 0.0 .and. ls.ne.mmax)then
            return
        else
c-----
c           adjust layering
c-----
             do 102 m = mmax,ls,-1
                d(m+1) = d(m)
                a(m+1) = a(m)
                b(m+1) = b(m)
                rho(m+1) = rho(m)
                qa(m+1) = qa(m)
                qb(m+1) = qb(m)
                etap(m+1) = etap(m)
                etas(m+1) = etas(m)
                frefp(m+1) = frefp(m)
                frefs(m+1) = frefs(m)
                bsh(m+1) = bsh(m)
                qbsh(m+1) = qbsh(m)
                rhosh(m+1) = rhosh(m)
  102       continue
            hsave = d(ls)
            d(ls) = hsave - dphh
            d(ls+1) = dphh
            ls = ls + 1
            mmax = mmax + 1
            if(d(mmax).lt.0.0)d(mmax)=0.0
        endif
        return
        end

        subroutine telop(tstar,T,geom,dt,alpha,npts,zts,dotel,
     1     tstaradj,dokjar,offset)
c-----
c       compute zts, the Fourier transform of the teleseism
c       response to include the time shift due to propagation and
c       the anelastic attenuation operation.
c
c       We follow 
c       Kjartansson, E., Constant Q-wave propagation and
c             Attenuation, J. Geophys. Res., Vol 84, 4737-4748, 1979
c
c       tstar     R*4    - T*
c       T         R*4    - travel time, required for Kjartansson operator
c       geom      R*4    - teleseismic geometrical spreading factor
c       dt        R*4    - samling interval
c       alpha     R*4    - time domain attenuation operator
c       npts      I*4    - number of data points
c       zts       C*4    - array of complex attenuation operators
c       dotel     L      - .true. comp[ute complex T* operator
c       tstaradj  R*4    - adjustment to get initial pulse
c       dokjar    L      - .true. use Kjartasson operator, else use
c                          Futterman
c       offset    R      - source pulse offset - note this is
c                          use here only for the telop correction
c----- 
        implicit none

        integer MAXPTS
        parameter (MAXPTS=16384)
        real tstar, T, dt, alpha, geom, tstaradj
        integer npts
        complex zts(MAXPTS)
        logical dotel, dokjar
        real offset 
c-----
c       internal variables
c-----
        integer n21, i
        real freq, df
        complex ztmp
        real qi
        real fac, dfac, depmax
        integer ifound

        complex z(MAXPTS)
        real x(MAXPTS)
        real x0, y0, x1, y1, xi


        n21 = npts/2 + 1
        df = 1./(npts*dt)
        if(tstar.gt.0.0)then
        qi = tstar/T
        else
        qi = 0.0
        endif


        if(.not. dotel)then
c-----
c       turn off for testing
c-----
             do i=1,n21
                  zts(i) = cmplx(1.0,0.0)
             enddo
             T = 0.0
             tstaradj = 0.0
        else
              do i=1,n21
                 freq = (i-1)*df
                 if(freq.lt.df)freq = 0.01 * df
                 call causlq(freq,alpha,qi,T,ztmp,dokjar)
                 zts(i) = ztmp*geom
              enddo
c-----
c       compute the time delay
c           define a zero phase source pulse
c-----
            do i=1,npts
                z(i) = cmplx(0.0,0.0)
                if(i.eq.1)then
                    z(i) = cmplx(0.50/dt,0.0)
                else if(i.eq.2)then
                    z(i) = cmplx(0.25/dt,0.0) * exp(-alpha*dt)
                else if(i.eq.npts)then
                    z(i) = cmplx(0.25/dt,0.0) * exp(+alpha*dt)
                endif
            enddo
            call zfour(z,npts,-1,dt,df)

              do i=1,n21
                 freq = (i-1)*df
                 if(freq.lt.df)freq = 0.01 * df
                 call causlq(freq,alpha,qi,T,ztmp,dokjar)
                 z(i) = z(i) * ztmp
              enddo
c-----
c            apply the time offset
c-----
             do i=1,n21
                 freq = (i-1)*df
                 fac = 6.2831853*freq*(-offset+T)
                 ztmp = cmplx(cos(fac), sin(fac) )
                 z(i) = z(i) * ztmp
             enddo
c-----
c            force a real time series
c-----
             do i=1,n21
                  if(i.gt.1)then
                      z(npts -i + 2) = conjg(z(i))
                  endif
             enddo
             call zfour(z,npts,+1,dt,df)
c-----
c            undamp
c-----
             fac = exp(alpha*offset)
             dfac = exp(alpha*dt)
             depmax = -1.0e+37
             do i=1,npts
                   z(i) = z(i) * fac
                   fac = fac * dfac
                   x(i) = real(z(i))
                   if(x(i).gt.depmax)then
                       depmax = x(i)
                   endif
             enddo
c-----
c            attempt to get the initial time by finding the first point
c                   >= depmax/100.0
c-----
            tstaradj = 0
            ifound = 0
            tstaradj = 0.0
            do i=1,npts-1
               if(x(i+1).ge.0.01*depmax .and. ifound.eq.0)then
                     ifound = 1
c-----
c                    do a 2 point interpolation
c-----
                     x0 = i
                     y0 = x(i)
                     x1 = i+1
                     y1 = x(i+1)
                     if(y1 .ne. y0)then
                     xi = x0 + ( 0.01*depmax - y0)*( x1-x0)/(y1-y0)
                     tstaradj = - offset + (xi-1)*dt 
                     endif
               endif
            enddo
        endif

        
        return
        end

        subroutine bsort(nn,x,isign)
c-----
c   http://linux.wku.edu/~lamonml/algor/sort/bubble.html
c   transliterated to FORTRAN
c-----
c       do bubble sort.
c       isign=+1  increase   =-1 decrease.
c
        real x(nn)
        real temp
        integer i,j
        do i= nn ,1, -1
             do j=1,i-1,1
             x0 = x(j+1)-x(j)
             if(isign.gt.0)x0 = - x0
             if( x0 .gt.0.0)then
                   temp = x(j)
                   x(j) = x(j+1)
                   x(j+1) = temp
             endif
             enddo
        enddo
      
        return
        end


        subroutine usage(ostr)
c------
c       write out program syntax
c-----
        character ostr*(*)
        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)
        if(ostr.ne. ' ' )then
            lostr = lgstr(ostr)
            write(LER,*)ostr(1:lostr)
        endif
        write(LER,*)'USAGE: ',
     1  'hudson96 [-TEL modeltel] [-SRC modelsrc] [-REC modlrec] ',
     2  '[-HS hs ]'
        WRITE(LER,*)
     1  '  [ -O offset ] [-P] [-S] [-GCARC gcarc] [-DT dt] [-NPTS npts]'
        WRITE(LER,*)
     3  '  [-NOSRC] [-NOREC] [-NOTEL] [-TSTAR tstar] [-ZDOWN]'
        WRITE(LER,*)
     4  '  [-TTONLY] [-HSRC hsrc] [-HREC hrec] [-?] [-h] [-V|Q]'
        write(LER,*)
     1  '-TEL modeltel (default none) Earth velocity model'
        write(LER,*)
     1  '-SRC modelsrc (default modeltel) Velocity model for source'
        write(LER,*)
     1  '-REC modelrec (default modeltel) Velocity model for receiver'
        write(LER,*)
     1  '-HS hs        (default 10 km) Source depth '
        write(LER,*)
     1  '-GCARC gcarc  (default 30 degrees) Arc distance '
        write(LER,*)
     1  '-NPTS npts    (default 1024) Number of points in time series'
        write(LER,*)
     1  '-DT dt        (default 0.05 sec) Sampling interval'
        write(LER,*)
     1  '-O  offset    (default 10.0 s for P, 20 s for S) Time offset ',
     1  'before signal '
        write(LER,*)
     1  '-P            (default true)  make P synthetic '
        write(LER,*)
     1  '-S            (default false)  make S synthetic '
        write(LER,*)
     1  '-NOSRC        (default false) Turn off source input'
        write(LER,*)
     1  '-NOREC        (default false) Turn off receiver part'
        write(LER,*)
     1  '-F            (default false) Use Futterman operator'
        write(LER,*)
     1  '                     else use Kjartansson (1979)'
        write(LER,*)
     1  '-NOTEL        (default false) Turn off teleseism ',
     2  ' geometrical'
        write(LER,*)
     1  '-TSTAR tstar  (default -12345.) If >=0 use this T* instead'
        write(LER,*)
     2  '               of one computed from velocity model'
        write(LER,*)
     1  '-TTONLY       (default false) Output only travel times'
        write(LER,*)
     1  '                     spreading and attenuation'
        write(LER,*)
     1  '-ZSRC zsrc    (default see below',
     2     ' Maximum thickness of Src model merged with'
        write(LER,*)
     1  '                     teleseismic model '
        write(LER,*)
     1  '-ZREC zrec    (default see below',
     2     ' Maximum thickness of Rec model merged with'
        write(LER,*)
     1  '                     teleseismic model '
        write(LER,*) 
     1 '              If halfspace of Src/Rec model has thickness',
     1     ' of 0 km, ' 
        write(LER,*) 
     1 '              the halfspace is extended to the depth of',
     2     ' zsrc/zrec km. '
        write(LER,*)
     1 '              The default of 100 km for the source and ',
     3     ' 60 km '
        write(LER,*)
     1 '              for the receiver models is adequate for crustal',
     2 ' models' 
        write(LER,*)
     1  '-ZDOWN        (default false) Uz is positive down. The',
     2  ' default is positive up,'
        write(LER,*)
     1  '                     the convention of earthquake seismology'
        write(LER,*)
     1  '-?            Display this usage message'
        write(LER,*)
     1  '-h            Display this usage message'
        write(LER,*)
     1  '-V            (default true ) verbose output'
        write(LER,*)
     1  '-Q            (default flase ) quiet output'
        write(LER,*)
     1  '-----------------------------------------------'
        write(LER,*)
     1  'The program creates the file hspec96.grn for hpulse96.'
        write(LER,*)
     1  'The program creates the files hudsonsrc.mod and',
     2  ' hudsonrec.mod,'
        write(LER,*)
     1  'which are the complete models (model96 format) used ',
     2  'for the src and rec regions from the surface to interior.'
        write(LER,*)
     1  'these can be plotted using shwmod96.'
        write(LER,*)


        stop
        end


        subroutine srclyr(depth,lmax,dph)
        implicit none
        real depth, dph
        integer lmax
        integer LER, LIN, LOT
        parameter (LER=0, LIN=5, LOT=6)
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        integer m
        real dep
c-----
c       Find source/receiver boundary. It is assumed that
c       it will lie upon a boundary
c
c       lmax = source layer 
c            = 0 is the free surface 
c       depth = source depth 
c       dph = height of  source above lmax + 1 interface 
c       lmax = 0 is the free surface 
c-----
        if(depth.le.0.0)then
            lmax = 1
            dph = 0.0
        else
            dep = 0.0 
            do 100 m = 2,mmax
                dep = dep + d(m-1) 
                dph = dep - depth 
                lmax = m 
                if(abs(dph).lt. 0.0001*d(m-1) .or.
     1              abs(dph).lt.1.0e-6)go to 101
  100       continue 
  101   continue 
        endif
        return 
        end 

        subroutine modcpy(totmp) 
        logical totmp
c-----
c       copy model to temporary array
c-----
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/modelt/dt(NL),at(NL),bt(NL),rhot(NL),
     1      mmaxt,qat(NL),qbt(NL),etapt(NL),etast(NL),
     1      frefpt(NL),frefst(NL)
c-----
c       copy to temporary array
c-----
        if(totmp)then
            do 20 i = 1,mmax 
                dt(i) = d(i)
                at(i) = a(i)
                bt(i) = b(i)
                rhot(i) = rho(i)
                qat(i) = qa(i)
                qbt(i) = qb(i)
                etapt(i) = etap(i)
                etast(i) = etas(i)
                frefpt(i) = frefp(i)
                frefst(i) = frefs(i)
   20       continue 
            mmaxt = mmax
        else
            do 30 i = 1,mmaxt 
                d(i) = dt(i)
                a(i) = at(i)
                b(i) = bt(i)
                rho(i) = rhot(i)
                qa(i) = qat(i)
                qb(i) = qbt(i)
                etap(i) = etapt(i)
                etas(i) = etast(i)
                frefp(i) = frefpt(i)
                frefs(i) = frefst(i)
   30       continue 
            mmax = mmaxt
        endif
        return 
        end 

        subroutine gcmdln(dt, gcarc, hs, npts,
     1    modeltel,modelsrc,modelrec,offset,
     2    dosrc,dorec,dotel,dop,dokjar,utstar,dottonly,zsrc,zrec,
     3    verbose,zdown)
c-----
c       parse command line arguments
c       requires subroutine mgtarg() and function mnmarg()
c-----
        implicit none

        integer LER,LOT,LIN
        parameter(LER=0, LIN=5, LOT=6)

        character modeltel*(*)
        character modelsrc*(*)
        character modelrec*(*)

        real dt, gcarc, hs, offset,utstar,zsrc,zrec
        integer npts
        logical dosrc,dorec,dotel,dop,dokjar,dottonly,verbose
        logical zdown

        integer lnblnk

c-----
c       internal variables
c-----
        integer i, nmarg
        integer mnmarg
        character*100 name
        logical ext

        logical loffsetdefault

c-----
c       initialize
c-----
        modeltel = ' '
        modelsrc = ' '
        modelrec = ' '

        hs = 0.0
        dt = 1
        npts = 1024
        gcarc = 50
        offset = 10
        dosrc = .true.
        dorec = .true.
        dotel = .true.
        dop = .true.
        dokjar = .true.
        loffsetdefault = .true.
        utstar = -12345.
        dottonly = .false.
        zsrc = 100.0
        zrec =  60.0
        verbose = .true.
        zdown = .false.
   

        nmarg=mnmarg()
        i = 0
   11   i = i + 1
             if(i.gt.nmarg)goto 13
             call mgtarg(i,name)
             if(name(1:4).eq.'-TEL')then
                   i = i + 1
                   call mgtarg(i,modeltel)
             else if(name(1:2).eq.'-S' .and.
     1             name(1:4).ne.'-SRC')then
                   dop = .false.
             else if(name(1:4).eq.'-SRC')then
                   i = i + 1
                   call mgtarg(i,modelsrc)
             else if(name(1:2).eq.'-P')then
                   dop = .true.
             else if(name(1:4).eq.'-REC')then
                   i = i + 1
                   call mgtarg(i,modelrec)
             else if(name(1:3).eq.'-DT')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,f10.0)')dt
             else if(name(1:3).eq.'-HS')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,f10.0)')hs
             else if(name(1:2).eq.'-G')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,f10.0)')gcarc
             else if(name(1:6).eq.'-TSTAR')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,f10.0)')utstar
             else if(name(1:5).eq.'-ZSRC')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,f10.0)')zsrc
             else if(name(1:5).eq.'-ZREC')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,f10.0)')zrec
             else if(name(1:6).eq.'-ZDOWN')then
                   zdown = .true.
             else if(name(1:7).eq.'-TTONLY')then
                   dottonly = .true.
             else if(name(1:2).eq.'-O')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,f10.0)')offset
                   if(offset.le.0.0)then
                         offset = 10.0
                         loffsetdefault = .false.
                   else
                         loffsetdefault = .true.
                   endif
             else if(name(1:2).eq.'-N' .and.
     1      name(1:3).ne.'-NO')then
                   i = i + 1
                   call mgtarg(i,name)
                   read(name,'(bn,i10)')npts
             else if(name(1:6).eq.'-NOSRC' )then
                   dosrc = .false.
             else if(name(1:6).eq.'-NOREC' )then
                   dorec = .false.
             else if(name(1:6).eq.'-NOTEL' )then
                   dotel = .false.
             else if(name(1:2).eq.'-F')then
                   dokjar = .false.
             else if(name(1:2).eq.'-?')then
                call usage(' ')
             else if(name(1:2).eq.'-h')then
                call usage(' ')
             else if(name(1:2).eq.'-V')then
                verbose = .true.
             else if(name(1:2).eq.'-Q')then
                verbose = .false.
             endif

             go to 11
   13   continue
c-----
c       safety checks
c-----
        if(modeltel .eq. ' ')call usage(' ')
        if(modelsrc .eq. ' ')modelsrc = modeltel
        if(modelrec .eq. ' ')modelrec = modeltel
        if (zsrc.le.0.0)zsrc = 100.0
        if (zrec.le.0.0)zrec =  60.0
        inquire(file = modeltel, exist=ext)
        if( .not. ext)then
            write(LER,*)' Model file is not valid: ',
     1      modelrec(1:lnblnk(modeltel))
            call usage(' ')
        endif
        inquire(file = modelsrc, exist=ext)
        if( .not. ext)then
            write(LER,*)' Source region file is not valid: ',
     1      modelrec(1:lnblnk(modeltel))
            call usage(' ')
        endif
        inquire(file = modelrec, exist=ext)
        if( .not. ext)then
            write(LER,*)' Receiver region file is not valid: ',
     1      modelrec(1:lnblnk(modeltel))
            call usage(' ')
        endif
         

        if(hs.lt.0.0 .or. hs.gt.750)then
            call usage(' Only depth between 0 - 750 km permitted')
        endif
        if( loffsetdefault .and. .not. dop )then
            offset = 20.0
        endif
        
        return
        end

        subroutine getrec(dt,alp,rayp,npts,recz,recr,rect,
     1      modelrec,hrtel,recdelay,recdist,rectime,dorec,
     2      pvelrec,svelrec,densrec,pincident)
c-----
c       get surface displacements in receiver crust due tro a
c       unit amplitude incident P or S plane wave in the halfspace
c-----
        implicit none
c-----
c       routine calling arguments
c-----
        real dt,alp,rayp,pvelrec,svelrec,densrec
        integer npts
        complex recz(*), recr(*), rect(*)
        character modelrec*(*)
        real hrtel, recdelay
        real recdist,rectime
        logical dorec
        logical pincident

c-----
c
c-----
        common/damp/alpha,ieqex
            real alpha
            integer ieqex
c-----
c       ray parameter values
c-----
        common/c/pmin,pmax,dp,pcntrl
        real pmin,pmax,dp,pcntrl
c-----
c       receiver model parameters
c-----
        integer LER, LIN, LOT
        parameter(LER=0, LIN=5, LOT=6)
c-----
c       earth model information
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep

        integer iunit, iiso, iflsph, idimen, icnvel, ierr
        character title*80
        common/lyrctl/lyrins
        logical lyrins
c-----
c       internal variables
c-----
        integer i, n21
        real freq, df
        complex Z, R, T
        integer lmaxr
        real dphr
        real eta
 
        alpha = alp
        pmin = rayp
        pmax = rayp
        dp = rayp
        pcntrl = -1.0
        lyrins = .true.

        n21 = npts/2 + 1
        df = 1.0/(npts*dt)
c-----
c       get the model at the receiver
c-----
        call getmod(1,modelrec,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
        if(ierr .lt. 0)stop
c-----
c       make sure that we use 1/Q
c-----
        do i=1,mmax
            if(qa(i).lt.0.0)qa(i) = 0.0
            if(qb(i).lt.0.0)qb(i) = 0.0
            if(qa(i) .gt. 1.0)qa(i) = 1.0/qa(i)
            if(qb(i) .gt. 1.0)qb(i) = 1.0/qb(i)
            if(frefp(i) .le. 0.0)frefp(i) = 1.0
            if(frefs(i) .le. 0.0)frefs(i) = 1.0
        enddo
c-----
c       insert a boundary at hrtel and determine the index of the layer here
c       output for safety as a debug
c-----
C        WRITE(6,*)'hrtel,refdep',hrtel,refdep
        call insert(hrtel+refdep)
        call dezero()
        call srclyr(hrtel+refdep,lmaxr,dphr)
c-----
c       Change the maximum number of layers and
c       place teleseism model parameters into new base
c       layer - note in future preserve the Q etc
c-----
        mmax = lmaxr
        a(mmax)   = pvelrec
        b(mmax)   = svelrec
        rho(mmax) = densrec
c-----
c       should also carry the Q etc
c-----
        WRITE(LOT,'(a)')
     1 '  LAYER   H(km)      PVel      SVel      Dens     QPi       QSi'
        do i=1,lmaxr
              if(i.lt.lmaxr)then
              WRITE(LOT,'(i5,4f10.2,2f10.4)')i,d(i),a(i),b(i),rho(i)
     1            ,qa(i),qb(i)
              else
              WRITE(LOT,'(i5,4f10.2,2f10.4)')-i,d(i),a(i),b(i),rho(i)
     1            ,qa(i),qb(i)
              endif
        enddo
c-----
c       compute the vertical tau delay
c-----
        recdelay = 0.0
        recdist  = 0.0
        rectime  = 0.0
        do i=1,lmaxr-1
           if(pincident)then
                eta = sqrt(1.0 - rayp*rayp*a(i)*a(i))
                recdelay = recdelay + d(i)*eta/a(i)
                recdist  = recdist  + d(i)*rayp*a(i)/eta
                rectime  = rectime  + d(i)/(a(i)*eta)
           else
                eta = sqrt(1.0 - rayp*rayp*b(i)*b(i))
                recdelay = recdelay + d(i)*eta/b(i)
                recdist  = recdist  + d(i)*rayp*b(i)/eta
                rectime  = rectime  + d(i)/(b(i)*eta)
           endif
        enddo
        WRITE(LOT,*)' recdelay  :',recdelay, 
     1     ' s  (used for receiver synthetic time shift)'
        WRITE(LOT,*)' recdist   :',recdist , 
     1     ' km (used for mantle path distance correction)'
        WRITE(LOT,*)' rectime   :',rectime , 
     1     ' s  (used for travel time computation)' 
    

c-----
c       compute the response renaming the excit of hrftn96 to exitpw
c       for incident P
c-----
        do i=1,n21
             freq = (i-1)*df
             if(freq.lt.df) freq = 0.01*df
                 call excitpw(freq,pincident,rayp,Z,R,T,pvelrec,svelrec)
             recz(i) = Z
             recr(i) = R
             rect(i) = T
        enddo
c-----
c       turn off for testing
c-----
        if(.not. dorec)then
             recdelay = 0.0
             rectime = 0.0
             do i=1,n21
                   recz(i) = cmplx(1.0,0.0)
                   recr(i) = cmplx(1.0,0.0)
                   rect(i) = cmplx(1.0,0.0)
             enddo
        endif
        return
        end

        subroutine getsrc(dt,alp,rayp,npts,srcp,modelsrc,hs,hstel,
     1      srcdelay,srcdist,srctime,dosrc,vsa,vsb,vsr,pteleseismic)
        implicit none
c-----
c       routine calling arguments
c-----
        real dt,alp,rayp,vsa,vsb,vsr
        integer npts
        character modelsrc*(*)
        real hs, hstel, srcdelay
        real srcdist, srctime
        logical dosrc
        logical pteleseismic

        integer MAXPTS
        parameter (MAXPTS=16384)
        complex srcp(21,MAXPTS)

c-----
c
c-----
        common/damp/alpha,ieqex
            real alpha
            integer ieqex
        common/jout/jsrc(21) , jbdrys, jbdryh
            integer jsrc, jbdrys, jbdryh
c-----
c       ray parameter values
c-----
        common/c/pmin,pmax,dp,pcntrl
        real pmin,pmax,dp,pcntrl
c-----
c       receiver model parameters
c-----
        integer LER, LIN, LOT
        parameter(LER=0, LIN=5, LOT=6)
c-----
c       earth model information
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
        integer mmax
        common/depref/refdep
        real refdep

        integer iunit, iiso, iflsph, idimen, icnvel, ierr
        character title*80

        common/lyrctl/lyrins
        logical lyrins
c-----
c       control of wavefield
c-----

        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar

        common/earth/radius
        real radius

c-----
c       internal variables
c-----
        integer i, n21, j
        real freq, df
        integer lmaxb, lmaxs
        real dphr
        real eta
        complex*16 gg(21)
        real sp, cp
        complex zstapha
        complex*16 zeye
        complex*16 ztmp
 
c-----
c       set the parameters normally used by hspec96p
c-----
        zeye = dcmplx(0.0d+00,1.00d+00)
        alpha = alp
        pmin = rayp
        pmax = rayp
        dp = rayp
        pcntrl = -1
        ieqex = 2
        jbdrys = 1
        jbdryh = 0
        lyrins = .true.

        dokjar = .true.


        n21 = npts/2 + 1
        df = 1.0/(npts*dt)
c-----
c       get the model at the source
c-----
        call getmod(1,modelsrc,mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
        if(ierr .lt. 0)stop
c-----
c       make sure that we use 1/Q
c-----
        do i=1,mmax
            if(qa(i).lt.0.0)qa(i) = 0.0
            if(qb(i).lt.0.0)qb(i) = 0.0
            if(qa(i) .gt. 1.0)qa(i) = 1.0/qa(i)
            if(qb(i) .gt. 1.0)qb(i) = 1.0/qb(i)
            if(frefp(i) .le. 0.0)frefp(i) = 1.0
            if(frefs(i) .le. 0.0)frefs(i) = 1.0
        enddo
        call modcpy(.true.)
c-----
c       check model for inconsistencies
c-----
        call chkmod()
c-----
c       insert a boundary at hstel and determine the index of the layer here
c       output for safety as a debug
c-----
        call insert(hstel+refdep)
c-----
c       get the index for the source depth
c-----
        call insert(hs+refdep)
        call dezero()
c-----
c       determine position of source and receive in new layers
c-----
        call srclyr(hstel+refdep,lmaxb,dphr)
        call srclyr(hs+refdep,lmaxs,dphr)
c-----
c       Change the maximum number of layers and
c       place teleseism model parameters into new base
c       layer - note in future preserve the Q etc
c-----
        mmax = lmaxb
        a(mmax)   = vsa
        b(mmax)   = vsb
        rho(mmax) = vsr
c-----
c       should also carry the Q etc
c-----
        WRITE(LOT,'(a)')
     1 '  LAYER   H(km)      PVel      SVel      Dens     QPi       QSi'
        do i=1,lmaxb
              if(i.ge.lmaxs .and. i.le.lmaxb)then
              WRITE(LOT,'(i5,4f10.2,2f10.4)') i,d(i),a(i),b(i),rho(i),
     1          qa(i),qb(i)
              else
              WRITE(LOT,'(i5,4f10.2,2f10.4)')-i,d(i),a(i),b(i),rho(i),
     1          qa(i),qb(i)
              endif
        enddo
c-----
c       transform the spherical model to a flat model
c-----
        if(iflsph.ne.0)then
            call adosph()
        endif
c-----
c       compute the vertical tau delay
c       For an explosion in water, there can be a P-S conversion 
c       at the fluid solid boundary that can generate S in the solid
c       So that is a special case
c-----
        srcdelay = 0.0
        srcdist  = 0.0
        srctime  = 0.0
        do i=lmaxs,lmaxb-1
           if(pteleseismic)then
                eta = sqrt(1.0 - rayp*rayp*a(i)*a(i))
                srcdelay = srcdelay + d(i)*eta/a(i)
                srcdist  = srcdist  + d(i)*rayp*a(i)/eta
                srctime  = srctime  + d(i)/(a(i)*eta)
           else
                if(b(i).lt.0.001*a(i))then
                    eta = sqrt(1.0 - rayp*rayp*a(i)*a(i))
                    srcdelay = srcdelay + d(i)*eta/a(i)
                    srcdist  = srcdist  + d(i)*rayp*a(i)/eta
                    srctime  = srctime  + d(i)/(a(i)*eta)
                else
                    eta = sqrt(1.0 - rayp*rayp*b(i)*b(i))
                    srcdelay = srcdelay + d(i)*eta/b(i)
                    srcdist  = srcdist  + d(i)*rayp*b(i)/eta
                    srctime  = srctime  + d(i)/(b(i)*eta)
                endif
           endif
        enddo
        WRITE(LOT,*)' srcdelay  :',srcdelay, 
     1     ' s  (used for source synthetic time shift)'
        WRITE(LOT,*)' srcdist   :',srcdist , 
     1     ' km (used for mantle path distance correction)'
        WRITE(LOT,*)' srctime   :',srctime , 
     1     ' s  (used for travel time computation)' 
        if(pteleseismic)then
             sp = rayp*a(lmaxb)
             cp = sqrt(1.0 -sp*sp)
        else
             sp = rayp*b(lmaxb)
             cp = sqrt(1.0 -sp*sp)
        endif

c-----
c       now that we have the material properties at the source and
c       receiver depths, apply the sphericity correction
cDEBUG
c-----
        do i=1,n21
             freq = (i-1)*df
             if(freq.lt.df) freq = 0.01*df
c-----
c            get the stationary phase factor
c            (omega/V) cos theta
c-----
             if(pteleseismic)then
                  zstapha = 
     1               cmplx(6.2831853*freq,-alpha)*cp/a(lmaxb)
             else
                  zstapha = 
     1               cmplx(6.2831853*freq,-alpha)*cp/b(lmaxb)
             endif
c-----
c            get the medium response for the Green's function integrands.
c            These are the Uz, Uz and Uphi terms from the propagator 
c            matrices. 
c-----
             call excitsrc(freq,rayp,gg,lmaxs,pteleseismic)
c-----
c            Now apply the method of stationary phase to
c            INT U(k,omega) Jn(kr) exp(-nu_V z) k dk
c-----
             gg( 1) =  gg( 1) * zstapha * zeye
             gg( 2) =  gg( 2) * zstapha 
             gg( 3) =  gg( 3) * zstapha * (-1.0d+00)
             gg( 4) =  gg( 4) * zstapha * zeye
             gg( 5) =  gg( 5) * zstapha * (-zeye)
             gg( 6) =  gg( 6) * zstapha * (-1.0d+00)
             gg( 7) =  gg( 7) * zstapha * zeye
             gg( 8) =  gg( 8) * zstapha 
             gg( 9) =  gg( 9) * zstapha * zeye
             gg(10) =  gg(10) * zstapha 
             gg(11) =  gg(11) * zstapha * (-1.0d+00)
             gg(12) =  gg(12) * zstapha * zeye
             gg(13) =  gg(13) * zstapha * zeye
             gg(14) =  gg(14) * zstapha * (-1.0d+00)
             gg(15) =  gg(15) * zstapha * zeye

c-----
c            The gg now are the far-field Green's functions as a function of
c            frequency. The corresponence is
c            1-ZDD   2-RDD   3-ZDS   4-RDS   5-ZSS   6-RSS
c            7-ZEX   8-REX   9-ZVF  10-RVF  11-ZHF  12-RHF  
c           13-TDS  14-TSS  15-THF
c-----
c       Finally, convert the uz, ur and ut to the teleseismic
c       Ap, Asv and Ash. For compatibility with the CPS hspec96 code
c       no longer distinguish between Z and R, but these values must be defined
c
c       If the downgoing ray makes an angle these with respect to the downward
c       positive z-axis, and if uz is positive down, ur positive away from source,
c       Then the direction of the positive unit P and SV vectors are
c             e_P  =\sin \theta e_r + \cos \theta e_z
c             e_SV =\cos \theta e_r - \sin \theta e_z
c
c       and the far-field displacements are
c             U    = ur e_r + uz e_z             
c-----
             do j=1,21
                srcp(j,i) = cmplx(0.0,0.0)
             enddo
             if(pteleseismic)then
C               GG(2) = - GG(2)
C               GG(8) = - GG(8)
C               GG(10) = - GG(10)
c-----
c                 DD
c-----
                 ZTMP       = (gg(1)*( cp) + gg(2)*sp)
                 srcp(1,i)  = cmplx(sngl(dreal(ZTMP)),sngl(dimag(ZTMP)))
                 srcp(2,i)  = srcp(1,i)
c-----
c                 DS
c-----
                 ZTMP       = (gg(3)*( cp) + gg(4)*sp)
                 srcp(3,i)  = cmplx(sngl(dreal(ZTMP)),sngl(dimag(ZTMP)))
                 srcp(4,i)  = srcp(3,i)
c-----
c                 SS
c-----
                 ZTMP       = (gg(5)*( cp) + gg(6)*sp)
                 srcp(6,i)  = cmplx(sngl(dreal(ZTMP)),sngl(dimag(ZTMP)))
                 srcp(7,i)  = srcp(6,i)
c-----
c                 EX
c-----
                 ZTMP       = (gg(7)*( cp) + gg(8)*sp)
                 srcp(9,i)  = cmplx(sngl(dreal(ZTMP)),sngl(dimag(ZTMP)))
                 srcp(10,i) = srcp(9,i)
c-----
c                 VF
c-----
                 ZTMP       = (gg(9)*( cp) + gg(10)*sp)
                 srcp(11,i) = cmplx(sngl(dreal(ZTMP)),sngl(dimag(ZTMP)))
                 srcp(12,i) = srcp(11,i)
c-----
c                 hF
c-----
                 ZTMP       = (gg(11)*( cp) + gg(12)*sp)
                 srcp(13,i) = cmplx(sngl(dreal(ZTMP)),sngl(dimag(ZTMP)))
                 srcp(14,i) = srcp(13,i)
             else
CRBH  amplitude OK but R flipped and Z not should be opposite
CRBH  because of zup vs down
C               GG(2) = - GG(2)
C               GG(8) = - GG(8)
C               GG(10) = - GG(10)

c-----
c                 DD
c-----
                 ZTMP       = (gg(1)*(-sp) + gg(2)*cp)
                 srcp(1,i)  = cmplx(sngl(dreal(ZTMP)),sngl(dimag(ZTMP)))
                 srcp(2,i)  = srcp(1,i)
c-----
c                 DS
c-----
                 ZTMP       = (gg(3)*(-sp) + gg(4)*cp)
                 srcp(3,i)  = cmplx(sngl(dreal(ZTMP)),sngl(dimag(ZTMP)))
                 srcp(4,i)  = srcp(3,i)
                 ZTMP       = gg(13)
                 srcp(5,i)  = cmplx(sngl(dreal(ZTMP)),sngl(dimag(ZTMP)))
c-----
c                 SS
c-----
                 ZTMP       = (gg(5)*(-sp) + gg(6)*cp)
                 srcp(6,i)  = cmplx(sngl(dreal(ZTMP)),sngl(dimag(ZTMP)))
                 srcp(7,i)  = srcp(6,i)
                 ZTMP       = gg(14)
                 srcp(8,i)  = cmplx(sngl(dreal(ZTMP)),sngl(dimag(ZTMP)))
c-----
c                 EX
c-----
                 ZTMP       = (gg(7)*(-sp) + gg(8)*cp)
                 srcp(9,i)  = cmplx(sngl(dreal(ZTMP)),sngl(dimag(ZTMP)))
                 srcp(10,i) = srcp(9,i)
c-----
c                 VF
c-----
                 ZTMP       = (gg(9)*(-sp) + gg(10)*cp)
                 srcp(11,i) = cmplx(sngl(dreal(ZTMP)),sngl(dimag(ZTMP)))
                 srcp(12,i) = srcp(11,i)
c-----
c                 HF
c-----
                 ZTMP       = (gg(11)*(-sp) + gg(12)*cp)
                 srcp(13,i) = cmplx(sngl(dreal(ZTMP)),sngl(dimag(ZTMP)))
                 srcp(14,i) = srcp(13,i)
                 ZTMP       = gg(15)
                 srcp(15,i) = cmplx(sngl(dreal(ZTMP)),sngl(dimag(ZTMP)))
             endif
               
        enddo
c-----
c       turn off for testing
c-----
        if(.not. dosrc)then
             srcdelay = 0.0
             srctime  = 0.0
             do j=1,21
                  do i=1,n21
                   srcp(j,i) = cmplx(1.0,0.0)
                  enddo
             enddo
        endif
        return
        end


        subroutine excitsrc(freq,rayp,gg,lmaxs,pteleseismic)
        implicit none
c-----
c       freq  - R frequency
c       rayp  - R ray parameter. If model in is km then rayp in s/km
c       gg    -   complex array of Green's function integrands
c       lmaxs - I index of source layer. Source is at the top of this layer
c       pteleseismic - L  .true. P down from source, .false. S down from source
c-----
c     subroutine arguments
c-----
        real freq, rayp
        complex*16 gg(21)
        integer lmaxs      
        logical pteleseismic
c-----
c     global arguments
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d,a,b,rho,qa,qb,etap,etas,frefp,frefs
        common/modlly/mmax
        integer mmax
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
        common/lyrctl/lyrins
        logical lyrins
c-----
c       internal parameters
c----
        real omega
        complex*16 wvno,om, wvno2, om2
        complex*16 spsv(6,4), ssh(6,2)
        integer nval(6)
c-----
c       vertical wavenumbers for the halfsapce
c-----
        complex*16  ran, rbn
        complex*16 xka, xkb
        complex*16 atna, atnb
        integer iwat

        integer i
  
        complex*16 xsh(2,2), rsh(2,2)
        complex*16 cdx(6,6), cdr(6,6), zpsv(4,4)

        complex*16 XZ(4)
        complex*16 KP, KSV
 
        integer ityp, jdisc
        real*8 exrpsv, exzpsv, exxpsv, exrsh, exsh
        real*8 exppsvcor, expshcor

c-----
c     spsv - P-SV source terms
c     ssh  - SH terms
c     nval - Bessel function order for far-field term
c    
c     indexing 
c       i = Green's function      
c     1   DD      
c     2   DS
c     3   SS
c     4   EX
c     5   VF
c     6   HF
c
c       j = Delta B component
c       P-SV                            SH
c     1 dUr  2 dUz 3 dTz 4 dTr   1 dUphi  2 dTphi
c-----


        omega=6.2831853*freq
        om=dcmplx(dble(omega),-dble(alpha))
        om2 = om * om
        wvno=dcmplx(dble(rayp),0.0d+00)*om
        wvno2 = wvno*wvno
c-----
c       evaluate the propagator matrices and then combine here according to the
c       source
c-----
        call aten(om,qa(mmax),qb(mmax),xka,xkb,alpha,a(mmax),
     1      b(mmax),atna,atnb,iwat,
     2      frefp(mmax),frefs(mmax))
        ran=CDSQRT(wvno2-xka*xka)
        rbn=CDSQRT(wvno2-xkb*xkb)
c-----
c       get the source terms for the  Grene's functions
c----- 
        call srcexc(lmaxs,om,wvno,spsv,ssh,nval)
c-----
c       for conversion of the U to displacment u and the u to amplitude,
c       get the ran, rbn vertical wavenumbers for the halfspace
c-----
c-----
c       get the matrices
c-----
        call gets(om,om2,wvno,wvno2,mmax,lmaxs,xsh,rsh,cdx,cdr,zpsv,
     1      exrpsv, exzpsv, exxpsv, exrsh, exsh)
        exppsvcor = dexp(exzpsv+exxpsv-exrpsv)
        expshcor  = dexp(exsh-exrsh)

c-----
c       apply the source terms to define the Green funcition integrands
c       SH
c           First get the Kd
c-----
            gg(13) = (xsh(2,1)*ssh(2,1) + xsh(2,2)*ssh(2,2))/rsh(2,2)
            gg(14) = (xsh(2,1)*ssh(3,1) + xsh(2,2)*ssh(3,2))/rsh(2,2)
            gg(15) = (xsh(2,1)*ssh(6,1) + xsh(2,2)*ssh(6,2))/rsh(2,2)
c-----
c          convert the Kd to Uphi
c-----
            gg(13) = gg(13) * wvno * expshcor
            gg(14) = gg(14) * wvno * expshcor
            gg(15) = gg(15) * wvno * expshcor
c-----
c       P-SV
c           First get the Kd
c-----
c-----
c       form the 4x1 vectors for KPD and KSVD - this is just the X Z
c       then dot with the source vector
c-----
        if(pteleseismic)then
c-----
c          Solve for KP down
c-----
           XZ(1) =                      CDX(6,1)*ZPSV(2,4) 
     1                     + CDX(6,2)*ZPSV(3,4) + CDX(6,3)*ZPSV(4,4)
           XZ(2) = -CDX(6,1)*ZPSV(1,4)                      
     1                     + CDX(6,4)*ZPSV(3,4) + CDX(6,5)*ZPSV(4,4)
           XZ(3) = -CDX(6,2)*ZPSV(1,4) - CDX(6,4)*ZPSV(2,4) 
     1                                          + CDX(6,6)*ZPSV(4,4)
           XZ(4) = -CDX(6,3)*ZPSV(1,4) - CDX(6,5)*ZPSV(2,4) 
     1                     - CDX(6,6)*ZPSV(3,4)
           do i=1,4
                 XZ(i) = exppsvcor*XZ(i)/cdr(6,6)
           enddo

           do ityp = 1,6
                KP = dcmplx(0.0d+00,0.0d+00)
                do jdisc = 1,4
                    KP = KP + SPSV(ityp,jdisc)*XZ(jdisc)
                enddo
                gg(2*ityp-1) =  - ran*KP
                gg(2*ityp  ) =   wvno*KP
           enddo

        else
c-----
c          Solve for KSV down
c-----

           XZ(1) =                      CDX(6,1)*ZPSV(2,3) 
     1                     + CDX(6,2)*ZPSV(3,3) + CDX(6,3)*ZPSV(4,3)
           XZ(2) = -CDX(6,1)*ZPSV(1,3)                      
     1                     + CDX(6,4)*ZPSV(3,3) + CDX(6,5)*ZPSV(4,3)
           XZ(3) = -CDX(6,2)*ZPSV(1,3) - CDX(6,4)*ZPSV(2,3) 
     1                                          + CDX(6,6)*ZPSV(4,3)
           XZ(4) = -CDX(6,3)*ZPSV(1,3) - CDX(6,5)*ZPSV(2,3) 
     1                     - CDX(6,6)*ZPSV(3,3)
           do i=1,4
                 XZ(i) = exppsvcor*XZ(i)/cdr(6,6)
           enddo
           do ityp = 1,6
                KSV = dcmplx(0.0d+00,0.0d+00)
c-----
c               S X Z gives -KSV, 
c-----
                do jdisc = 1,4
                    KSV = KSV - SPSV(ityp,jdisc)*XZ(jdisc)
                enddo
                gg(2*ityp-1) =   wvno*KSV
                gg(2*ityp  ) = - rbn *KSV

           enddo
        endif
        return
        end
        subroutine excitpw(freq,dop,rayp,Z,R,T,pvelrec,svelrec)
        
c-----
c       get surfae displacements at reiver crust due to an
c       incident unit amplitude P or S wae
c-----
c       Given      -1 
c             R = E   A    .... A    for P-SV, and    (13.3.4)
c                  N   N-1       1
c                   -1
c             r =  e  a    .... a    for SH           (13.3.5)
c                   N  N-1       1
c-----
c                                             |12
c       R = (   R  ) (1/ (omega/VP_N) ) ( 1/(R| )   A (inc)
c                22                           |12    P
c       Z = (-i R  ) (1/ (omega/VP_N) ) ( 1/(R| )   A (inc) 
c                21                           |12    P
c                                                            (13.3.6)
c                                             |12
c       R = (-i  R  ) (1/ (omega/VS_N) ) ( 1/(R| )  A (inc)
c                12                           |12    SV 
c                                             |12
c       Z = (-  R  ) (1/ (omega/VS_N) ) ( 1/(R| )   A (inc)
c                11                           |12    SV
c        
c        
c       T = (1 / k r  )  A( inc)
c                   11    SH
c       where R, T, and Z are the surface horizontal (radial), transverse and vertical 
c       displacements
c-----
        implicit none
        real freq, omega, rayp,pvelrec,svelrec
        logical dop
        complex Z, R, T

        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
            real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmax
            integer mmax
c-----
c       matrix components in layers and boundaries saved
c-----
        common/damp/alpha,ieqex
            real alpha
            integer ieqex

         complex*16 r1212
c-----
c       internal variables
c-----
        complex*16 om, om2, wvno, wvno2
        complex*16 zeye
        complex*16 cei(6,6)
        complex*16 dz, dr, dt
c-----
c       P-SV quantities
c-----
        complex*16 rpsv(4,4), rsh(2,2)
        real*8 expsv, exsh, exa 


        omega = 6.2831853*freq
        om =  dcmplx(dble(omega), dble(-alpha))
        om2 = om * om
        wvno  = dcmplx(dble(rayp),0.0d+00)*om
        wvno2 = wvno * wvno
        zeye = dcmplx(0.0d+00,1.0d+00)
C       WRITE(6,*)'EPW-freq,om,rayp,wvno:',
C    1     freq,om,rayp,wvno
C       WRITE(6,*)'epw-freq,dop,rayp,pvelrec,svelrec,alpha:',
C    1   freq,dop,rayp,pvelrec,svelrec,alpha

        call getr(rpsv,rsh,expsv,exsh,mmax,wvno,om,om2,wvno2,cei,
     1     exa)
C       r1212 = rpsv(1,1)*rpsv(2,2) - rpsv(1,2)*rpsv(2,1)
C       WRITE(6,*)'r1212:',r1212
        r1212 = cei(1,1)
C       WRITE(6,*)'r1212:',r1212
C       WRITE(6,*)'excitpw:om,expsv,exsh,exapsv',
C    1     om,expsv,exsh,exa
C       WRITE(6,*)' r1212    :',r1212
C       WRITE(6,*)' rpsv(1,1):',rpsv(1,1)
C       WRITE(6,*)' rpsv(2,1):',rpsv(2,1)
C       WRITE(6,*)' rpsv(1,2):',rpsv(1,2)
C       WRITE(6,*)' rpsv(2,2):',rpsv(2,2)

C       WRITE(6,*)'excitpw:om,r1212,rpsv,rsh,expsv,expsh',
C    1     om,r1212,rpsv,rsh,expsv,exsh
C       WRITE(6,*)'excitpw:om,rsh,exsh',
C    1     om,exsh
C       WRITE(6,*)' rsh(1,1):',rsh(1,1)
C       WRITE(6,*)' rsh(2,1):',rsh(2,1)
C       WRITE(6,*)' rsh(1,2):',rsh(1,2)
C       WRITE(6,*)' rsh(2,2):',rsh(2,2)
        
c-----
c       compute the receiver function
c-----
        if(dop)then
c-----
c           Note Uz is positive down
c-----
            DZ = -exp(expsv-exa)*zeye*rpsv(2,1)/r1212
            DZ = DZ * pvelrec /om
c-----
c           Ur
c-----
            DR = exp(expsv-exa)*rpsv(2,2)/r1212
            DR = DR * pvelrec /om
c-----
c           Ut
c-----
            DT = dcmplx(0.0d+00,0.0d+00)
        else
c-----
c           UZ
c-----
            DZ = -exp(expsv-exa)*rpsv(1,1)/r1212
            DZ = DZ * svelrec /om
c-----
c           Ur
c-----
            DR = -exp(expsv-exa)*zeye*rpsv(1,2)/r1212
            DR = DR * svelrec /om
c-----
c           Ut 
c-----
                 DT = exp(-exsh)/(wvno*rsh(1,1))
        endif
        Z = cmplx(sngl(dreal(DZ)),sngl(dimag(DZ)))
        R = cmplx(sngl(dreal(DR)),sngl(dimag(DR)))
        T = cmplx(sngl(dreal(DT)),sngl(dimag(DT)))
        return
        end


        subroutine localmod(mmax,refdep,d,a,b,rho,qa,qb,etap,etas,
     1   frefp,frefs,
     1   Lmmax,Lrefdep,Ld,La,Lb,Lrho,Lqa,Lqb,Letap,Letas,Lfrefp,Lfrefs,
     2   tolocal)
c-----
c       put from getmod storage to local storage is tolocal == .true.
c       else get
c-----
        integer mmax
        real refdep, d(*), a(*), b(*), rho(*), qa(*), qb(*)
        real etap(*), etas(*), frefp(*), frefs(*)
        integer Lmmax
        real Lrefdep, Ld(*), La(*), Lb(*), Lrho(*), Lqa(*), Lqb(*)
        real Letap(*), Letas(*), Lfrefp(*), Lfrefs(*)

        logical tolocal

c-----
c     copy from getmod input to local storage
c-----
        if(tolocal)then
              Lmmax = mmax
              Lrefdep = refdep
              do i=1,Lmmax
                  Ld(i) = d(i)
                  La(i) = a(i)
                  Lb(i) = b(i)
                  Lrho(i) = rho(i)
                  Lqa(i) = qa(i)
                  Lqb(i) = qb(i)
                  Letap(i) = etap(i)
                  Letas(i) = etas(i)
                  Lfrefp(i) = frefp(i)
                  Lfrefs(i) = frefs(i)
              enddo
        else
c-----
c     copy form local storage to putmod arrays
c     also do QC so that VS != 0, and Qinv is output
c     assuming that Q is always greater than 1
c-----
              mmax = Lmmax
              refdep = Lrefdep
              do i=1,mmax
                  d(i) = Ld(i)
                  a(i) = La(i)
                  if(Lb(i).eq.0.0)then
                       b(i) = 0.0001
                  else
                       b(i) = Lb(i)
                  endif
                  rho(i) = Lrho(i)
                  if(Lqa(i).gt.1.0)then
                       qa(i) = 1.0/Lqa(i)
                  else
                       qa(i) = Lqa(i)
                  endif
                  if(Lqb(i).gt.1.0)then
                       qb(i) = 1.0/Lqb(i)
                  else
                       qb(i) = Lqb(i)
                  endif
                  etap(i) = Letap(i)
                  etas(i) = Letas(i)
                  frefp(i) = Lfrefp(i)
                  frefs(i) = Lfrefs(i)
              enddo
        endif
        return
        end


        subroutine modmerge(modeltel,modelsrc,modelrec,
     1      numbdy,layerdepth, zsrc, zrec)
c-----
c       merge the models to create hudsonsrc.mod and hudsonrec.mod
c-----
c-----
c       velocity model parameters
c------
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep
c-----
c       temporary storage for teleseism model
c-----
        common/telmod/Td(NL),Ta(NL),Tb(NL),Trho(NL),
     1      Tqa(NL),Tqb(NL),Tetap(NL),Tetas(NL),
     2      Tfrefp(NL), Tfrefs(NL)
        real Td,Ta,Tb, Trho,Tqa,Tqb,Tetap,Tetas,Tfrefp,Tfrefs
        real Trefdep
        integer Tmmax
c-----
c       temporary storage for src model
c-----
        common/srcmod/Sd(NL),Sa(NL),Sb(NL),Srho(NL),
     1      Sqa(NL),Sqb(NL),Setap(NL),Setas(NL),
     2      Sfrefp(NL), Sfrefs(NL)
        real Sd,Sa,Sb, Srho,Sqa,Sqb,Setap,Setas,Sfrefp,Sfrefs
        real Srefdep
        integer Smmax
c-----
c       temporary storage for rec model
c-----
        common/recmod/Rd(NL),Ra(NL),Rb(NL),Rrho(NL),
     1      Rqa(NL),Rqb(NL),Retap(NL),Retas(NL),
     2      Rfrefp(NL), Rfrefs(NL)
        real Rd,Ra,Rb, Rrho,Rqa,Rqb,Retap,Retas,Rfrefp,Rfrefs
        real Rrefdep
        integer Rmmax
c-----
c       temporary storage for new model
c-----
        common/tmpmod/Nd(NL),Na(NL),Nb(NL),Nrho(NL),
     1      Nqa(NL),Nqb(NL),Netap(NL),Netas(NL),
     2      Nfrefp(NL), Nfrefs(NL)
        real Nd,Na,Nb, Nrho,Nqa,Nqb,Netap,Netas,Nfrefp,Nfrefs
        real Nrefdep
        integer Nmmax
c-----
c       command line parameters
c-----
        character modeltel*(*), modelsrc*(*), modelrec*(*)
        integer numbdy
        real layerdepth(3*NL), zsrc, zrec
c-----
c       parameters for reading the model
c-----
        integer mmax, iunit, iiso, iflsph, idimen, icnvel
        integer ierr
        character title*80
c-----
c       internal routine  variables
c-----
        integer i
c-----
c       get the teleseism model and store separately
c-----
        call getmod(1,modeltel,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
        call localmod(mmax,refdep,d,a,b,rho,qa,qb,etap,etas,frefp,frefs,
     1   Tmmax,Trefdep,Td,Ta,Tb,Trho,Tqa,Tqb,Tetap,Tetas,Tfrefp,Tfrefs,
     2   .true.)
c-----
c       get the src model and store separately
c-----
        call getmod(1,modelsrc,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
        if(d(mmax).eq.0.0)d(mmax) = zsrc
        call localmod(mmax,refdep,d,a,b,rho,qa,qb,etap,etas,frefp,frefs,
     1   Smmax,Srefdep,Sd,Sa,Sb,Srho,Sqa,Sqb,Setap,Setas,Sfrefp,Sfrefs,
     2   .true.)
c-----
c       get the rec model and store separately
c-----
        call getmod(1,modelrec,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
        if(d(mmax).eq.0.0)d(mmax) = zrec
        call localmod(mmax,refdep,d,a,b,rho,qa,qb,etap,etas,frefp,frefs,
     1   Rmmax,Rrefdep,Rd,Ra,Rb,Rrho,Rqa,Rqb,Retap,Retas,Rfrefp,Rfrefs,
     2   .true.)
c-----
c       now create the combined src/teleseismic model
c       using ugly linear searches
c
c       since the combined model will have a number of layers >=
c       to any original model, we focus on filling the combined
c       model
c-----
        Nrefdep = Trefdep
        Nmmax = numbdy -1
c-----
        
        do i=1,numbdy-1
          if(layerdepth(i+1).le. zsrc)then
              call getindex(layerdepth(i),layerdepth(i+1),Sd,Smmax,k)
              if(k.gt.0)then
                Nd(i) = layerdepth(i+1) - layerdepth(i)
                Na(i) = Sa(k)
                Nb(i) = Sb(k)
                Nrho(i) = Srho(k)
                Nqa(i) = Sqa(k)
                Nqb(i) = Sqb(k)
                Netap(i) = Setap(k)
                Netas(i) = Setas(k)
                Nfrefp(i) = Sfrefp(k)
                Nfrefs(i) = Sfrefs(k)
              else
c-----
c      added 13 MAY 2015
c-----
                Nd(i) = layerdepth(i+1) - layerdepth(i)
                Na(i) = Sa(Smmax)
                Nb(i) = Sb(Smmax)
                Nrho(i) = Srho(Smmax)
                Nqa(i) = Sqa(Smmax)
                Nqb(i) = Sqb(Smmax)
                Netap(i) = Setap(Smmax)
                Netas(i) = Setas(Smmax)
                Nfrefp(i) = Sfrefp(Smmax)
                Nfrefs(i) = Sfrefs(Smmax)
              endif
          else
              call getindex(layerdepth(i),layerdepth(i+1),Td,Tmmax,k)
              if(k.gt.0)then
                Nd(i) = layerdepth(i+1) - layerdepth(i)
                Na(i) = Ta(k)
                Nb(i) = Tb(k)
                Nrho(i) = Trho(k)
                Nqa(i) = Tqa(k)
                Nqb(i) = Tqb(k)
                Netap(i) = Tetap(k)
                Netas(i) = Tetas(k)
                Nfrefp(i) = Tfrefp(k)
                Nfrefs(i) = Tfrefs(k)
              endif
          endif
        enddo
c-----
c   lets look at the result
c-----
C        WRITE(0,*)'MERGED TELESEISM/SRC MODEL TOP'
C        do i=1,25
C        WRITE(0,'(i5,4f10.4)')i,Nd(i),Na(i),Nb(i),NRho(i)
C        enddo
C        WRITE(0,*)'..... ......... ......... ......... .........'
C        do i=Nmmax -24, Nmmax
C        WRITE(0,'(i5,4f10.4)')i,Nd(i),Na(i),Nb(i),NRho(i)
C        enddo
c-----
c       put this into the array for output
c-----
        call localmod(mmax,refdep,d,a,b,rho,qa,qb,etap,etas,frefp,frefs,
     1   Nmmax,Nrefdep,Nd,Na,Nb,Nrho,Nqa,Nqb,Netap,Netas,Nfrefp,Nfrefs,
     2   .false.)
        iflsph = 1
        call putmod(1,'hudsonsrc.mod',Nmmax,
     1       'Merged teleseism/src/model',
     2       iunit,iiso,iflsph,idimen,icnvel,.false.)
c-----
c       now create the combined rec/teleseismic model
c       using ugly linear searches
c
c       since the combined model will have a number of layers >=
c       to any original model, we focus on filling the combined
c       model
c-----
        Nrefdep = Trefdep
        Nmmax = numbdy -1
c-----
        
        do i=1,numbdy-1
          if(layerdepth(i+1).le. zrec)then
              call getindex(layerdepth(i),layerdepth(i+1),Rd,Rmmax,k)
              if(k.gt.0)then
                Nd(i) = layerdepth(i+1) - layerdepth(i)
                Na(i) = Ra(k)
                Nb(i) = Rb(k)
                Nrho(i) = Rrho(k)
                Nqa(i) = Rqa(k)
                Nqb(i) = Rqb(k)
                Netap(i) = Retap(k)
                Netas(i) = Retas(k)
                Nfrefp(i) = Rfrefp(k)
                Nfrefs(i) = Sfrefs(k)
              else
c-----
c      added 13 MAY 2015
c-----
                Nd(i) = layerdepth(i+1) - layerdepth(i)
                Na(i) = Ra(Rmmax)
                Nb(i) = Rb(Rmmax)
                Nrho(i) = Rrho(Rmmax)
                Nqa(i) = Rqa(Rmmax)
                Nqb(i) = Rqb(Rmmax)
                Netap(i) = Retap(Rmmax)
                Netas(i) = Retas(Rmmax)
                Nfrefp(i) = Rfrefp(Rmmax)
                Nfrefs(i) = Sfrefs(Rmmax)
              endif
          else
              call getindex(layerdepth(i),layerdepth(i+1),Td,Tmmax,k)
              if(k.gt.0)then
                Nd(i) = layerdepth(i+1) - layerdepth(i)
                Na(i) = Ta(k)
                Nb(i) = Tb(k)
                Nrho(i) = Trho(k)
                Nqa(i) = Tqa(k)
                Nqb(i) = Tqb(k)
                Netap(i) = Tetap(k)
                Netas(i) = Tetas(k)
                Nfrefp(i) = Tfrefp(k)
                Nfrefs(i) = Tfrefs(k)
              endif
          endif
        enddo
c-----
c   lets look at the result
c-----
C        WRITE(0,*)'MERGED TELESEISM/REC MODEL TOP'
C        do i=1,20
C        WRITE(0,'(i5,4f10.4)')i,Nd(i),Na(i),Nb(i),NRho(i)
C        enddo
C        WRITE(0,*)'..... ......... ......... ......... .........'
C        do i=Nmmax -19, Nmmax
C        WRITE(0,'(i5,4f10.4)')i,Nd(i),Na(i),Nb(i),NRho(i)
C        enddo
c-----
c       put this into the array for output
c-----
        call localmod(mmax,refdep,d,a,b,rho,qa,qb,etap,etas,frefp,frefs,
     1   Nmmax,Nrefdep,Nd,Na,Nb,Nrho,Nqa,Nqb,Netap,Netas,Nfrefp,Nfrefs,
     2   .false.)
        iflsph = 1
        call putmod(1,'hudsonrec.mod',Nmmax,
     1       'Merged teleseism/rec/model',
     2       iunit,iiso,iflsph,idimen,icnvel,.false.)
        return
        end

        subroutine npow2(npts)
c-----
c       Given npts, determine the N=2**m such that N >= npts
c       return the new ntps
c-----
        integer*4 nsamp, npts
        nsamp = npts
        npts = 1
 1000   continue
        if(npts.ge.nsamp)return
        npts = 2*npts
        go to 1000
        end

        subroutine srclay(depth,lmax,dph)
        parameter (LER=0, LIN=5, LOT=6)
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        common/modlly/mmax
        common/lyrctl/lyrins
        logical lyrins
        if(.not.lyrins)then
            call modcpy(.false.)
            call insert(depth)
        endif
        call srclyr(depth,lmax,dph)
        return
        end


      subroutine srcexc(m,om,wvno,spsv,ssh,nval)
      implicit none

c-----
c     define the source terms
c-----
c     m    - source layer
c     om   - complex angulere frequency
c     wvno - wavenumber
c     spsv - P-SV source terms
c     ssh  - SH terms
c     nval - Bessel function order for far-field term
c     ra   - sqrt(k^2 - ka^2) for source
c     rb   - sqrt(k^2 - kb^2) for source
c    
c     indexing 
c       i         
c     1   DD      
c     2   DS
c     3   SS
c     4   EX
c     5   VF
c     6   HF
c
c       j
c       P-SV                            SH
c     1 dUr  2 dUz 3 dTz 4 dTr   1 dUphi  2 dTphi
c-----
      integer m
      complex*16 om, wvno
      complex*16 spsv(6,4), ssh(6,2)
      integer nval(6)
c-----
c     global parameters
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b,rho,qa,qb,etap,etas,frefp,frefs

        common/damp/alpha,ieqex
        real alpha
        integer ieqex

        common/lwater/lfluid
        logical lfluid

c-----
c     internal parameters
c-----
      real twopi
      complex*16 zone
      integer i, j
      integer iwat
      complex*16 xka,xkb,atna,atnb
      complex*16 mu, l3mu2, l2mu
c-----
c     initialize
c-----
      twopi = 6.2831853
      zone = dcmplx(1.0d+00,00d+00)
      do i=1,6
         do j=1,4
            spsv(i,j) = dcmplx(0.0d+00,0.0d+00)
         enddo
         do j=1,2
            ssh (i,j) = dcmplx(0.0d+00,0.0d+00)
         enddo
      enddo
       call aten(om,qa(m),qb(m),xka,xkb,
     1              alpha,a(m),b(m),atna,atnb,iwat,
     2              frefp(m),frefs(m))
      mu   = rho(m)*(atna*b(m))**2
      l2mu = rho(m)*(atna*a(m))**2
      l3mu2 = 3*l2mu -4*mu
     
c --- DD
      spsv(1,2) = 2*zone/(twopi*l2mu)
      spsv(1,4) = -wvno*l3mu2/(twopi*l2mu)
c --- DS
      spsv(2,1) = zone/(twopi*mu)
c --- SS
      spsv(3,4) = -wvno/twopi
c --- EX
      spsv(4,2)= zone/(twopi*l2mu)
      spsv(4,4)= 2*wvno*mu/(twopi*l2mu)
c --- VF
      spsv(5,3)= -zone/twopi
c --- HF
      spsv(6,4) = -zone/twopi

c --- DS
      ssh (2,1) = -zone/(twopi*mu)
c --- SS
      ssh (3,2) =  wvno/twopi
c --- HF
      ssh (6,2) =  zone/twopi
      nval(1) = 0
      nval(2) = 1
      nval(3) = 2
      nval(4) = 0
      nval(5) = 0
      nval(6) = 1
      return
      end
        subroutine defineprop(modeltel,modelsrc,modelrec,hs,
     1      hstel, hrtel, zsrc, zrec,
     2      vpsrc,vssrc,densrc, vprec,vsrec,denrec)
c-----
c       read the models
c       determine unique depths to interfaces
c       create the merged models
c
c       Input:
c       modeltel    C*   teleseismic model
c       modelsrc    C*   layered model in source region
c       modelrec    C*   layered model in recweiver region
c       hs          R    source depth
c       zsrc        R    minimum depth of source model
c       rsrc        R    minimum depth of receiver model
c
c       Outpur:
c       hstel       R    Depth of base in source region
c       hsrec       R    Depth of base in receiver region
c       vpsrc       R    P velocity of source  model at depth of hstel
c       vssrc       R    S velocity of source  model at depth of hstel
c       densrc      R    Density    of source  model at depth of hstel
c       vprec       R    P velocity of receiver  model at depth of hstel
c       vsrec       R    S velocity of receiver  model at depth of hstel
c       denrec      R    Density    of receiver  model at depth of hstel
c-----
        implicit none
c-----
c       routine parameters
c-----
        character modeltel*(*)
        character modelsrc*(*)
        character modelrec*(*)
        real hs
        real hstel, hrtel, zsrc, zrec
        real vpsrc, vssrc, densrc
        real vprec, vsrec, denrec
c-----
c       velocity model parameters
c------
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep

        integer mmax, iunit, iiso, iflsph, idimen, icnvel
        integer ierr
        character title*80

        real layerdepth(3*NL)

c-----
c       internal variables
c-----
        real depth
        integer i,j

c-----
c       determine the base of the receiver model if not modeltel
c       also it must be a minimum of zrec km deep which should be upper mantle
c-----
        if(modelrec.eq.modeltel)then
             hrtel = zrec
        else
             call getmod(1,modelrec,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
c-----
c            compute the depth to the base of the model
c-----
             depth = 0.0
             do i=1,mmax
                 depth = depth + d(i)
             enddo
             if(depth.lt.zrec)then
                   hrtel = zrec
             else
                   hrtel = depth 
             endif
        endif
c-----
c       determine the base of the source model if not modeltel
c       also it must be a minimum of zsrc km deep which should be upper mantle
c-----
        if(modelsrc.eq.modeltel)then
             hstel = zsrc
        else
             call getmod(1,modelsrc,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
c-----
c            compute the depth to the base of the model
c-----
             depth = 0.0
             do i=1,mmax
                 depth = depth + d(i)
             enddo
             if(depth.lt.zsrc)then
                   hstel = zsrc
             else
                   hstel = depth
             endif
        endif
        if(hstel.le.hs)then
             hstel = hs + zsrc
        endif
c----- 
c       hstel and hrtel are the depthss to the base of the
c       layered structyres at the source and receiver locations
c       
c-----
c            get the velocities and densities at depths of hstel and hrtel from
c            the teleseismic model
c-----
        call getbase(modeltel,hstel,hrtel,vpsrc,vssrc,densrc,
     1    vprec,vsrec,denrec) 
c-----
c       build the hudsonsrc.mod hudsonrec.mod model files
c       These must be spherical
c       These must have the same number of layers
c       read all models, compute depths, sort, make uniq
c
c   beware of refdep
c-----
        do i=1,3*NL
             layerdepth(i) = 0.0
        enddo
        j=1
c-----
c            get boundaries for teleseism model
c-----
             call getmod(1,modeltel,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
             depth = 0
             do i=1,mmax 
                   layerdepth(j) = depth
                   j = j + 1
                   depth = depth + d(i)
             enddo
             layerdepth(j) = depth
c-----
c            get boundaries for receiver shallow model
c-----
             call getmod(1,modelrec,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
             depth = 0
             j = j + 1
             do i=1,mmax 
                   layerdepth(j) = depth
                   j = j + 1
                   depth = depth + d(i)
             enddo
             layerdepth(j) = depth
             j = j + 1
             layerdepth(j) = zsrc
             j = j + 1
             layerdepth(j) = zrec
c-----
c            get boundaries for source shallow model
c-----
             call getmod(1,modelsrc,mmax,title,iunit,iiso,iflsph,idimen,
     1            icnvel, ierr,.false.)
             depth = 0
             j = j + 1
             do i=1,mmax 
                   layerdepth(j) = depth
                   j = j + 1
                   depth = depth + d(i)
             enddo
             layerdepth(j) = depth
c-----
c       now sort and get uniq
c-----
            call bsort(j,layerdepth,+1)

            call uniq(j,layerdepth)
c-----
c           output the layer thicknesses of the model
c-----
c-----
c       now do tel+rec
c       tel first and then top fill rec
c-----
        call modmerge(modeltel,modelsrc,modelrec, j, 
     1          layerdepth, zsrc, zrec)

        return
        end

      subroutine getbase(modeltel,hstel,hrtel,vpsrc,vssrc,densrc,
     1    vprec,vsrec,denrec)
      implicit none
c-----
c     subroutine arguments
c-----

      character modeltel*(*)
      real hstel,hrtel
      real vpsrc,vssrc,densrc
      real vprec,vsrec,denrec

c-----
c       velocity model parameters
c------
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b, rho,qa,qb,etap,etas,frefp,frefs
        common/depref/refdep
        real refdep

        integer mmax, iunit, iiso, iflsph, idimen, icnvel
        integer ierr
        character title*80
c-----
c       local variables
c-----
        real dlow, dhigh
        integer i
        

        call getmod(1,modeltel,mmax,title,iunit,iiso,iflsph,idimen,
     1        icnvel, ierr,.false.)
        if(iflsph.ne.0)then
            call adosph()
        endif

        dlow = 0.0
        do i=1,mmax
           dhigh = dlow + d(i)
           if(hstel.ge.dlow .and. hstel.le.dhigh)then
              vpsrc  = a(i)
              vssrc  = b(i)
              densrc = rho(i)
           endif
           if(hrtel.ge.dlow .and. hrtel.le.dhigh)then
              vprec  = a(i)
              vsrec  = b(i)
              denrec = rho(i)
           endif
           dlow = dhigh
        enddo
        

        return
        end
      subroutine getr(rpsv,rsh,expsv,exsh,mmax,wvno,om,om2,wvno2,cr,
     1     exapsv)
c-----
c     get the R matrix for P-SV and SH 
c            -1
C       R = E   a    ... a
c            N   N-1      1
c     for use with the receiver response. In terms of the incident potentials, we 
c     must solve the free surface problem
c
c     | Ku |      | U |
c     |    |  = R |   |
c     | Kd |      | 0 |
c
c     where the K and U matrices are 1x1 for and 2x1 for P-SV.  The R is 2x2 for SH
c     and 4x4 for P-SV. For the SH problem
c     we need r11 and for the P-SV problem R11, R12, R21, R22 and the 
c     sub-determinant  R11 R22 - R12 R21.    
c     rpsv = 4x4 P-SV R matrix
c     rsh  = 2x2 SH   R matrix
c     expsv - exponent term of the exp(expsv) factored out of the 4x4 R matrix for P-SV
c     expsh - exponent term of the exp(expsh) factored out of the 2x2 R matrix for SH
c     mmax  - layer index of the halfspace boundary for the receiver structure. This 
c          means that conversions from deeper layers will not be used to create the
c          receiver response
c     wvno - complex wavenumber
c     om   - complex angular frequency
c     om2  - om squared
c     wvno2  - wvno squared
c     cr  - first row  compound R matrix. cr(1,1) is the numerically
c          stable sub-determinant R11 R22 - R12 R21
c     exapsv - exponent term of the exp(exapsv) factoed out of the compount R matrix
c-----
      implicit none
c-----
c     subroutine arguments
c-----
      complex*16 rpsv(4,4), rsh(2,2)
      complex*16 wvno,om,wvno2,om2
      complex*16 cr(6,6)
      real*8 expsv, exsh, exapsv
      integer mmax
c-----
c     global parameters
c-----
        common/lwater/lfluid
        logical lfluid
c-----
c     internal variables
c-----
      integer m
      complex*16 apsv(4,4), ash(2,2)
      complex*16 tpsv(4,4), tsh(2,2)
      real*8 ex,exb,exa
      complex*16 tei(6,6)
      complex*16 capsv(6,6)
c-----
c     initialize matrices
c-----
      call identity2(rsh)
      call identity4(rpsv)
      call identity6(cr)

     
       

c-----
c      halfspace
c-----
       expsv = 0.0d+00
       exsh  = 0.0d+00
       exapsv = 0.0d+00
       call emati(mmax,wvno,om,om2,wvno2,rpsv,rsh)
       call ceinv(mmax,om,om2,wvno,wvno2,cr)
       

C      WRITE(6,*)'rpsv:',rpsv
C      WRITE(6,*)'rsh:',rsh
C      WRITE(6,*)'cr:',cr
c      now multiply my layer matrices
c      Now form R = E_N^{-1} a_{N-1} ... a_1
c               r = e_N^{-1} a_{N-1} ... a_1
c             C(R)= C(E_N^{-1}) C(a_{N_1}) ... C(a_1)
c       where C( ) is the compound matrix
c-----
       do m=mmax-1,1,-1
           call hask(m,om,om2,wvno,wvno2,apsv,ash,ex,exb)
           expsv = expsv + ex
           call multabc4(rpsv,apsv,tpsv)
           call copyab4(tpsv,rpsv)

           exsh  = exsh  + exb
           call multabc2(rsh,ash,tsh)
           call copyab2(tsh,rsh)

           call chask(m,capsv,om,wvno,wvno2,om2,exa)
           exapsv = exapsv + exa
           call multabc6(cr,capsv,tei)
           call copyab6(tei,cr)

       enddo
C      WRITE(6,*)'rpsv:',rpsv
C      WRITE(6,*)'rsh:',rsh
C      WRITE(6,*)'cr:',cr
C      WRITE(6,*)'cr(1,1):',cr(1,1)
       return
       end
      subroutine gets(om,om2,wvno,wvno2,mmax,lmaxs,xsh,rsh,cdx,cdr,zpsv,
     1     exrpsv, exzpsv, exxpsv, exrsh, exxsh)
c-----
c     get the matrices required for the source generation of
c     waves into the halfspace from a source in the layer stack. 
c     The equation to be solved is
c
c     | 0  |          | U |
c     |    |  = XS +R |   |
c     | Kd |          | 0 |
c
c     where the K and U matrices are 1x1 for SH and 2x1 for P-SV. The X and R 
c     are 2x2 for SH and 4x4 for P-SV. The S is the displacement-stress 
c     discontinuity at the source boundary and is 2x1 fore SH and 4x1 for P-SV.
c     Here
c          -1                         -1
c     R = E  a   ...  a      and X = E   a   ... a
c          N  N-1      1              N   N-1     lmaxs
c
c     To get the potentials in the halfspace, rerrange the equations to 
c
c     | U |    -1  -1      -1      | 0 |    -1  -1      -1
c     |   | = a   a   ... a    E   |   | - a   a   ... a        S   , or
c     | 0 |    1   2       N-1  N  | Kd|    1   2       lmaxs-1
c
c
c     | U |      | 0 |   
c     |   | = R' |   | - X' S  
c     | 0 |      | Kd|  
c
c     Focusing on the last row to solve for Kd gives
c
c          X' S
c           2j j  
c     Kd = ------  for SH. The P-SV problems requires the solution of
c          R'
c           22
c
c     |  0 |   | R'   R'   |  |     |   |  X'  S  |
c     |    |   |  33   34  |  | KdP |   |   3j  j |
c     |    | = |           |  |     | - |         |
c     |  0 |   | R'   R'   |  | KdSV|   |  X'  S  |
c     |    |   |  43   44  |  |     |   |   4j  j |
c
c     where the terms in j are summations.  To solve for the Kd the
c     sub-determinant  
c                     R'  R'  - R'  R'    is required. This
c                      33  44    34  43
c     is provided by the compound matrix of R', e.g., C(R')(6,6)
c     To start the multiplication, one just needs the 6th column
c     of C(E )
c           N
c-----
c     lmaxs - index of source layer where values will be saved
c-----
      implicit none
c-----
c     subroutine arguments
c-----
      complex*16 zpsv(4,4), rsh(2,2),xsh(2,2)
      complex*16 cdx(6,6), cdr(6,6)
      complex*16 wvno,om,wvno2,om2
      integer lmaxs
      integer mmax
      real*8 exrpsv, exzpsv, exxpsv, exrsh, exxsh
c-----
c     global parameters
c-----
        common/lwater/lfluid
        logical lfluid
c-----
c     internal variables
c-----
      integer m
      complex*16 apsvi(4,4), ashi(2,2)
      complex*16 capsvi(6,6)
      complex*16 tpsv(4,4), tsh(2,2)
      complex*16 tcpsv(6,6)
      real*8 ex,exb,exa
      complex*16 esh(2,2), epsv(4,4), cepsv(6,6)
c-----
c      initialize matrices
c-----
       call identity2(xsh)
       call identity2(rsh)
       call identity4(zpsv)
       call identity6(cdr)
       call identity6(cdx)
c-----
c      intialize normalixation factors
c-----
c NEED exponent for e and for r
c NEED exponent for X Z and R
       exrpsv = 0.0d+00
       exzpsv = 0.0d+00
       exxpsv = 0.0d+00
       exrsh  = 0.0d+00
       exxsh  = 0.0d+00


       do m=1,mmax-1
          call haski(m,om,om2,wvno,wvno2,apsvi,ashi,ex,exb)
          call chaski(m,capsvi,om,wvno,wvno2,om2,exa)
          if(m.lt.lmaxs)then
              call multabc6(cdx,capsvi,tcpsv)
              call copyab6(tcpsv,cdx)
              exxpsv = exxpsv + exa
              call anorm6(cdx,exxpsv)
              call multabc2(xsh,ashi,tsh)
              call copyab2(tsh,xsh)
              exxsh = exxsh + exb
              call anorm2(xsh,exxsh)
          else
              call multabc4(zpsv,apsvi,tpsv)
              call copyab4(tpsv,zpsv)
              call anorm4(zpsv,exzpsv)
              exzpsv = exzpsv + ex
          endif
          call multabc2(rsh,ashi,tsh)
          call copyab2(tsh,rsh)
          exrsh  = exrsh  + exb
          call anorm2(rsh,exrsh)
          call multabc6(cdr,capsvi,tcpsv)
          call copyab6(tcpsv,cdr)
          call anorm6(cdr,exrpsv)
          exrpsv = exrpsv + exa
       enddo
c-----
c       put in halfspace terms
c----- 
        call emat(mmax,wvno,om,om2,wvno2,epsv,esh)
c-----
c       SH e and r matrices
c-----
        call multabc2(rsh,esh,tsh)
        call copyab2(tsh,rsh)
c-----
c       P SV Z matrix
c-----
        call multabc4(zpsv,epsv,tpsv)
        call copyab4(tpsv,zpsv)
c-----
c       P SV compound R matrix
c-----

        call ce(mmax,om,om2,wvno,wvno2,cepsv)
        call multabc6(cdr,cepsv,tcpsv)
        call copyab6(tcpsv,cdr)

       return
       end
        subroutine aten(omega,qa,qb,xka,xkb,alpha,a,b,atna,atnb,iwat,
     2      frefp,frefs)
c-----
c       make velocities complex, using Futterman causality operator
c-----
        real*4 qa,qb,alpha,a,b
        complex*16 omega,at,atna,atnb,xka,xkb
        real*8 pi, om1p, om1s, oml, fac, pi2
        common/cntrl/ishank,hnkarg,dstcor,dokjar
        integer*4 dstcor
        real*4 hnkarg
        logical ishank,dokjar
        real*8 CDABS
        complex*16 CDLOG
c-----
c       reference frequency is fref hz
c-----
        om1p=6.2831853*frefp
        om1s=6.2831853*frefs
        pi2 = 1.5707963
        pi=3.1415927d+00
        if(dokjar)then
c-----
c       Kjartansson Constant Q, causal Q operator
c       Kjartansson, E. (1979). 
c          Constant Q-wave propagation and attenuation,
c       J. Geophys. Res. 84, 4737-4748.
c-----
            gama = atan(qa)/sngl(pi)
            gamb = atan(qb)/sngl(pi)
            if(gama.le.0.0)then
                atna = cmplx(1.0,0.0)
            else
                fac = pi2*gama
                rfac = sngl(dsin(fac)/dcos(fac))
                atna = dcmplx(1.0d+00,0.0d+00)/
     1              (( (omega/om1p)**dble(-gama) ) *
     2               dcmplx(1.0d+00,-dble(rfac)))
            endif
            if(b.gt.1.0e-04*a)then
                if(gamb.le.0.0)then
                    atnb = cmplx(1.0,0.0)
                else
                    fac = pi2*gamb
                    rfac = sngl(dsin(fac)/dcos(fac))
                    atnb = dcmplx(1.0d+00,0.0d+00)/
     1              (( (omega/om1s)**dble(-gamb) ) *
     2               dcmplx(1.0d+00,-dble(rfac)))
                endif
            endif
        else
c-----
c       Futterman Causal Q
c-----
c           low frequency cutoff is 0.01 hz
c-----
            oml=0.062831853d+00
            atna=dcmplx(1.0d+00,0.0d+00)
            atnb=dcmplx(1.0d+00,0.0d+00)
            if(qa.gt.0.0)then
              at=dcmplx(0.0d+00,0.0d+00)
              if(CDABS(omega).gt.oml) at=CDLOG(omega/om1p)/pi
              if(CDABS(omega).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1p*fac))/pi
              endif
              atna=(1.+dble(qa)*at+dcmplx(0.0d+00,dble(qa/2.)))
            endif
            if(qb.gt.0.0)then
              at=dcmplx(0.0d+00,0.0d+00)
              if(CDABS(omega).gt.oml) at=CDLOG(omega/om1s)/pi
              if(CDABS(omega).le.oml) then
                fac=dsqrt(oml*oml + dble(alpha*alpha))/oml
              at=CDLOG(dcmplx(dble(oml),-dble(alpha))/(om1s*fac))/pi
              endif
               atnb=(1.+dble(qb)*at+dcmplx(0.0d+00,dble(qb/2.)))
            endif
        endif
        xka=omega/(dble(a)*atna)
        if(b.le.1.0e-04*a)then
            iwat = 1
            xkb = dcmplx(0.0d+00,0.0d+00)
        else
            iwat = 0
            xkb=omega/(dble(b)*atnb)
        endif
        return
        end

c-----
      subroutine ce(m,om,om2,wvno,wvno2,h)
c-----
c                       
c     last column of reduced compound E
c                                       m
c-----
      implicit none
c-----
c     command line arguments
c-----
      integer m
      complex*16 om,om2,wvno,wvno2,h(6,6)
c-----
c     global parameters
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b,rho,qa,qb,etap,etas,frefp,frefs

        common/modlly/mmax
        integer mmax

        common/damp/alpha,ieqex
        real alpha
        integer ieqex
        integer i,j

c-----
c     internal parameters
c-----
      complex*16 ra,rb,gam,gamm1,xka,xkb
      complex*16 atna, atnb
      integer iwat

      do i=1,6
         do j=1,6
            h(i,j) = dcmplx(0.0d+00,0.0d+00)
         enddo
      enddo

        call aten(om,qa(m),qb(m),xka,xkb,alpha,a(m),
     1      b(m),atna,atnb,iwat,
     2      frefp(m),frefs(m))
        ra=CDSQRT(wvno2-xka*xka)
        rb=CDSQRT(wvno2-xkb*xkb)
        gam = dble(b(m))*wvno/om
        gam = gam * atnb
        gam = 2.0d+00 * (gam * gam)
        gamm1 = gam - dcmplx(1.0d+00,0.0d+00)
c-----
c     5x1 reduced compound H=E  matrix
c     related to the original 6x1 submatrix by c     [11,21,31,41,51,61]=[11',21',31'/2,-31'/2,41',51']
c-----

        h(1,6)=(wvno2-ra*rb)
        h(2,6)=-dble(rho(m))*om2*rb
        h(3,6)=rho(m)*om2*
     1      (-gam*(ra*rb/wvno) +wvno*gamm1)
        h(4,6) =- h(3,6)
        h(5,6)=-dble(rho(m))*om2*ra
        h(6,6)=dble(rho(m)*rho(m))*om2*om2*
     1      (-gam*gam*ra*rb/wvno2+gamm1*gamm1)
c-----
c      only those above are required. The rest are givenhere for completeness
c-----
        h(1,1) = wvno2-ra*rb
        h(1,2) = -2.*wvno*ra
        h(1,3) = wvno2 + ra*rb
        h(1,4) = -h(1,3)
        h(1,5) = 2*wvno*rb
        h(1,6)=(wvno2-ra*rb)

        h(2,1) = rho(m)*om2*rb
        h(2,2) = dcmplx(0.0d+00,0.0d+00)
        h(2,3) = - rho(m)*om2*rb
        h(2,4) = - rho(m)*om2*rb
        h(2,5) =  dcmplx(0.0d+00,0.0d+00)
        h(2,6)=-dble(rho(m))*om2*rb

        h(3,1)=rho(m)*om2*
     1      (-gam*(ra*rb/wvno) +wvno*gamm1)
        h(3,2) = -2*rho(m)*om2*ra*gam
        h(3,3) = rho(m)*om2*
     1      ( gam*(ra*rb/wvno) +wvno*gamm1)
        h(3,4) = -h(3,3)
        h(3,5) = 2*rho(m)*om2*rb*gamm1
        h(3,6)=rho(m)*om2*
     1      (-gam*(ra*rb/wvno) +wvno*gamm1)
    
        h(4,1) = rho(m)*om2*
     1      ( gam*(ra*rb/wvno) -wvno*gamm1)
        h(4,2) = 2*rho(m)*om2*ra*gamm1
        h(4,3) = rho(m)*om2*
     1      (-gam*(ra*rb/wvno) -wvno*gamm1)
        h(4,4) = -h(4,3)
        h(4,5) = -2*rho(m)*om2*gam*rb
        h(4,6) = - h(3,6)

        h(5,1) = - rho(m)*om2*ra
        h(5,2) = dcmplx(0.0d+00,0.0d+00)
        h(5,3) = - rho(m)*om2*ra
        h(5,4) = - rho(m)*om2*ra
        h(5,5) = dcmplx(0.0d+00,0.0d+00)
        h(5,6) = dble(rho(m))*om2*ra

        h(6,1) = dble(rho(m)*rho(m))*om2*om2*
     1      (-gam*gam*ra*rb/wvno2+gamm1*gamm1)
        h(6,2) = -2.*dble(rho(m)*rho(m))*om2*om2*
     1         gam*gamm1*ra/wvno
        h(6,3) = dble(rho(m)*rho(m))*om2*om2*
     1     (gamm1*gamm1 +gam*gam*ra*rb/wvno2)

        h(6,4) = - h(6,3)
        h(6,5) = 2*dble(rho(m)*rho(m))*om2*om2*
     1      gam*gamm1*rb/wvno
        h(6,6)=dble(rho(m)*rho(m))*om2*om2*
     1      (-gam*gam*ra*rb/wvno2+gamm1*gamm1)
      return
      end
      subroutine ceinv(m,om,om2,wvno,wvno2,g)
c-----
c                                    -1
c     first row of reduced compound E
c                                    m
c-----
      implicit none
c-----
c     command line arguments
c-----
      integer m
      complex*16 om,om2,wvno,wvno2,g(6,6)
c-----
c     global parameters
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b,rho,qa,qb,etap,etas,frefp,frefs

        common/modlly/mmax
        integer mmax

        common/damp/alpha,ieqex
        real alpha
        integer ieqex

c-----
c     internal parameters
c-----
      complex*16 ra,rb,gam,gamm1,xka,xkb
      complex*16 atna, atnb
      integer iwat
      integer i,j

        call aten(om,qa(m),qb(m),xka,xkb,alpha,a(m),
     1      b(m),atna,atnb,iwat,
     2      frefp(m),frefs(m))
        ra=CDSQRT(wvno2-xka*xka)
        rb=CDSQRT(wvno2-xkb*xkb)
        gam = dble(b(m))*wvno/om
        gam = gam * atnb
        gam = 2.0d+00 * (gam * gam)
        gamm1 = gam - dcmplx(1.0d+00,0.0d+00)
c-----
c     5x1 reduced compound G=E inverse matrix
c     related to the original 6x1 submatrix by
c     [11,12,13,14,15,16]=[11',12',13',-13',14',15']
c-----
        do i=1,6
           do j=1,6
             g(i,j)=dcmplx(0.0d+00,0.0d+00)
           enddo
        enddo

        g(1,1)=dble(rho(m)*rho(m))*om2*om2*
     1      (-gam*gam*ra*rb+wvno2*gamm1*gamm1)
        g(1,2)=-dble(rho(m))*(wvno2*ra)*om2
        g(1,3)=  dble(rho(m))*( gam*ra*rb-wvno2*gamm1)
     1      *om2*wvno
        g(1,4)= -g(1,3)
        g(1,5)=dble(rho(m))*(wvno2*rb)*om2
        g(1,6)=wvno2*(wvno2-ra*rb)

        g(1,1) = 0.25*g(1,1)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        g(1,2) = 0.25*g(1,2)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        g(1,3) = 0.25*g(1,3)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        g(1,4) = 0.25*g(1,4)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        g(1,5) = 0.25*g(1,5)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
        g(1,6) = 0.25*g(1,6)/(-rho(m)*rho(m)*om2*om2*wvno2*ra*rb)
      return
      end

c-----
c       6x6 compound matrices for P-SV
c-----
        subroutine chask(m,ca,om,wvno,wvno2,om2,exa)
        implicit none
c-----
c       subroutine arguments
c-----
        integer m
        complex*16 ca(6,6)
        complex*16 om,wvno, wvno2, om2
        real*8 exa
c-----
c     global parameters
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b,rho,qa,qb,etap,etas,frefp,frefs

        common/ ovrflw / a0,cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz
        real *8 a0
        complex*16 cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz

        common/damp/alpha,ieqex
        real alpha
        integer ieqex

        common/lwater/lfluid
        logical lfluid
c-----
c       interval variables
c-----

        complex*16 xka,xkb,ra,rb
        complex*16 atna,atnb
        complex*16 h
        integer iwat
        complex*16 p,q
        real*8 exb, ex
        complex*16 w,x,y,z,cosp,cosq,yl,zl,cosql
        integer i,j



        complex*16 gam,gam2,gamm1,gamm2,a0c,xz2,wy2,temp
        complex*16 cqww2, cqxw2, g1wy2, gxz2, g2wy2, g2xz2
        complex*16 gg1, a0cgg1
        complex*16 zrho, zrho2

       complex*16 zone
c----- 
c      Structure of 6x6 compound matrix  
c        A11     A12     A13    -A13     A15     A16
c        A21     A22     A23    -A23     A25     A15
c        A31     A32     A33    1-A33   -A23    -A13
c       -A31    -A32    1-A33    A33     A23     A13
c        A51     A52    -A32     A32     A22     A12
c        A61     A51    -A31     A31     A21     A11
c-----
c       this will be multipled on the left by the G matrix
c
c       [ G11   G12 G13 -G13    G15 G16 ]
c
c-----
c       or on the right by
c
c       [ H11   H21 H31 -H31    H51 H61  ] ^T
c-----


       zone = dcmplx(1.0d+00, 0.0d+00)
       call aten(om,qa(m),qb(m),xka,xkb,
     1              alpha,a(m),b(m),atna,atnb,iwat,
     2              frefp(m),frefs(m))
       h =(dble(rho(m)*b(m)*b(m))*atnb*atnb)
       gam=dble(b(m))*(wvno/om)
       gam = gam * atnb
       gam = dcmplx(2.0d+00,0.0d+00)*gam*gam
       gamm1 = gam - zone
       ra=CDSQRT(wvno2-xka*xka)
       rb=CDSQRT(wvno2-xkb*xkb)
       p=ra*dble(d(m))
       q=rb*dble(d(m))
       call var(p,q,ra,rb,w,x,y,z,cosp,cosq,
     1     ex,exa,exb,yl,zl,cosql,iwat)
         zrho = dcmplx(dble(rho(m)),0.0d+00)
C       WRITE(6,*)'chask: ex,exa,exb:',ex,exa,exb
c-----
c       For PSV exp(exa) is factored out of all
c       where exa is Re(ra h)
c       cosp = cosh(ra h)
c       w    = sinh(ra h)/ra
c       x    = ra sinh(ra h)
c       cosq = cosh(rb h)
c       y    = sinh(rb h)/rb
c       z    = rb sinh(rb h)
c       shorthand for some PSV combinations
c
c       The compount matrices need the following terms
c       which are in common ovrflw. Note that the 
c       exp(exa) is factored out where the exa is
c       Re ( (ra+rb)*h)
c       cpcq   = cosh(ra h) cosh(rb h )
c       cpy    = cosh(ra h) sinh (rb h)/rb
c       cpz    = cosh(ra h) rb sinh(rb h)
c       cqw    = cosh(rb h) sinh(ra h)/ra
c       cqx    = cosh(rb h) ra sinh(ra h)
c       xy     = ra sinh(ra h ) sinh(rb h )/rb
c       xz     = ra sinh(ra h ) rb sinh(rb h )
c       wy     = sinh(ra h)/ra  sinh (rb h)/rb
c       wz     = sinh(ra h)/ra  rb sinh(rb h )
c       a0      = exp(-exa)
c-----
      if(iwat.eq.1)then
c-----
c       fluid layer
c-----
        do i=1,6
           do j=1,6
              ca(i,j) = dcmplx(0.0d+00,0.0d+00)
           enddo
        enddo
        ca(1,1) = cosp
        ca(1,2) = - x /(zrho*om2)
        ca(2,1) = - w *zrho*om2
        ca(2,2) = cosp
        ca(3,3) = dcmplx(0.0d+00,0.0d+00)
        ca(4,4) = dcmplx(0.0d+00,0.0d+00)
        ca(5,5) = cosp
        ca(5,6) = - x /(zrho*om2)
        ca(6,5) = - w *zrho*om2
        ca(6,6) = cosp
      else
c-----
c       elastic layer
c-----
         zrho2= zrho*zrho
         gam2  = gam*gam
         gamm1 = gam-1.
         gamm2 = gamm1*gamm1
         cqww2 = cqw * wvno2
         cqxw2 = cqx / wvno2
         gg1 = gam*gamm1
         a0c  = dcmplx(2.0d+00,0.0d+00)*(dcmplx(a0,0.0d+00)-cpcq)
         xz2  = xz/wvno2
         gxz2 = gam*xz2
         g2xz2 = gam2 * xz2
         a0cgg1 = a0c*(gam+gamm1)
         wy2  = wy*wvno2
         g2wy2 = gamm2 * wy2
         g1wy2 = gamm1 * wy2
c-----
c       OK by symmetry
c----
         temp = a0c*gg1 + g2xz2 + g2wy2
         ca(1,1)=  cpcq-temp
         ca(1,2)=  (-cqx + wvno2*cpy)/(zrho*om2)
         temp = dcmplx(0.5d+00,0.0d+00)*a0cgg1 + gxz2 + g1wy2
         ca(1,3)=  wvno*temp/(zrho*om2)
         ca(1,4)= -ca(1,3)
         ca(1,5)=  (-cqww2+cpz)/(zrho*om2)
         temp = wvno2*(a0c + wy2) + xz
         ca(1,6)= -temp/(zrho2*om2*om2)
   
         ca(2,1)=  (-gamm2*cqw + gam2*cpz/wvno2)*zrho*om2
         ca(2,2)=  cpcq
         ca(2,3)=  (gamm1*cqww2 - gam*cpz)/wvno
         ca(2,4)= -ca(2,3)
         ca(2,5)= -wz
         ca(2,6)=  ca(1,5)
   
         temp =dcmplx(0.5d+00,0.0d+00)*a0cgg1*gg1 
     1          + gam2*gxz2 + gamm2*g1wy2
         ca(3,1)= -temp*zrho*om2/wvno
         ca(3,2)= -wvno*(gam*cqxw2 - gamm1*cpy)
         temp = a0c*gg1 + g2xz2 + g2wy2
         ca(3,3)=  a0 + temp
         ca(3,4)=  a0 - ca(3,3)
         ca(3,5)= -ca(2,3)
         ca(3,6)= -ca(1,3)
   
         ca(4,1)= -ca(3,1)
         ca(4,2)= -ca(3,2)
         ca(4,3)=  a0 - ca(3,3)
         ca(4,4)=  ca(3,3)
         ca(4,5)=  ca(2,3)
         ca(4,6)=  ca(1,3)
   
         ca(5,1)=  (-gam2*cqxw2 + gamm2*cpy)*zrho*om2
         ca(5,2)= -xy
         ca(5,3)= -ca(3,2)
         ca(5,4)=  ca(3,2)
         ca(5,5)=  ca(2,2)
         ca(5,6)=  ca(1,2)
   
         temp = gamm2*(a0c*gam2 + g2wy2) + gam2*g2xz2
         ca(6,1)= -zrho2*om2*om2*temp/wvno2
         ca(6,2)=  ca(5,1)
         ca(6,3)= -ca(3,1)
         ca(6,4)=  ca(3,1)
         ca(6,5)=  ca(2,1)
         ca(6,6)=  ca(1,1)
      endif
      return
      end

        subroutine chaski(m,ca,om,wvno,wvno2,om2,exa)
        implicit none
c-----
c       inverse reduced P-SV compound matrix
c       subroutine arguments
c-----
        integer m
        complex*16 ca(6,6)
        complex*16 om,wvno, wvno2, om2
        real*8 exa
c-----
c       rather than calculate a new matrix, use the property
c               -1
c       that C(A (h))  = C(a(-h))
c-----
        call chask(m,ca,om,wvno,wvno2,om2,exa)
        ca(1,2) = - ca(1,2)
        ca(1,5) = - ca(1,5)
        ca(2,1) = - ca(2,1)
        ca(2,3) = - ca(2,3)
        ca(2,4) = - ca(2,4)
        ca(2,6) = - ca(2,6)
        ca(3,2) = - ca(3,2)
        ca(3,5) = - ca(3,5)
        ca(4,2) = - ca(4,2)
        ca(4,5) = - ca(4,5)
        ca(5,1) = - ca(5,1)
        ca(5,3) = - ca(5,3)
        ca(5,4) = - ca(5,4)
        ca(5,6) = - ca(5,6)
        ca(6,2) = - ca(6,2)
        ca(6,5) = - ca(6,5)

        return
        end

        subroutine emat(m,wvno,om,om2,wvno2,epsv,esh)
c-----
c       evaluate E matrix for the halfspace
c-----
        implicit none
        integer m
        complex*16 wvno,om,om2,wvno2,epsv(4,4),esh(2,2)

c-----
c       global veriables
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b,rho,qa,qb,etap,etas,frefp,frefs
        common/modlly/mmax
        integer mmax
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
c-----
c       internal veriables
c-----
        complex*16 xka,xkb,ra,rb,gam,gamm1
        complex*16 atna,atnb
        complex*16 h
        integer iwat

        call aten(om,qa(m),qb(m),xka,xkb,alpha,a(m),
     1      b(m),atna,atnb,iwat,
     2      frefp(m),frefs(m))
        ra=CDSQRT(wvno2-xka*xka)
        rb=CDSQRT(wvno2-xkb*xkb)
        gam = dble(b(m))*wvno/om
        gam = gam * atnb
        gam = 2.0d+00 * (gam * gam)
        gamm1 = gam - dcmplx(1.0d+00,0.0d+00)
        h =(dble(rho(m)*b(m)*b(m))*atnb*atnb)
c-----
c       P-SV (8.4.4)
c-----
        epsv(1,1) =  wvno
        epsv(1,2) =  rb
        epsv(1,3) =  wvno
        epsv(1,4) = -rb
        epsv(2,1) =  ra
        epsv(2,2) =  wvno
        epsv(2,3) = -ra
        epsv(2,4) =  wvno
        epsv(3,1) =  rho(m)*om2*gamm1
        epsv(3,2) =  rho(m)*om2*gam*rb/wvno
        epsv(3,3) =  rho(m)*om2*gamm1
        epsv(3,4) = -rho(m)*om2*gam*rb/wvno
        epsv(4,1) =  rho(m)*om2*gam*ra/wvno
        epsv(4,2) =  rho(m)*om2*gamm1
        epsv(4,3) = -rho(m)*om2*gam*ra/wvno
        epsv(4,4) =  rho(m)*om2*gamm1
c-----
c       SH (8.3.2)
c-----
        esh(1,1) = wvno
        esh(1,2) = wvno
        esh(2,1) =  h*wvno*rb
        esh(2,2) = -h*wvno*rb
        return
        end
        
        subroutine emati(m,wvno,om,om2,wvno2,epsvi,eshi)
c-----
c       evaluate inverse E matrix for the halfspace
c-----
        implicit none
        integer m
        complex*16 wvno,om,om2,wvno2,epsvi(4,4),eshi(2,2)
c-----
c       global variables
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b,rho,qa,qb,etap,etas,frefp,frefs
        common/modlly/mmax
        integer mmax
        common/damp/alpha,ieqex
        real alpha
        integer ieqex
c-----
c       internal veriables
c-----
        complex*16 xka,xkb,ra,rb,gam,gamm1
        complex*16 atna,atnb
        complex*16 h
        integer iwat

        call aten(om,qa(m),qb(m),xka,xkb,alpha,a(m),
     1      b(m),atna,atnb,iwat,
     2      frefp(m),frefs(m))
        ra=CDSQRT(wvno2-xka*xka)
        rb=CDSQRT(wvno2-xkb*xkb)
        gam = dble(b(m))*wvno/om
        gam = gam * atnb
        gam = 2.0d+00 * (gam * gam)
        gamm1 = gam - dcmplx(1.0d+00,0.0d+00)
        h =(dble(rho(m)*b(m)*b(m))*atnb*atnb)
c-----
c       P-SV (8.4.5)
c-----
        epsvi(1,1) =  0.5*gam/wvno
        epsvi(1,2) = -0.5*gamm1/ra
        epsvi(1,3) = -0.5/(rho(m)*om2)
        epsvi(1,4) =  0.5*wvno/(rho(m)*om2*ra)

        epsvi(2,1) = -0.5*gamm1/rb
        epsvi(2,2) =  0.5*gam/wvno
        epsvi(2,3) =  0.5*wvno/(rb*rho(m)*om2)
        epsvi(2,4) = -0.5/(rho(m)*om2)

        epsvi(3,1) =  0.5*gam/wvno
        epsvi(3,2) =  0.5*gamm1/ra
        epsvi(3,3) = -0.5/(rho(m)*om2)
        epsvi(3,4) = -0.5*wvno/(rho(m)*om2*ra)

        epsvi(4,1) =  0.5*gamm1/rb
        epsvi(4,2) =  0.5*gam/wvno
        epsvi(4,3) = -0.5*wvno/(rb*rho(m)*om2)
        epsvi(4,4) = -0.5/(rho(m)*om2)
c-----
c       SH (8.3.2)
c-----
        eshi(1,1) =  dcmplx(0.5d+00, 0.0d+00)/wvno
        eshi(1,2) =  dcmplx(1.0d+00,0.0d+00)/(2.*h*wvno*rb)
        eshi(2,1) =  dcmplx(0.5d+00, 0.0d+00)/wvno
        eshi(2,2) = -dcmplx(1.0d+00,0.0d+00)/(2.*h*wvno*rb)
        return
        end

        subroutine var(p,q,ra,rb,w,x,y,z,cosp,cosq,ex,
     1      exa,exl,yl,zl,cosql,iwat)
c     not modified for negative p,q
c     this assumes that real p and real q have same signs
        common/ovrflw/a0,cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz
        complex*16 cpcq,cpy,cpz,cqw,cqx,xy,xz,wy,wz
        complex*16 p,q,ra,rb,w,x,y,z,cosp,cosq
        complex*16 yl,zl,cosql
        complex*16 eqp,eqm,epp,epm,sinp,sinq
        real *8 a0,pr,pi,qr,qi,fac,qmp,ex,exa,exl
c-----
c       form terms such as cos(p), sin(p), cos(q), sin(q)
c       and cos(p)*cos(q)
c
c       Introduce a factorization of exponentials to
c       make a pseudo floating point system
c
c       ex is the exponent in cosp
c       exl is the exponent in cosq for SH
c       exa is the exponent in cosp*cosq
c-----
        real*8 DREAL
      ex=0.0d+00
      exl = 0.0d+00
      a0=0.0d+00
      pr=dreal(p)
      pi=dimag(p)
      epp=dcmplx(dcos(pi),dsin(pi))/2.
      epm=dconjg(epp)
      ex=pr
      fac=0.0
      if(pr.lt.15.) fac=dexp(-2.*pr)
C     EX=0.0
C     FAC=1.0
      cosp=epp + fac*epm
      sinp=epp - fac*epm
      w=sinp/ra
      x=ra*sinp
        if(iwat.eq.1)then
c-----
c       fluid layer
c-----
            a0 = 1.0d+00
            exa = ex
            cosq = 1.0d+00
            y = 0.0d+00
            z = 0.0d+00
            cosql = 1.0d+00
            yl = 0.0d+00
            zl = 0.0d+00
            exl = 0.0d+00
        else
c-----
c       elastic layer
c-----
            qr=dreal(q)
            qi=dimag(q)
            eqp=dcmplx(dcos(qi),dsin(qi))/2.
            eqm=dconjg(eqp)
            exl=qr
            fac=0.0d+00
            if(qr.lt.15.) fac=dexp(-2.*qr)
C           EXL=0.0
C           FAC=1.0
            cosql=eqp + fac*eqm
            sinq=eqp - fac*eqm
            yl=sinq/rb
            zl=rb*sinq
c-----
c       form factors for compound P-SV matrix
c-----
            exa=pr + qr
C           EXA=0.0
            cpcq=cosp*cosql
            cpy=cosp*yl
            cpz=cosp*zl
            cqw=cosql*w
            cqx=cosql*x
            xy=x*yl
            xz=x*zl
            wy=w*yl
            wz=w*zl
            fac=0.0d+00
            qmp=qr-pr
            if(qmp.gt.-40.) fac=dexp(qmp)
C           QMP=0.0
C           FAC=1.0
            cosq=cosql*fac
            y=fac*yl
            z=fac*zl
            fac=0.0d+00
            if(exa.lt.60.) a0=dexp(-exa)
            FAC=1.0
        endif
        return
        end


        subroutine hask(m,om,om2,wvno,wvno2,apsv,ash,expsv,exsh)
c-----
c       evaluate propagator matrix for layer m  (Appendix A for P-SV and (8.3.3)for SH
c-----
        implicit none
        integer m
        complex*16 om,om2,wvno,wvno2,apsv(4,4),ash(2,2)
        real*8 expsv, exsh
c-----
c       global veriables
c-----
        integer NL
        parameter(NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL),
     2      frefp(NL), frefs(NL)
        real d,a,b,rho,qa,qb,etap,etas,frefp,frefs
        common/modlly/mmax
        integer mmax
        common/damp/alpha,ieqex
        real alpha
        integer ieqex

        common/lwater/lfluid
        logical lfluid
c-----
c       internal veriables
c-----
        complex*16 xka,xkb,ra,rb,gam,gamm1
        complex*16 atna,atnb
        complex*16 h
        integer iwat
        complex*16 p,q
        real*8 exa, exb, ex
        complex*16 w,x,y,z,cosp,cosq,yl,zl,cosql
        complex*16 cpq, gcpq, zw2, gzw2, g1w, g1y, gx
        real*8 zrho
        integer i,j

        call aten(om,qa(m),qb(m),xka,xkb,alpha,a(m),
     1      b(m),atna,atnb,iwat,
     2      frefp(m),frefs(m))
        ra=CDSQRT(wvno2-xka*xka)
        rb=CDSQRT(wvno2-xkb*xkb)
        gam = dble(b(m))*wvno/om
        gam = gam * atnb
        gam = 2.0d+00 * (gam * gam)
        gamm1 = gam - dcmplx(1.0d+00,0.0d+00)
        h =(dble(rho(m)*b(m)*b(m))*atnb*atnb)
        p=ra*dble(d(m))
        q=rb*dble(d(m))
        zrho = rho(m)
        call var(p,q,ra,rb,w,x,y,z,cosp,cosq,
     1              ex,exa,exb,yl,zl,cosql,iwat)
c-----
c       For PSV exp(exa) is factored out of all
c       where exa is Re(ra h)
c       cosp = cosh(ra h)
c       w    = sinh(ra h)/ra
c       x    = ra sinh(ra h)
c       cosq = cosh(rb h)
c       y    = sinh(rb h)/rb
c       z    = rb sinh(rb h)
c       shorthand for some PSV combinations
c
c       For SH exp(exl) is factored out
c       cosql = cosh(rh h)
c       yl    = sinh(rb h)/rb
c       zl    = rb sinh(rb h)
c
c-----
        if(iwat.eq.1)then
c           fluid
            do i=1,4
              do j=1,4
                 apsv(i,j) = dcmplx(0.0d+00,0.0d+00)
              enddo
            enddo
            apsv(1,1) =   dcmplx(1.0d+00,0.0d+00)
            apsv(4,4) =   dcmplx(1.0d+00,0.0d+00)
            apsv(2,2) =   cosp
            apsv(3,3) =   cosp
            apsv(2,3) = - x / (zrho * om2 )
            apsv(3,2) = - zrho * om2 * w
         
            ash(1,1) =    dcmplx(1.0d+00,0.0d+00)
            ash(2,1) =    dcmplx(0.0d+00,0.0d+00)
            ash(1,2) =    dcmplx(0.0d+00,0.0d+00)
            ash(2,2) =    dcmplx(1.0d+00,0.0d+00)
            expsv = ex
            exsh  = exb
        else
c           solid
            expsv = ex
            cpq = cosp-cosq
            gcpq = gam*cpq
            zw2 = z/wvno2
            gzw2 = gam*zw2
            g1w = gamm1*w
            g1y = gamm1*y
            gx = gam*x
            apsv(1,1)=   gcpq + cosq
            apsv(1,2)=   wvno*(-g1w+gzw2)
            apsv(1,3)= - wvno * cpq/(zrho*om2)
            apsv(1,4)=   (wvno2*w-z)/(zrho*om2)
            apsv(2,1)=   (gx - wvno2*g1y)/wvno
            apsv(2,2)= - gcpq + cosp
            apsv(2,3)=   (-x+wvno2*y)/(zrho*om2)
            apsv(2,4)= - apsv(1,3)
            apsv(3,1)=   zrho*om2*gamm1*gcpq/wvno
            apsv(3,2)=   zrho*om2*((-gamm1*g1w)+(gam*gzw2))
            apsv(3,3)=   apsv(2,2)
            apsv(3,4)= - apsv(1,2)
            apsv(4,1)=   zrho*om2*(((gam*gx)/wvno2) - (gamm1*g1y))
            apsv(4,2)= - apsv(3,1)
            apsv(4,3)= - apsv(2,1)
            apsv(4,4)=   apsv(1,1)

            exsh  = exb
            ash(1,1) = cosql
            ash(2,1) = zl*h
            ash(1,2) = yl/h
            ash(2,2) = cosql
        endif

        return
        end
        subroutine haski(m,om,om2,wvno,wvno2,apsvi,ashi,expsv,exsh)
c-----
c       evaluate inverse propagator matrix for layer m  (Appendix A for P-SV and (8.3.3)for SH
c       this is easily done by changin the sign of all sinh terms in the hask routine
c-----
        implicit none
        integer m
        complex*16 om,om2,wvno,wvno2,apsvi(4,4),ashi(2,2)
        real*8 expsv, exsh

        call hask(m,om,om2,wvno,wvno2,apsvi,ashi,expsv,exsh)
c-----
c       change the sign of the odd i+j
c-----
        apsvi(1,2) = - apsvi(1,2)
        apsvi(1,4) = - apsvi(1,4)
        apsvi(2,1) = - apsvi(2,1)
        apsvi(2,3) = - apsvi(2,3)
        apsvi(3,2) = - apsvi(3,2)
        apsvi(3,4) = - apsvi(3,4)
        apsvi(4,1) = - apsvi(4,1)
        apsvi(4,3) = - apsvi(4,3)

        ashi(1,2) = - ashi(1,2)
        ashi(2,1) = - ashi(2,1)
        return
        end
c-----
c      simple matrix routines
c-----
      subroutine  copyab2(a,b)
c-----
c     copy 2x2 matrix a to matrix b
c-----
      implicit none
      integer DIM
      parameter (DIM=2)
c-----
c     procedure arguments
c-----
      complex*16 a(DIM,DIM), b(DIM,DIM)
c-----
c     internal variables
c-----
      integer i,j
      do j=1,DIM
           do i=1,DIM
              b(i,j) = a(i,j)
           enddo
      enddo
      return
      end

      subroutine  copyab4(a,b)
c-----
c     copy 4x4 matrix a to matrix b
c-----
      implicit none
      integer DIM
      parameter (DIM=4)
c-----
c     procedure arguments
c-----
      complex*16 a(DIM,DIM), b(DIM,DIM)
c-----
c     internal variables
c-----
      integer i,j
      do j=1,DIM
           do i=1,DIM
              b(i,j) = a(i,j)
           enddo
      enddo
      return
      end

      subroutine  copyab6(a,b)
c-----
c     copy 6x6 matrix a to matrix b
c-----
      implicit none
      integer DIM
      parameter (DIM=6)
c-----
c     procedure arguments
c-----
      complex*16 a(DIM,DIM), b(DIM,DIM)
c-----
c     internal variables
c-----
      integer i,j
      do j=1,DIM
           do i=1,DIM
              b(i,j) = a(i,j)
           enddo
      enddo
      return
      end

      subroutine  multabc2(a,b,c)
c-----
c     evaluate 2x2 multiplication c=a b
c-----
      implicit none
      integer DIM
      parameter (DIM=2)
c-----
c     procedure arguments
c-----
      complex*16 a(DIM,DIM), b(DIM,DIM), c(DIM,DIM)
c-----
c     internal variables
c-----
      integer i,j,k
      do j=1,DIM
           do i=1,DIM
                c(i,j) = dcmplx(0.0d+00,0.0d+00)
                do k=1,DIM
                     c(i,j) = c(i,j) + a(i,k)*b(k,j)
                enddo
           enddo
      enddo
      return
      end

      subroutine  multabc4(a,b,c)
c-----
c     evaluate 4x4 multiplication c=a b
c-----
      implicit none
      integer DIM
      parameter (DIM=4)
c-----
c     procedure arguments
c-----
      complex*16 a(DIM,DIM), b(DIM,DIM), c(DIM,DIM)
c-----
c     internal variables
c-----
      integer i,j,k
      do j=1,DIM
           do i=1,DIM
                c(i,j) = dcmplx(0.0d+00,0.0d+00)
                do k=1,DIM
                     c(i,j) = c(i,j) + a(i,k)*b(k,j)
                enddo
           enddo
      enddo
      return
      end

      subroutine  multabc6(a,b,c)
c-----
c     evaluate 6x6 multiplication c=a b
c-----
      implicit none
      integer DIM
      parameter (DIM=6)
c-----
c     procedure arguments
c-----
      complex*16 a(DIM,DIM), b(DIM,DIM), c(DIM,DIM)
c-----
c     internal variables
c-----
      integer i,j,k
      do j=1,DIM
           do i=1,DIM
                c(i,j) = dcmplx(0.0d+00,0.0d+00)
                do k=1,DIM
                     c(i,j) = c(i,j) + a(i,k)*b(k,j)
                enddo
           enddo
      enddo
      return
      end


      subroutine cmultrow6(a,b,c)
      implicit none
      integer DIM
      parameter (DIM=6)
c-----
c     form c = a b where
c        a (1x6), a(6x6) and C(1X6)
C-----
c     subroutine arguments
c-----
       complex*16 a(1,DIM), b(DIM,DIM), c(1,DIM)
c-----
c      internal parameters
c-----
       integer i,j
       do i=1,DIM
          c(1,i) = dcmplx(0.0d+00,0.0d+00)
          do j=1,DIM
              c(1,i) = c(1,i) + a(1,j)*b(j,i)
          enddo
       enddo
       return
       end

      subroutine cmultcol6(a,b,c)
      implicit none
      integer DIM
      parameter (DIM=6)
c-----
c     form c = a b where
c        a (6x6), a(bx1) and C(6X1)
C-----
c     subroutine arguments
c-----
       complex*16 a(DIM,DIM), b(DIM,1), c(DIM,1)
c-----
c      internal parameters
c-----
       integer i,j
       do i=1,DIM
          c(i,1) = dcmplx(0.0d+00,0.0d+00)
          do j=1,DIM
              c(i,1) = c(i,1) + a(i,j)*b(j,1)
          enddo
       enddo
       return
       end



        subroutine identity4(a)
c-----
c       create a 4x4 identity matrix
c-----
        implicit none
        integer NDIM
        parameter (NDIM=4)
        complex*16 a(NDIM,NDIM)

        integer i,j
        do i=1,NDIM
           do j=1,NDIM
              if(i.ne.j)then
                  a(i,j) = dcmplx(0.0d+00,0.0d+00)
              else
                  a(i,j) = dcmplx(1.0d+00,0.0d+00)
              endif
           enddo
         enddo
         return
         end

        subroutine identity6(a)
c-----
c       create a 6x6 identity matrix
c-----
        implicit none
        integer NDIM
        parameter (NDIM=6)
        complex*16 a(NDIM,NDIM)

        integer i,j
        do i=1,NDIM
           do j=1,NDIM
              if(i.ne.j)then
                  a(i,j) = dcmplx(0.0d+00,0.0d+00)
              else
                  a(i,j) = dcmplx(1.0d+00,0.0d+00)
              endif
           enddo
         enddo
         return
         end

        subroutine identity2(a)
c-----
c       create a 2x2 identity matrix
c-----
        implicit none
        integer NDIM
        parameter (NDIM=2)
        complex*16 a(NDIM,NDIM)

        integer i,j
        do i=1,NDIM
           do j=1,NDIM
              if(i.ne.j)then
                  a(i,j) = dcmplx(0.0d+00,0.0d+00)
              else
                  a(i,j) = dcmplx(1.0d+00,0.0d+00)
              endif
           enddo
         enddo
         return
         end
      subroutine htrav(x,hs,hr,dop,tfirst,rayp,tstar)
      implicit none
      real x,hs, hr,tfirst,rayp,tstar
      logical dop
c-----
c     read the hudsonsrc.mod and hudsonrec.mod
c     and calculate first arrival travel times for
c     refractions
c     It is assumed that both models have the same number
c     of layers and that at depth they are identical. This
c     is ensured by the way these files are created. 
c     However to ensure this the source and receiver depths must be
c     placed into each model
c
c     also comptue the total tstar
c     Note the tak135sph.mod approximated the core as a fluid with S velocity.
c     There is no inner core.  Thus when computing arrivals, another
c     constraint is that the bottom velocity must be greater than
c     all of  the velocitie from depths hs inlthe source model
c     and hr in the receiver model to that depth
c-----
        integer NL
        parameter (NL=200)
        common/isomod/d(NL),a(NL),b(NL),rho(NL),
     1      qa(NL),qb(NL),etap(NL),etas(NL), 
     2      frefp(NL), frefs(NL)
        real d, a, b, rho, qa, qb, etap, etas, frefp, frefs
        common/modlly/mmmax
        integer mmmax
        common/depref/refdep
        real refdep
        common/shwave/bsh(NL), rhosh(NL), qbsh(NL)
        real bsh, rhosh, qbsh

        integer mmax

        integer iunit, iiso, iflsph, idimen, icnvel, ierr
        character title*80

        common/earth/radius
        real radius
c-----
c       local variables
c-----
        double precision tds, tdr
        double precision tss, tsr
        real dphs, dphr
        integer lmaxs, lmaxr, lmax, lmin
        double precision ts
        integer i,m
        real vraymx, vavg, qiavg
        double precision tints, tintr, time
        double precision s, c
        integer bottom
        logical docheckok
        double precision sumxs, sumxr
        double precision tsssv,tsrsv,sumxrsv,sumxssv

        real h(NL),v(2,NL),qi(2,NL)
c-----
c       initialize
c-----
        radius = 6371.0
        tstar = 0.0
c-----
c       get model on source side
c-----
             call getmod(1,'hudsonsrc.mod',mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
             mmmax = mmax
             call insert(hs)
             call insert(hr)
             call srclyr(hs,lmaxs,dphs)
             tds = radius*alog(radius/(radius-hs))
             call adosph()
c-----
c            save the source side
c-----
             do i=1,mmax
                  h(i) = d(i)
                  if(dop)then
                        v(1,i)  = a(i)
                        qi(1,i) = qa(i)
                        if(qi(1,i).gt.1.0)then
                           qi(1,i) = 1.0/qi(1,i)
                        endif
                  else
                        if(b(i).lt.0.001*a(i))then
c-----
c                            In fluid layer permit P is fluid as
c                           part of S wavefield
c-----
                             v(1,i)  = a(i)
                             qi(1,i) = qa(i)
                        else
                             v(1,i)  = b(i)
                             qi(1,i) = qb(i)
                        endif
                        if(qi(1,i).gt.1.0)then
                           qi(1,i) = 1.0/qi(1,i)
                        endif
                  endif
             enddo
c-----
c       get model on receiver side
c-----
             call getmod(1,'hudsonrec.mod',mmax,title,iunit,iiso,iflsph,
     1          idimen,icnvel,ierr,.false.)
             mmmax = mmax
             call insert(hs)
             call insert(hr)
             call srclyr(hr,lmaxr,dphr)
             tdr = radius*alog(radius/(radius-hr))
             call adosph()
c-----
c            save the receiver side
c-----
             do i=1,mmax
                  if(dop)then
                        v(2,i)  = a(i)
                        qi(2,i) = qa(i)
                        if(qi(2,i).gt.1.0)then
                           qi(2,i) = 1.0/qi(2,i)
                        endif
                  else
                        if(b(i).lt.0.001*a(i))then
c-----
c                            In fluid layer permit P is fluid as
c                           part of S wavefield
c-----
                             v(2,i)  = a(i)
                             qi(2,i) = qa(i)
                        else
                             v(2,i) =  b(i)
                             qi(2,i) = qb(i)
                        endif
                        if(qi(2,i).gt.1.0)then
                           qi(2,i) = 1.0/qi(2,i)
                        endif
                  endif
             enddo
c-----
c        Calculate the downward refraction time
c        Start with properties of deepest of source and receiver depths, e.g.,
c        lmax = max0(lmaxs,lmaxr)
c        Then loop over possible refractors
c        always check that the possible refractor velocity is 
c        greater than the layer velocities at source or receiver
c        vraymax=max1(v(1,laxms),v(2,lmaxr))
c-----
        vraymx = max1(v(1,lmaxs),v(2,lmaxr))
        lmax   = max0(lmaxs,lmaxr)
        lmin   = min0(lmaxs,lmaxr)
c-----
c       here the refraction is at the top of layer m
c       First check if this can be a refractor with velocity > vraymx
c       Then compute the shortest distance for a refraction,
c          which is SUM h_i tan \theta_i down to the refractor
c          and require x >= SUM h_i tan \theta_i 
c       If this passes, then
c          t = SUM h_i cos \theta_i / V_i + x /V_m
c       Then save the minimum time
c
c       Because it may be possible that the v(1,lmax) != v(2,lmax)
c          HOPE that the minimum time is due to a deeper layer.
c          This is a limitation of this approach
c          
c-----
        tfirst = 1.0e+30
        do m=lmax+1,mmax
            vavg = 0.5*(v(1,m)+v(2,m))
            qiavg = 0.5*(qi(1,m) + qi(2,m))
c-----
c           loop over layers above
c----- 
            tintr = 0.0d+00
            tints = 0.0d+00
            sumxr = 0.0
            sumxs = 0.0
            if(docheckok(vavg,v,lmin,m))then
               tss = 0.0
               tsr = 0.0
               do i=lmin,m-1
c-----
c              all layer velocities down to m-1 must be
c              less than v(m)
c-----
                  if(i.ge.lmaxr)then
                     s = v(2,i)/vavg
                     c = dsqrt(1.0+00-s*s)
                     tintr = tintr + h(i)*c/v(2,i)
                     tsr   = tsr   + qi(2,i)*h(i)/(c*v(2,i))
                     sumxr = sumxr + h(i) * s/c
       
                  endif
                  if(i.ge.lmaxs)then
                     s = v(1,i)/vavg
                     c = dsqrt(1.0+00-s*s)
                     tints = tints + h(i)*c/v(1,i)
                     tss   = tss   + qi(1,i)*h(i)/(c*v(1,i))
                     sumxs = sumxs + h(i) * s/c
                  endif
               enddo
           
c-----
c              time computed using refraction rule. The horizontal offset
c              of the ray from the upper layering is icorporated into
c              the intercept times
c-----
               time = tints + tintr + x/vavg
               if(time.lt.tfirst.and.x.gt.(sumxr+sumxs))then
                   tfirst = sngl(time )
                   rayp = 1.0/vavg
c-----
c              to compute T* we need the actual length along the halfspace
c              and thus must account for the horizontal offsets of the layers
c-----
                ts = tss + tsr  + qiavg*(x-sumxs-sumxr)/vavg
                   tstar = sngl(ts)
                   bottom = m
                tsssv = tss
                tsrsv = tsr
                sumxrsv = sumxr
                sumxssv = sumxs
               endif
            endif
        enddo

        return
        end

        logical function docheckok(vavg,v,lmin,m)
        implicit none
        real vavg
        integer lmin,m
        integer NL
        parameter (NL=200)
        real v(2,NL)

        integer i
        docheckok = .true.
        do i=lmin,m-1
           if(vavg.lt.v(1,i) .or. vav g.lt.v(2,i))then
               docheckok = .false.
           endif
        enddo
        return
        end
        subroutine firstarr(r,hs,hr,dop,time,vprec,vsrec,denrec,
     1      vpsrc, vssrc, densrc, rayp, geom, tstar, dogeom)
        implicit none
c-----
c       r   R   Epicentral distance
c       hs  R   Source depth
c       hr  R   Receiver depth
c       dop L   .true. if P or .false. if S
c       time    R   First arrival time
c       vprec    R   Velocity of P wave at receiver
c       vsrec    R   Velocity of S wave at receiver
c       denrec     R   Density at receiver
c       vpsrc R   P-wave velocity at source
c       vssrc R   S-wave velocity at source
c       densrc R   Density at source
c       rayp R   Ray parameter in sec/km
c       geom R   geometrical spreading factor
c       tstar R   geometrical spreading factor
c       dogeom L .true. compute geometrical spreading. Only do this for
c                teleseisms
c-----
        real r, hs, hr, time, vprec, vsrec,denrec, vpsrc, vssrc, densrc
        real rayp, geom, tstar
        logical  dogeom
c-----
c-----
c       internal variables
c-----
        real tp, tm, rp(2)
        real tstarp, tstarm
        real dpdx, deginrad, sinis, sinir, cosir, cosis, sindeg
        real vs, vr, rs, rr
        real dr, fac
        logical dop

        
        common/earth/radius
        real radius

        common/depref/refdep
        real refdep

        common/control/verbose
        logical verbose

c-----
c       compute the travel time
c-----
        call htrav(r,hs,hr,dop,time,rayp,tstar)

c-----
C       compute the geometrical spreading
C
C       since a flat earth is always used, we use the flat layered
C       geometrical spreading from Officer
C       The geometrical spreading is dimensionless and gives the decrease in amplitude from
C       a distance of 1 km from the source to the receiver
C                        2                            2
C            ( rhos Vs  a sin Is Vs                  d T      )
C       sqrt |  -------------------------------     -----     |
C            |                                         2      |
C            ( rhor Vr sin DEG  cos Ir Rs Cos Is     dx       )
C
C       which is dimensionally 1/km
C
C       where p (sec/km) = dT/dx
C       a = radius of sphere about source - we use 1.0 km
C       Is= incident angle at source
C       Ir= incident angle at receiver
C       Rs= distance from center of Earth to source
C       Rr= distance from center of Earth to receiver
C       DEG=epicental distance in degrees
C       rhos and vs are the density and wave velocity at the source depth
c       rhor and vr are the density and wave velocity at the receiver depth
c
c       To get the dp/dx, we determine p at different distances, e.g.,
c       at r -dr and r + dr , and then form
c       an interpolating polynomial
c
        dr = 500.0

     
        if(dogeom)then
        call htrav(r-dr,hs,hr,dop,tm,rp(1),tstarm)
        call htrav(r+dr,hs,hr,dop,tp,rp(2),tstarp)


              dpdx = abs(rp(1) - rp(2))/(dr - ( - dr))
              deginrad = r/radius
              sindeg = sin(deginrad)
  
              if(dop)then
                  vs = vpsrc
                  vr = vprec
                  rs = densrc
                  rr = denrec
              else 
                  vs = vssrc
                  vr = vsrec
                  rs = densrc
                  rr = denrec
              endif
              sinis = rayp*vs
              cosis = sqrt(1.0-sinis*sinis)
              sinir = rayp*vr
              cosir = sqrt(1.0-sinir*sinir)

              fac = (rs*vs*sinis*vs*dpdx)/
     1           (rr*vr*sindeg*cosir*(radius-hs)*cosis)
              geom = sqrt(abs(fac))
C           if(verbose)then
C                 WRITE(0,*)'Subroutine frstar:'
C                 WRITE(0,*)'rs     :',rs
C                 WRITE(0,*)'vs     :',vs
C                 WRITE(0,*)'hs     :',hs
C                 WRITE(0,*)'rr     :',rr
C                 WRITE(0,*)'vr     :',vr
C                 WRITE(0,*)'sinis  :',sinis
C                 WRITE(0,*)'cosis  :',cosis
C                 WRITE(0,*)'cosir  :',cosir
C                 WRITE(0,*)'rp(1)  :',rp(1)
C                 WRITE(0,*)'rp(2)  :',rp(2)
C                 WRITE(0,*)'dr     :',dr
C                 WRITE(0,*)'dpdx   :',dpdx
C                 WRITE(0,*)'sindeg :',sindeg
C                 WRITE(0,*)'radius :',radius
C                 WRITE(0,*)'geom  :',geom 
C           endif
      else
c-----
c     default to have a defined return value
c-----
              geom = 1.0
      endif
        call htrav(r,hs,hr,dop,time,rayp,tstar)
        return
        end

      subroutine anorm2(a,ex)
      implicit none
c-----
c     normalize the propagator matrix
c     and increment the exponent
c-----
      integer DIM
      parameter (DIM=2)
c-----
c     subroutine variables
c-----
      complex*16 a(DIM,DIM)
      real*8 ex

      real*8 test,testt,x,y,fac,xnorm

c-----
c     internal variables
c-----
      integer i,j

c-----
c     carefully determine the size of the
c     matrix by identifying the largest value
c     in terms of absolute value
c
c     To get absolute value of z = x + i y
c     first identify the largest of all x and y
c     as testt
c-----
        test = 0.0D+00
        testt = 0.0D+00
        do  j=1,DIM
            do  i = 1,DIM
            if(dabs(dreal(a(i,j))).gt.testt)
     1          testt =dabs(dreal(a(i,j)))
            if(dabs(dimag(a(i,j))).gt.testt)
     1          testt =dabs(dimag(a(i,j)))
            enddo
        enddo
c-----
c      instead of forming |z|^2 = x^2 + y^2
c      divide by the testt (never zero) to have a set of numbers
c      near unity. Save the largest |z|^2
c-----
        if(testt.lt.1.0e-30)testt=1.0
        do  j=1,DIM
            do  i =1,DIM
                x=dreal(a(i,j))/testt
                y=dimag(a(i,j))/testt
                fac = x*x + y*y
                if(test.lt.fac) test = fac
            enddo
        enddo
        test = testt*dsqrt(test)
c-----
c       test is the normalization factor but it must never be
c       zero
c-----
        if(test.lt.1.0d-30) test=1.0
        xnorm = 1./test
c-----
c     apply the normalization to the matrix
c-----
      do i=1,DIM
         do j=1,DIM
            a(i,j) = a(i,j)*xnorm
         enddo
      enddo
      ex = ex + dlog(test)
      return
      end

  
      subroutine anorm4(a,ex)
      implicit none
c-----
c     normalize the propagator matrix
c     and increment the exponent
c-----
      integer DIM
      parameter (DIM=4)
c-----
c     subroutine variables
c-----
      complex*16 a(DIM,DIM)
      real*8 ex

      real*8 test,testt,x,y,fac,xnorm

c-----
c     internal variables
c-----
      integer i,j

c-----
c     carefully determine the size of the
c     matrix by identifying the largest value
c     in terms of absolute value
c
c     To get absolute value of z = x + i y
c     first identify the largest of all x and y
c     as testt
c-----
        test = 0.0D+00
        testt = 0.0D+00
        do  j=1,DIM
            do  i = 1,DIM
            if(dabs(dreal(a(i,j))).gt.testt)
     1          testt =dabs(dreal(a(i,j)))
            if(dabs(dimag(a(i,j))).gt.testt)
     1          testt =dabs(dimag(a(i,j)))
            enddo
        enddo
c-----
c      instead of forming |z|^2 = x^2 + y^2
c      divide by the testt (never zero) to have a set of numbers
c      near unity. Save the largest |z|^2
c-----
        if(testt.lt.1.0e-30)testt=1.0
        do  j=1,DIM
            do  i =1,DIM
                x=dreal(a(i,j))/testt
                y=dimag(a(i,j))/testt
                fac = x*x + y*y
                if(test.lt.fac) test = fac
            enddo
        enddo
        test = testt*dsqrt(test)
c-----
c       test is the normalization factor but it must never be
c       zero
c-----
        if(test.lt.1.0d-30) test=1.0
        xnorm = 1./test
c-----
c     apply the normalization to the matrix
c-----
      do i=1,DIM
         do j=1,DIM
            a(i,j) = a(i,j)*xnorm
         enddo
      enddo
      ex = ex + dlog(test)
      return
      end

  
  
      subroutine anorm6(a,ex)
      implicit none
c-----
c     normalize the propagator matrix
c     and increment the exponent
c-----
      integer DIM
      parameter (DIM=6)
c-----
c     subroutine variables
c-----
      complex*16 a(DIM,DIM)
      real*8 ex

      real*8 test,testt,x,y,fac,xnorm

c-----
c     internal variables
c-----
      integer i,j

c-----
c     carefully determine the size of the
c     matrix by identifying the largest value
c     in terms of absolute value
c
c     To get absolute value of z = x + i y
c     first identify the largest of all x and y
c     as testt
c-----
        test = 0.0D+00
        testt = 0.0D+00
        do  j=1,DIM
            do  i = 1,DIM
            if(dabs(dreal(a(i,j))).gt.testt)
     1          testt =dabs(dreal(a(i,j)))
            if(dabs(dimag(a(i,j))).gt.testt)
     1          testt =dabs(dimag(a(i,j)))
            enddo
        enddo
c-----
c      instead of forming |z|^2 = x^2 + y^2
c      divide by the testt (never zero) to have a set of numbers
c      near unity. Save the largest |z|^2
c-----
        if(testt.lt.1.0e-30)testt=1.0
        do  j=1,DIM
            do  i =1,DIM
                x=dreal(a(i,j))/testt
                y=dimag(a(i,j))/testt
                fac = x*x + y*y
                if(test.lt.fac) test = fac
            enddo
        enddo
        test = testt*dsqrt(test)
c-----
c       test is the normalization factor but it must never be
c       zero
c-----
        if(test.lt.1.0d-30) test=1.0
        xnorm = 1./test
c-----
c     apply the normalization to the matrix
c-----
      do i=1,DIM
         do j=1,DIM
            a(i,j) = a(i,j)*xnorm
         enddo
      enddo
      ex = ex + dlog(test)
      return
      end

  
