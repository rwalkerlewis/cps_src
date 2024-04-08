        subroutine gtrofa(rho,vp)
        implicit none
        real rho, vp
c-----
c       rho - density in gm/cm^3
c       vp  - P-wave velocity in km/s
c-----
c    CHANGES:
c      15 AUG 2013 - incremented arrays to be able use Vp up to 12 km/sec
c               this is essentially at depth of about 900 km 
c      29 JUL 2022 Changed to use
c           Brocher,, T.M. (2005). Empirical relations 
c           between elastic wavesppeds and density in the Earth's crust, 
c           Bull. Seism.Soc.Am, 95, 2081-2092.
c           Note: Brocher extended past Vp=8.5 Vs=4.7 to handle upper
c           mantle Vp upto 12 (900 km) to roughly fit values in PREM
c-----
        if(vp.le.8.5)then
           rho = vp*(1.6612+vp*
     1     (-0.4721+vp*(0.0671+vp*(-0.0043+vp*0.000106))))
        else
           rho = 3.45 + ( vp - 8.5)*1.32 /3.5
        endif
        return
        end

        subroutine gtaofb(vp,vs)
        implicit none
        real vs, vp
c-----
c       vp  - P-wave velocity in km/s
c       vs  - S-wave velocity in km/s
c-----
c    CHANGES:
c      29 JUL 2022 Changed to use
c           Brocher,, T.M. (2005). Empirical relations 
c           between elastic wavesppeds and density in the Earth's crust, 
c           Bull. Seism.Soc.Am, 95, 2081-2092.
c           Note: for Vs > 4.84, e.g., Vp = 8.5, results are extrapolated using PREM
c             there is a problemwith this since the PREM includes effect of partial melp
c             in upper mantle
c-----
        if(vs.le.4.84)then
            vp=0.9408+vs*(2.094+vs*(-0.8206+vs*(0.2683+vs*(-0.0251))))
        else
            vp = 8.50 + (vs - 4.84)*3.50/1.83
        endif
        return
        end

     
