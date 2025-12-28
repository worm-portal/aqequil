      subroutine pcorrm(adh,adhh,adhv,al10,aphi,avcnst,bdh,bdhh,bdhv,
     $ bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,
     $ dvfe,iopg,ipch,ipchmx,ipcv,ipcvmx,nopgmx,presg,press,rcnstv,
     $ tempk,xhfe,xlke,xvfe)
c
c     This subroutine makes pressure corrections for miscellaneous
c     thermodynamic functions, such as Debye-Huckel parameters.
c     It normally corrects for pressures off the standard T-P grid
c     (e.g., from "presg" to "press"). However, it can be used to
c     correct from any pressure to another by substituting the value
c     of the former for "presg". This is done is EQ6 to correct for
c     changing pressure along an isothermal reaction path. In order
c     to do this, the pressure derivatives of the parameters in
c     question are also corrected to the pressure of interest. The
c     derivatives of highest order are necessarily treated as constants,
c     so there is no correction in this case.
c
c     Note the following technical references:
c
c       Ananthaswamy, J., and Atkinson, G., 1984, Thermodynamics of
c         concentrated electrolyte mixtures. 4. Pitzer-Debye-Huckel
c         limiting slopes for water from 0 to 100 C and from 1 atm
c         to 1 kbar, Journal of Chemical and Engineering Data,
c         v. 29, p. 81-87.
c
c       Bradley, D. J., and Pitzer, K. S., 1979, Thermodynamics of
c         electrolytes. 12. Dielectric properties of waer and
c         Debye-Huckel parameters to 350 C and 1 kbar, Journal
c         of Physical Chemistry, v. 83, p. 1599-1603.
c
c       Helgeson, H. C., and Kirkham, D. H., 1974b, Theoretical
c         prediction of the thermodynamic behavior of aqueous
c         electrolytes at high pressures and temperatures: II.
c         Debye-Huckel parameters for activity coefficients and
c         relative partial molal properties, American Journal
c         of Science, v. 274, p. 1199-1261.
c
c     EQLIB/pcorrx.f performs the same function for the standard state
c     thermodynamic functions, such as the equilibrium constants.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c      avcnst = 2.303 RT, with R in units of bar-cm3/mol-K
c
c     Principal output:
c
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipchmx,ipcvmx,nopgmx
c
      integer iopg(nopgmx)
c
      integer ipch,ipcv
c
      real*8 dadhh(ipchmx),dadhv(ipcvmx),dbdhh(ipchmx),dbdhv(ipcvmx),
     $ dbdth(ipchmx),dbdtv(ipcvmx),dhfe(ipchmx),dvfe(ipcvmx)
c
      real*8 adh,adhh,adhv,al10,aphi,avcnst,bdh,bdhh,bdhv,bdot,bdoth,
     $ bdotv,presg,press,rcnstv,tempk,xhfe,xlke,xvfe
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ipc,n
c
      real*8 a22rt,avrt,dp,xx
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
c     a22rt  = 2(2.303)RT (R in P-V units)
c     avrt   = RT (R in P-V units)
c
      a22rt = 2*avcnst
      avrt = rcnstv*tempk
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the pressure difference.
c
      dp = press - presg
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     This section addresses various parameters at the Gibbs energy
c     level (e.g., A(gamma,10), B(gamma)). It is divided into two
c     parts, one for first-order corrections only, the other for
c     corrections of order 2 or higher. This division is made for
c     speed, as often the data will only be available for first-order
c     corrections.
c
      if (ipcv .eq. 0) then
c
c       Make only first-order corrections.
c
c       The A(gamma,10) and A(phi) parameters. Note that A(V) (adhv)
c       is not dA(gamma,10)/dp nor dA(phi)/dp, though it does contain
c       these derivatives. Also, Helgeson and Kirkham's (1974b)
c       definition of A(V) is not identical to Bradley and Pitzer's
c       (1979). See the discussion of Ananthaswamy and Atkinson (1984).
c       Here if iopg(1) is less than or equal to 0, we follow Helgeson
c       and Kirkham's definition, otherwise, Bradley and Pitzer's.
c
        if (iopg(1) .le. 0) then
c
c         Helgeson and Kirkham (1974b), eq 48:
c
c           A(V) = -2(2.303)RT (dA(gamma,10)/dp)
c
c         Compare Ananthaswamy and Atkinson (1984), table 1:
c
c           A(V) = -6RT (dA(phi)/dp)
c
          xx = -adhv*dp/a22rt
          adh = adh + xx
          aphi = adh*al10/3.
        else
c
c         Ananthaswamy and Atkinson (1984), table 1:
c
c           A(V) = -4RT (dA(phi)/dp)
          xx = -adhv*dp/(4.*avrt)
          aphi = aphi + xx
          adh= 3.*aphi/al10
        endif
c
c       The B(gamma) parameter.
c
c         Helgeson and Kirkham (1974b), eq 49:
c
c           B(V) = 2(2.303)RT (dB(gamma)/dp)
c
c       Note that B(gamma) is carried multiplied by 10-8,
c       while B(V) is carried multiplied by 10-6.
c
        xx = 0.01*bdhv*dp/a22rt
        bdh = bdh + xx
c
c       The B-dot parameter.
c
c         by analogy to Helgeson and Kirkham (1974b), eq 49:
c
c           B-dot(V) = 2(2.303)RT (dB-dot/dp)
c
        xx = bdotv*dp/a22rt
        bdot = bdot + xx
c
c       The log K for the "Eh" reaction.
c
        xx = -xvfe*dp/avcnst
        xlke = xlke + xx
      endif
c
      if (ipcv .gt. 0) then
c
c       Make only corrections of order 2 or higher.
c
c       The A(gamma,10) and A(phi) parameters.
c
        xx = adhv*dp
        do ipc = 1,ipcv
          n = ipc + 1
          xx = xx - ( dadhv(ipc)*(dp**n) )/fctrl(n)
        enddo
        if (iopg(1) .le. 0) then
          xx = -xx/a22rt
          adh = adh + xx
          aphi = adh*al10/3.
        else
          xx = -xx/(4.*avrt)
          aphi = aphi + xx
          adh= 3.*aphi/al10
        endif
c
c       The B(gamma) parameter.
c
        xx = bdhv*dp
        do ipc = 1,ipcv
          n = ipc + 1
          xx = xx - ( dbdhv(ipc)*(dp**n) )/fctrl(n)
        enddo
        xx = 0.01*xx/a22rt
        bdh = bdh + xx
c
c       Correct B-dot to the current pressure, using a correction of
c       order 2 or higher.
c
        xx = bdotv*dp
        do ipc = 1,ipcv
          n = ipc + 1
          xx = xx - ( dbdtv(ipc)*(dp**n) )/fctrl(n)
        enddo
        xx = xx/a22rt
        bdot = bdot + xx
c
c       The log K for the "Eh" reaction.
c
        xx = -xvfe*dp
        do ipc = 1,ipcv
          n = ipc + 1
          xx = xx - ( dvfe(ipc)*(dp**n) )/fctrl(n)
        enddo
        xx = xx/avcnst
        xlke = xlke + xx
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     This section addresses enthalpy-related parameters.
c
      if (ipch .gt. 0) then
c
c       The A(H) parameter.
c
        xx = 0.
        do ipc = 1,ipch
          xx = xx  + ( dadhh(ipc)*(dp**ipc) )/fctrl(ipc)
        enddo
        adhh = adhh + xx
c
        do ipc = 1,ipch - 1
c
c         The pressure derivatives of the A(H) parameter.
c
          xx = 0.
          do i = ipc + 1,ipch
            n = i - ipc
            xx = xx  + ( dadhh(i)*(dp**n) )/fctrl(n)
          enddo
          dadhh(ipc) = dadhh(ipc) + xx
        enddo
c
c       The B(H) parameter.
c
        xx = 0.
        do ipc = 1,ipch
          xx = xx  + ( dbdhh(ipc)*(dp**ipc) )/fctrl(ipc)
        enddo
        bdhh = bdhh + xx
c
        do ipc = 1,ipch - 1
c
c         The pressure derivatives of the B(H) parameter.
c
          xx = 0.
          do i = ipc + 1,ipch
            n = i - ipc
            xx = xx  + ( dbdhh(i)*(dp**n) )/fctrl(n)
          enddo
          dbdhh(ipc) = dbdhh(ipc) + xx
        enddo
c
c       The B-dot(H) parameter.
c
        xx = 0.
        do ipc = 1,ipch
          xx = xx  + ( dbdth(ipc)*(dp**ipc) )/fctrl(ipc)
        enddo
        bdoth = bdoth + xx
c
        do ipc = 1,ipch - 1
c
c         The pressure derivatives of the B-dot(H) parameter.
c
          xx = 0.
          do i = ipc + 1,ipch
            n = i - ipc
            xx = xx  + ( dbdth(i)*(dp**n) )/fctrl(n)
          enddo
          dbdth(ipc) = dbdth(ipc) + xx
        enddo
c
c       The enthalpy of the "Eh" reaction.
c
        xx = 0.
        do ipc = 1,ipch
          xx = xx  + ( dhfe(ipc)*(dp**ipc) )/fctrl(ipc)
        enddo
        xhfe = xhfe + xx
c
        do ipc = 1,ipch - 1
c
c         The pressure derivatives of the enthalpy of the "Eh" reaction.
c
          xx = 0.
          do i = ipc + 1,ipch
            n = i - ipc
            xx = xx  + ( dhfe(i)*(dp**n) )/fctrl(n)
          enddo
          dhfe(ipc) = dhfe(ipc) + xx
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     This section addresses volume-related parameters.
c
      if (ipcv .gt. 0) then
c
c       The A(V) parameter.
c
        xx = 0.
        do ipc = 1,ipcv
          xx = xx  + ( dadhv(ipc)*(dp**ipc) )/fctrl(ipc)
        enddo
        adhv = adhv + xx
c
        do ipc = 1,ipcv - 1
c
c         The pressure derivatives of the A(V) parameter.
c
          xx = 0.
          do i = ipc + 1,ipcv
            n = i - ipc
            xx = xx  + ( dadhv(i)*(dp**n) )/fctrl(n)
          enddo
          dadhv(ipc) = dadhv(ipc) + xx
        enddo
c
c       The B(V) parameter.
c
        xx = 0.
        do ipc = 1,ipcv
          xx = xx  + ( dbdhv(ipc)*(dp**ipc) )/fctrl(ipc)
        enddo
        bdhv = bdhv + xx
c
        do ipc = 1,ipcv - 1
c
c         The pressure derivatives of the B(V) parameter.
c
          xx = 0.
          do i = ipc + 1,ipcv
            n = i - ipc
            xx = xx  + ( dbdhv(i)*(dp**n) )/fctrl(n)
          enddo
          dbdhv(ipc) = dbdhv(ipc) + xx
        enddo
c
c       The B-dot(V) parameter.
c
        xx = 0.
        do ipc = 1,ipcv
          xx = xx  + ( dbdtv(ipc)*(dp**ipc) )/fctrl(ipc)
        enddo
        bdotv = bdotv + xx
c
        do ipc = 1,ipcv - 1
c
c         The pressure derivatives of the B-dot(V) parameter.
c
          xx = 0.
          do i = ipc + 1,ipcv
            n = i - ipc
            xx = xx  + ( dbdtv(i)*(dp**n) )/fctrl(n)
          enddo
          dbdtv(ipc) = dbdtv(ipc) + xx
        enddo
c
c       The volume of the "Eh" reaction.
c
        xx = 0.
        do ipc = 1,ipcv
          xx = xx  + ( dvfe(ipc)*(dp**ipc) )/fctrl(ipc)
        enddo
        xvfe = xvfe + xx
c
        do ipc = 1,ipcv - 1
c
c         The pressure derivatives of the volume of the "Eh" reaction.
c
          xx = 0.
          do i = ipc + 1,ipcv
            n = i - ipc
            xx = xx  + ( dvfe(i)*(dp**n) )/fctrl(n)
          enddo
          dvfe(ipc) = dvfe(ipc) + xx
        enddo
      endif
c
      end
