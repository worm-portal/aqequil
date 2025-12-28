      subroutine pcorrx(avcnst,dhfs,dvfs,ipch,ipchmx,ipcv,ipcvmx,
     $ nbasp,nbt,nbtmax,ndrsr,nst,nstmax,presg,press,xhfs,
     $ xlks,xvfs)
c
c     This subroutine makes pressure corrections for equilibrium
c     constants and related thermodynamic functions. It normally
c     corrects for pressures off the standard T-P grid (e.g., from
c     "presg" to "press"). However, it can be used to correct from
c     any pressure to another by substituting the value of the former
c     for "presg". This is done in EQ6 to correct for changing pressure
c     along an isothermal reaction path. In order to do this, the
c     pressure derivatives of enthalpy and volume functions (dhfs,dvfs)
c     are also corrected to the pressure of interest. The derivatives of
c     highest order are necessarily treated as constants, so there
c     is no correction in this case.
c
c     EQLIB/pcorrm.f performs the same function for some miscellaneous
c     thermodynamic functions, such as Debye-Huckel
c     parameters.
c
c     This subroutine is called by:
c
c       EQ3NR/arrset.f
c       EQ3NR/eq3nr.f
c       EQ6/absswb.f
c       EQ6/eq6.f
c       EQ6/tpadv.f
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
      integer ipchmx,ipcvmx,nbtmax,nstmax
c
      integer ipch,ipcv,nbt,nst
c
      integer nbasp(nbtmax),ndrsr(2,nstmax)
c
      real*8 dhfs(ipchmx,nstmax),dvfs(ipcvmx,nstmax),xhfs(nstmax),
     $ xlks(nstmax),xvfs(nstmax)
c
      real*8 avcnst,presg,press
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ipc,n,nb,ns,nt
c
      real*8 dp,xx
c
      real*8 fctrl
c
c-----------------------------------------------------------------------
c
c     Calculate the pressure difference.
c
      dp = press - presg
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ipcv .eq. 0) then
c
c       Correct the log K values to the current pressure, using
c       a first order (constant volume of reaction) correction.
c       This is expected to be the most common order for pressure
c       corrections, hence this special block has been written
c       for the sake of efficiency.
c
        do ns = 1,nst
          if (xlks(ns).gt.-9999999. .and. xlks(ns).lt.9999999.) then
            xx = -xvfs(ns)*dp/avcnst
            xlks(ns) = xlks(ns) + xx
          endif
        enddo
      endif
c
      if (ipcv .gt. 0) then
c
c       Correct the log K values to the current pressure, using
c       a correction of order 2 or higher. The coding in this block
c       could be used to handle the order 1 case. However, that
c       would not be quite as efficient as the above block.
c
        do ns = 1,nst
          if (xlks(ns).gt.-9999999. .and. xlks(ns).lt.9999999.) then
            xx = -xvfs(ns)*dp
            do ipc = 1,ipcv
              n = ipc + 1
              xx = xx - ( dvfs(ipc,ns)*(dp**n) )/fctrl(n)
            enddo
            xx = xx/avcnst
            xlks(ns) = xlks(ns) + xx
          endif
        enddo
      endif
c
      if (ipcv .ge. 0) then
c
c       Uncorrect the log K values for strict basis species,
c       as these are fixed at zero and the volume function array
c       contains the standard partial molar volume, not the
c       standard partial molar volume of reaction.
c
        do nb = 1,nbt
          ns = nbasp(nb)
          nt = ndrsr(2,ns) - ndrsr(1,ns) + 1
          if (nt .lt. 2) then
            xlks(ns) = 0.
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ipch .gt. 0) then
c
c       Correct the enthalpy function values to the current pressure.
c
        do ns = 1,nst
c
          xx = 0.
          do ipc = 1,ipch
            xx = xx  + ( dhfs(ipc,ns)*(dp**ipc) )/fctrl(ipc)
          enddo
          xhfs(ns) = xhfs(ns) + xx
c
          do ipc = 1,ipch - 1
c
c           Correct the pressure derivatives of the enthalpy functions
c           to the current pressure.
c
            xx = 0.
            do i = ipc + 1,ipch
              n = i - ipc
              xx = xx  + ( dhfs(i,ns)*(dp**n) )/fctrl(n)
            enddo
            dhfs(ipc,ns) = dhfs(ipc,ns) + xx
          enddo
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ipcv .gt. 0) then
c
c       Correct the volume function values to the current pressure.
c
        do ns = 1,nst
c
          xx = 0.
          do ipc = 1,ipcv
            xx = xx  + ( dvfs(ipc,ns)*(dp**ipc) )/fctrl(ipc)
          enddo
          xvfs(ns) = xvfs(ns) + xx
c
          do ipc = 1,ipcv - 1
c
c           Correct the pressure derivatives of the volume functions
c           to the current pressure.
c
            xx = 0.
            do i = ipc + 1,ipcv
              n = i - ipc
              xx = xx  + ( dvfs(i,ns)*(dp**n) )/fctrl(n)
            enddo
            dvfs(ipc,ns) = dvfs(ipc,ns) + xx
          enddo
        enddo
      endif
c
      end
