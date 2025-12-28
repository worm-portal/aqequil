      subroutine evdatr(adhfs,advfs,axhfs,axlks,axvfs,dhfs,dvfs,
     $ ipch,ipchmx,ipcv,ipcvmx,narxmx,narxt,nst,nstmax,ntpr,ntprmx,
     $ tempc,xhfs,xlks,xvfs)
c
c     This subroutine evaluates equilibrium constants and related
c     thermodynamic functions as functions of temperature and
c     standard grid pressure. To make pressure corrections off the
c     standard T-P grid for these data, use EQLIB/pcorrx.f.
c
c     This subroutine is called by:
c
c       EQLIB/absswa.f
c       EQLIB/evdata.f
c       EQ3NR/eq3nr.f
c       EQ6/absswb.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
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
      integer ipchmx,ipcvmx,narxmx,nstmax,ntprmx
c
      integer narxt(ntprmx)
c
      integer ipch,ipcv,nst,ntpr
c
      real*8 adhfs(narxmx,ntprmx,ipchmx,nstmax),
     $ advfs(narxmx,ntprmx,ipcvmx,nstmax),axhfs(narxmx,ntprmx,nstmax),
     $ axlks(narxmx,ntprmx,nstmax),axvfs(narxmx,ntprmx,nstmax),
     $ dhfs(ipchmx,nstmax),dvfs(ipcvmx,nstmax),xhfs(nstmax),
     $ xlks(nstmax),xvfs(nstmax)
c
      real*8 tempc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ipc,ns
c
      real*8 prop
c
c-----------------------------------------------------------------------
c
c     Compute the log K values for all reactions.
c
      do ns = 1,nst
c
c       Compute the log K values on the standard grid.
c
        if (axlks(1,ntpr,ns) .lt. 9999999.) then
c
c         Calling sequence substitutions:
c           axlks for arr
c           ns for k
c           nstmax for nmax
c
          call evdat3(axlks,ns,nstmax,narxmx,narxt,ntpr,ntprmx,prop,
     $    tempc)
          xlks(ns) = prop
        else
          xlks(ns) = 9999999.
        endif
c
        if (ipch .ge. 0) then
c
c         Compute the enthalpy function values on the standard grid.
c
          if (axhfs(1,ntpr,ns) .lt. 9999999.) then
c
c           Calling sequence substitutions:
c             axhfs for arr
c             ns for k
c             nstmax for nmax
c
            call evdat3(axhfs,ns,nstmax,narxmx,narxt,ntpr,ntprmx,
     $      prop,tempc)
            xhfs(ns) = prop
          else
            xhfs(ns) = 9999999.
          endif
c
          do ipc = 1,ipch
c
c           Compute the pressure derivatives of the enthalpy function
c           values on the standard grid.
c
c           Calling sequence substitutions:
c             adhfs for arr
c             ipchmx for ipcxmx
c             ns for k
c             nstmax for nmax
c
            call evdat4(adhfs,ipc,ipchmx,ns,nstmax,narxmx,narxt,ntpr,
     $      ntprmx,prop,tempc)
            dhfs(ipc,ns) = prop
          enddo
        endif
c
        if (ipcv .ge. 0) then
c
c         Compute the volume function values on the standard grid.
c
          if (axvfs(1,ntpr,ns) .lt. 9999999.) then
c
c           Calling sequence substitutions:
c             axvfs for arr
c             ns for k
c             nstmax for nmax
c
            call evdat3(axvfs,ns,nstmax,narxmx,narxt,ntpr,ntprmx,
     $      prop,tempc)
            xvfs(ns) = prop
          else
            xvfs(ns) = 9999999.
          endif
c
          do ipc = 1,ipcv
c
c           Compute the pressure derivatives of the volume function
c           values on the standard grid.
c
c           Calling sequence substitutions:
c             advfs for arr
c             ipcvmx for ipcxmx
c             ns for k
c             nstmax for nmax
c
            call evdat4(advfs,ipc,ipcvmx,ns,nstmax,narxmx,narxt,ntpr,
     $      ntprmx,prop,tempc)
            dvfs(ipc,ns) = prop
          enddo
        endif
c
      enddo
c
      end
