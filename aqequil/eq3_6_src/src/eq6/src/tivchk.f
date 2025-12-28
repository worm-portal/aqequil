      subroutine tivchk(deltim,delxi,qtvchk,time1,time0,timemx,
     $ tiplol,tiplot,tiprnl,tiprnt,tolxst)
c
c     This subroutine checks to make sure that the calculated
c     time does not exceed any specified limits such as the maximum
c     time. This routine should be called only if delxi is less than
c     or equal to the minimum step size, dlxmin.
c
c     This subroutine is called by:
c
c       EQ6/path.f
c       EQ6/eqshel.f
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
      logical qtvchk
c
      real(8) deltim,delxi,time1,time0,timemx,tiplol,tiplot,tiprnl,
     $ tiprnt,tolxst
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c     None
c
c-----------------------------------------------------------------------
c
      qtvchk = .false.
c
c     Check the next print point in time.
c
      if (((time1 - tiprnt)/tiprnt) .gt. tolxst) then
        time1 = tiprnt
        deltim = tiprnt - time0
        qtvchk = .true.
      endif
c
c     Check the next print point in log time.
c
      if (((time1 - tiprnl)/tiprnl) .gt. tolxst) then
        time1 = tiprnl
        deltim = tiprnl - time0
        qtvchk = .true.
      endif
c
c     Check the next plot point in time.
c
      if (((time1 - tiplot)/tiplot) .gt. tolxst) then
        time1 = tiplot
        deltim = tiplot - time0
        qtvchk = .true.
      endif
c
c     Check the next plot point in log time.
c
      if (((time1 - tiplol)/tiplol) .gt. tolxst) then
        time1 = tiplol
        deltim = tiplol - time0
        qtvchk = .true.
      endif
c
c     Check the maximum time.
c
      if (((time1 - timemx)/timemx) .gt. tolxst) then
        time1 = timemx
        deltim = timemx - time0
        qtvchk = .true.
      endif
c
      end
