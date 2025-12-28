      subroutine indatd(arr,ipc,ipcx_asv,nad1,narxt,narx_asv,ntprt,
     $ ntpr_asv,ux24)
c
c     This subroutine reads from the unformatted data file (unit number
c     nad1) the 3D array arr, which contains coefficients of
c     interpolating polynomials representing a thermodynamic property
c     as a function of temperature. There are ntprt temperature
c     ranges, and a set of narxt(ntpr) coefficients for the ntpr-th
c     temperature range for the ipc-th entity. When this subroutine is
c     called, another array name is usually substituted for arr and
c     a corresponding dimensioning variable is substituted for ipcx_asv.
c
c     Compare with EQLIB/indatc.f, in which arr is a 2D array.
c
c     This subroutine is called by:
c
c       EQLIB/indata.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nad1   = unit number of the data file
c
c     Principal output:
c
c       arr    = 3D array of polynomial coefficients
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipcx_asv,narx_asv,ntpr_asv
c
      integer nad1
c
      integer narxt(ntpr_asv)
c
      integer ipc,ntprt
c
      real*8 arr(narx_asv,ntpr_asv,ipcx_asv)
c
      character*24 ux24
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,ntpr
c
c-----------------------------------------------------------------------
c
c     The content of ux24 may be useful in debugging, but has no
c     other usage.
c
      read (nad1) ux24
      do ntpr = 1,ntprt
        read (nad1) (arr(n,ntpr,ipc), n = 1,narxt(ntpr))
      enddo
c
      end
