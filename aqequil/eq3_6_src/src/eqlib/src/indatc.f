      subroutine indatc(arr,nad1,narx_asv,narxt,ntpr_asv,ntprt,ux24)
c
c     This subroutine reads from the unformatted data file (unit number
c     nad1) the 2D array arr, which contains coefficients of
c     interpolating polynomials representing a thermodynamic property
c     as a function of temperature. There are ntprt temperature
c     ranges, and a set of narxt(ntpr) coefficients for the ntpr-th
c     temperature range. When this subroutine is called, another array
c     name is usually substituted for arr.
c
c     Compare with EQLIB/indatd.f, in which arr is a 3D array.
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
c       arr    = 2D array of polynomial coefficients
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer narx_asv,ntpr_asv
c
      integer nad1
c
      integer narxt(ntpr_asv)
c
      integer ntprt
c
      real*8 arr(narx_asv,ntpr_asv)
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
        read (nad1) (arr(n,ntpr), n = 1,narxt(ntpr))
      enddo
c
      end
