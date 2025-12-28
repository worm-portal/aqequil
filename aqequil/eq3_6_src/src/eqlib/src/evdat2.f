      subroutine evdat2(arr,narxmx,narxt,ntpr,ntprmx,prop,tempc)
c
c     This subroutine evaluates a thermodynamic property as a function
c     of temperature, using an interpolating polynomial whose
c     coefficients are stored in a 2D array arr. The second dimension
c     of this array corresponds to a temperature range.
c
c     Compare with EQLIB/evdat3.f, in which arr is a 3D array, and
c     EQLIB/evdat4.f, in which this is a 4D array.
c
c     This subroutine is called by:
c
c       EQLIB/alters.f
c       EQLIB/evdata.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       arr    = two dimensional array of polynomial coefficients
c                  describing some thermodynamic function
c       tempc  = temperature, C
c       ntpr   = temperature range flag
c       narxmx = first dimension of the arr array, the maximum number
c                  of coefficients per temperature range
c       narxt  = array of numbers of coefficients in temperature
c                  ranges
c       ntprmx = second dimension of the arr array, the number
c                  of temperature ranges.
c
c     Principal output:
c
c       prop   = the calculated property
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer narxmx,ntprmx
c
      integer narxt(ntprmx)
c
      integer ntpr
c
      real*8 arr(narxmx,ntprmx)
      real*8 prop,tempc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,nn,nt
c
c-----------------------------------------------------------------------
c
c     Evaluate the polynomial.
c
      prop = 0.
      nt = narxt(ntpr)
      do nn = 1,nt
        n = nt + 1 - nn
        prop = arr(n,ntpr) + tempc*prop
      enddo
c
      end
