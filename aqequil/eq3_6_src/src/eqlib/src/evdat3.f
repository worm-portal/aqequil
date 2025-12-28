      subroutine evdat3(arr,k,nmax,narxmx,narxt,ntpr,ntprmx,prop,tempc)
c
c     This subroutine evaluates a thermodynamic property as a function
c     of temperature, using an interpolating polynomial whose
c     coefficients are stored in a 3D array arr. The second dimension
c     of this array corresponds to a temperature range.
c
c     Compare with EQLIB/evdat2.f, in which arr is a 2D array, and
c     EQLIB/evdat4.f, in which this is a 4D array.
c
c
c     This subroutine is called by:
c
c       EQLIB/alters.f
c       EQLIB/evdatr.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       k      = index for the third dimension of array arr
c       arr    = two dimensional array of polynomial coefficients
c                  describing some thermodynamic function
c       tempc  = temperature, C
c       ntpr   = temperature range flag
c       narxmx = first dimension of the arr array, the maximum number
c                  of coefficients per temperature range
c       narxt  = array of numbers of coefficients in temperature
c                  ranges
c       nmax   = third dimension of the arr array (can be ipchmx,
c                  ipcvmx, ngtmax, nmtmax, or ngtmax)
c
c     Prinicpal output:
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
      integer nmax,narxmx,ntprmx
c
      integer narxt(ntprmx)
c
      integer k,ntpr
c
      real*8 arr(narxmx,ntprmx,nmax)
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
        prop = arr(n,ntpr,k) + tempc*prop
      enddo
c
      end
