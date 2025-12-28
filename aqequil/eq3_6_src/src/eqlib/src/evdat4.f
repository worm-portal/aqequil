      subroutine evdat4(arr,ipc,ipcxmx,k,nmax,narxmx,narxt,ntpr,
     $ ntprmx,prop,tempc)
c
c     This subroutine evaluates a thermodynamic property as a function
c     of temperature, using an interpolating polynomial whose
c     coefficients are stored in a 4D array arr. The second dimension
c     of this array corresponds to a temperature range. The third
c     dimension of this array usually corresponds to an order parameter,
c     as for pressure correction.
c
c     Compare with EQLIB/evdat3.f, in which arr is a 3D array, and
c     EQLIB/evdat3.f, in which this is a 3D array.
c
c     This subroutine is called by:
c
c       EQLIB/evdatr.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       k      = the index for the fourth dimension of the arr array
c       arr    = two dimensional array of polynomial coefficients
c                  describing some thermodynamic function
c       ipcxmx = the third dimension of the arr array
c       ipc    = the index for the third dimension of the arr array
c       tempc  = temperature, C
c       ntpr   = temperature range flag
c       narxmx = first dimension of the arr array, the number
c                  of coefficients per temperature range
c       narxt  = array of numbers of coefficients in temperature
c                  ranges
c       ntprmx = second dimension of the arr array, the number
c                  of temperature ranges.
c       nmax   = third dimension of the arr array (can be ngtmax,
c                  nmtmax, or ngtmax)
c
c     Prinicpal output:
c
c       prop   = the calculated property array
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipcxmx,nmax,narxmx,ntprmx
c
      integer narxt(ntprmx)
c
      integer ipc,k,ntpr
c
      real*8 arr(narxmx,ntprmx,ipcxmx,nmax),prop
      real*8 tempc
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
        prop = arr(n,ntpr,ipc,k) + tempc*prop
      enddo
c
      end
