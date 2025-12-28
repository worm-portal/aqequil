      real*8 function tlg(x)
c
c     This subroutine is the function log10, except that the log10(0.)
c     is returned with a value of -99999.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       x      = argument
c
c     Output:
c
c       tlg    = log10 x (except that tlg(0.) = -99999.)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      real*8 x
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c       None
c
c-----------------------------------------------------------------------
c
      if (x .eq. 0.) then
        tlg = -99999.
      else
        tlg = log10(x)
      endif
c
      end
