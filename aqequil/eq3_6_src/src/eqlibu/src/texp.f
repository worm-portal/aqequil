      real*8 function texp(x)
c
c     This subroutine is the function 10**x. The argument is tested in
c     order to avoid overflow.
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
c       texp   = 10**x (maximum returned value is 1.e+125)
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
c     Test for overflow. The recommended truncation limit on the
c     argument x is 125. On most systems, the real*8 exponent
c     range is at least 308, which provides sufficient cushion.
c     In considering how much cushion is enough, remember that
c     exponentiated numbers may be subsequently be squared (in which
c     case 1.e+125 gives rise to 1.e+250) or otherwise raised to some
c     positive power. If your machine has a smaller but still acceptable
c     exponent range (no less than 100), you should adjust the the limit
c     used below accordingly. For example, if your machine's exponent
c     range is 100, try a truncation limit of about 1.e+40. Your
c     machine's exponent range is actually computed by EQLIBU/flpars.f.
c
      if (x .ge. 125.) then
c
c       Avoid overflow.
c
com     BEGIN_MACHINE_DEPENDENT_CODE
com
com       Note: Some compilers may not allow use of the construction
com       "1.e+125" in the following setting. They may infer from the
com       "e" part that this is a REAL*4 number, and conclude that the
com       exponent is out of range. The "d" construction is appropriate
com       for any 32-bit machine. It shouldn't cause a problem on a
com       64-bit machine, such as a CRAY. If it should, however, change
com       to the "e" construction.
com
          texp = 1.d+125
cxx       texp = 1.e+125
com
com     END_MACHINE_DEPENDENT_CODE
c
      else
        texp = 10.**x
      endif
c
      end
