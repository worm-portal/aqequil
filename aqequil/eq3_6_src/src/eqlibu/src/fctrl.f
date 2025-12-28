      real*8 function fctrl(i)
c
c     This subroutine returns the factorial function. The factorials
c     from 0! to 10! are programmed into a data statement.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       i     = argument
c
c     Output:
c
c       fctrl = factorial
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer i
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ifac(0:10),ifactr,j
c
c-----------------------------------------------------------------------
c
c     Set low factorials.
c
      data (ifac(j), j = 0,10) /1,1,2,6,24,120,720,5040,40320,
     $ 362880,3628800/
c
c-----------------------------------------------------------------------
c
      if ( i.ge.0 .and. i.le.10) then
        fctrl = ifac(i)
      elseif (i.gt.10) then
        ifactr = ifac(10)
        do j = 10 + 1,i
          ifactr = ifactr*j
        enddo
        fctrl = ifactr
      endif
c
      end
