      subroutine gntpr(ntpr,ntprmx,ntprt,tempc,tempcu)
c
c     This subroutine finds the value of ntpr, the temperature range
c     flag.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c       EQ6/tpadv.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ntprt  = the number of temperature ranges
c       tempc  = temperature, C
c       tempcu = array of temperatures, C, defining the upper boundaries
c                  of the temperature ranges
c
c     Principal output:
c
c       ntpr   = temperature range flag
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ntprmx
c
      integer ntpr,ntprt
c
      real*8 tempcu(ntprmx)
c
      real*8 tempc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c     None
c
c-----------------------------------------------------------------------
c
c     Determine the temperature range flag (ntpr).
c
      do ntpr = 1,ntprt
        if (tempc .le. tempcu(ntpr)) go to 100
      enddo
      ntpr = ntprt
  100 continue
c
      end
