      subroutine getlu(nlu,nerr)
c
c     This subroutine finds a currently unused unit number.
c
c     This subroutine is called by:
c
c       EQLIBU/openin.f
c       EQLIBU/openou.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       None
c
c     Output:
c
c       nlu    = first currently unused unit number
c       nerr   = error flag:
c                  = 0   Okay
c                  = 1   Error
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nlu,nerr
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer iumax,iumin
c
      logical qopen
c
c-----------------------------------------------------------------------
c
      data iumax,iumin /40,0/
c
c-----------------------------------------------------------------------
c
      nerr = 0
c
c     Loop through all valid file numbers, beginning with the largest.
c
      do nlu = iumax,iumin,-1
        inquire(unit=nlu,opened=qopen)
        if (.not.qopen) go to 999
      enddo
c
      nerr = 1
c
  999 continue
      end
