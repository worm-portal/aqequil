      subroutine dadjdx(delxi,dlxmin,iodb,nodbmx,noutpt,qadjdx)
c
c     This subroutine determines whether or not delxi can be reduced
c     to satisfy some criterion, such as the pH not exceeding the
c     requested maximum value.
c
c     This subroutine is called by:
c
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer nodbmx
c
      integer noutpt
c
      integer iodb(nodbmx)
c
      logical qadjdx
c
      real*8 delxi,dlxmin
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c     None
c
c-----------------------------------------------------------------------
c
      qadjdx = .false.
      if (delxi .le. dlxmin) then
c
c       The step size is already at the minimum value.
c       Do not cut it.
c
        if (iodb(1) .gt. 0) write (noutpt,1100)
 1100    format(3x,'The step size will not be cut because it is',
     $   ' already at the',/5x,'minimum value.',/)
      else
c
c       Set up to go back and cut the step size to satisfy the
c       accuracy criterion for calculating the point at which the
c       event of concern occurs.
c
        qadjdx = .true.
      endif
c
      end
