      subroutine rijust(ustr)
c
c     This subroutine right-justifies the non-blank portion of the
c     string ustr.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       ustr   = the input string variable
c
c     Output:
c
c       ustr   = the output string variable
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      character*(*) ustr
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,jj,jbl,j2,nchars
c
      integer ilnobl
c
c-----------------------------------------------------------------------
c
c     Get the length of the string variable.
c
      nchars = len(ustr)
c
c     Get the position of the last non-blank character and the number
c     of blanks on the right-hand-side.
c
      j2 = ilnobl(ustr)
      jbl = nchars - j2
c
      if (jbl .gt. 0) then
        do jj = j2,1,-1
          j = jj + jbl
          ustr(j:j) = ustr(jj:jj)
        enddo
        do j = 1,jbl
          ustr(j:j) = ' '
        enddo
      endif
c
      end
