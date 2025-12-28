      subroutine cejust(ustr)
c
c     This subroutine centers the non-blank portion of the string ustr.
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
      integer j,jj,jbl,jbll,jblr,jshl,jshr,j1,j2,nchars
c
      integer ifnobl,ilnobl
c
c-----------------------------------------------------------------------
c
c     Get the length of the string variable.
c
      nchars = len(ustr)
c
c     Get the position of the first non-blank character and the number
c     of blanks on the left-hand-side.
c
      j1 = ifnobl(ustr)
      jbll = j1 - 1
c
c     Get the position of the last non-blank character and the number
c     of blanks on the right-hand-side.
c
      j2 = ilnobl(ustr)
      jblr = nchars - j2
c
      jbl = (jbll + jblr)/2
      jshl = jbll - jbl
c
      if (jshl .gt. 0) then
c
c       Shift left.
c
        do jj = j1,nchars
          j = jj - jshl
          ustr(j:j) = ustr(jj:jj)
        enddo
        do j = nchars - jshl,nchars
          ustr(j:j) = ' '
        enddo
      endif
c
      jshr = - jshl
      if (jshr .gt. 1) then
        do jj = j2,1,-1
          j = jj + jshr
          ustr(j:j) = ustr(jj:jj)
        enddo
        do j = 1,jshr
          ustr(j:j) = ' '
        enddo
      endif
c
      end
