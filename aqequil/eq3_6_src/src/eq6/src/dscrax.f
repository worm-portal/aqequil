      subroutine dscrax(nf1,nf2,nllnmx,ulinex)
c
c     This subroutine descrambles a file of tables whose lines are
c     interspersed, but which are marked 'a', 'b', 'c', etc., in
c     column one. The contents of the scrambled file are copied to
c     the descrambled file as the descrambling takes place. The
c     descrambled file must already be open. The maximum line length
c     of this file is nllnmx characters.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       nf1    = unit number of the scrambled file
c       nf2    = unit number of the descrambled file
c       nllnmx = the maximum character length of a line on the
c                descrambled files; the maximum character length
c                for the scrambled file is nllnmx + 1
c       ulinex = variable holding a line
c
c     Output:
c
c       None
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nf1,nf2
c
      integer nllnmx
c
      character(len=nllnmx) ulinex
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,i1,i2,j1
c
      character(len=1) uc,ux
c
c-----------------------------------------------------------------------
c
      i1 = ichar('a')
      i2 = ichar('Z')
      if (i2 .lt. i1) then
        i1 = ichar('A')
        i2 = ichar('z')
      endif
c
      rewind nf2
      do i = i1,i2
        uc = char(i)
        rewind nf1
  100   read (nf1,'(a1,a)',end=110) ux,ulinex
        if (ux .eq. uc) then
          j1 = len_trim(ulinex)
          write (nf2,'(a)') ulinex(1:j1)
        endif
        go to 100
c
  110   continue
      enddo
c
      end
