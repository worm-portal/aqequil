      subroutine dscram(nf1,nf2)
c
c     This subroutine unscrambles a file of tables whose lines are
c     interspersed, but which are marked 'a', 'b', 'c', etc., in
c     column one. The contents of the scrambled file are copied to
c     the unscrambled file as the unscrambling takes place. The
c     unscrambled file must already be open. The record length of
c     this file will is set at 128 characters.
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
c       nf2   = unit number of the unscrambled file
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
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,i1,i2,j
c
      character*128 uline
      character*1 uc,ux
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
        do j = 1,10000
          read (nf1,1000,end=100) ux,uline
 1000     format(a1,a128)
          if (ux .eq. uc) write (nf2,1005) uline
 1005     format(a128)
        enddo
  100   continue
      enddo
c
      end
