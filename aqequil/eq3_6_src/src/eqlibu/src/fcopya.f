      subroutine fcopya(nf1,nf2)
c
c     This subroutine appends the contents of the file whose unit number
c     is nf1 to the file whose unit number is nf2. The line length
c     is assumed to be 128 characters. The second file  must already
c     be open and correctly positioned.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       nf1    = unit number of the first file
c       nf2    = unit number of the second file
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
      integer j
c
      character*128 uline
c
c----------------------------------------------------------------------
c
      rewind nf1
      do j = 1,10000
        read (nf1,1000,end=999) uline
 1000   format(a128)
        write (nf2,1000) uline
      enddo
c
  999 continue
      end
