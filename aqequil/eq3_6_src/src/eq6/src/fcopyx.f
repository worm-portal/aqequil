      subroutine fcopyx(nf1,nf2,nllnmx,ulinex)
c
c     This subroutine appends the contents of the file whose unit
c     number is nf1 to the file whose unit number is nf2. The line
c     length is assumed to be nllnmx characters. The second file
c     must already be open and correctly positioned.
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
c       nllnmx = the maximum character length of a line on the
c                scrambled and descrambled files
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
      integer j1
c
      character(len=nllnmx) ulinex
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j
c
c----------------------------------------------------------------------
c
      rewind nf1
  100 read (nf1,'(a)',end=110) ulinex
      j1 = len_trim(ulinex)
      write (nf2,'(a)') ulinex(1:j1)
      go to 100
c
  110 continue
c
      end
