      subroutine stripl(nin,nout)
c
c     This subroutine copies the file whose unit number is "nin" to
c     that whose unit number is "nout". Lines beginning with an asterix
c     are not copied. Lines exceeding a length of 80 characters are
c     truncated.
c
c     This subroutine is called by:
c
c       EQPT/ofiles.f
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       nin    = unit number of the file to be stripped
c       nout   = unit number of the stripped file
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
      integer nin,nout
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2
c
      integer ilnobl
c
      character*80 uibuf
      character*1 ustar,ux
c
c-----------------------------------------------------------------------
c
      data ustar  /'*'/
c
c-----------------------------------------------------------------------
c
   10 read (nin,1000,end=999) uibuf
 1000 format(a80)
c
c     Skip lines with an asterisk in column 1.
c
      ux = uibuf(1:1)
      if (ux .eq. ustar) go to 10
c
      j2 = ilnobl(uibuf)
      if (j2 .le. 0) j2 = 1
      write (nout,1010) uibuf(1:j2)
 1010 format(a)
      go to 10
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
