      subroutine rdd2l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,uline2,ulscr)
c
c     This subroutine reads a "2 line" block containing an unknown
c     header or tag string from an EQ3/6 input in menu-style ("D")
c     format. The first line contains the header and data. The second
c     line is a separator line containing a single field of dashes.
c
c     EQLIBU/rdd2lh.f is very similar to the present subroutine. It
c     differs in that the header is required to match a defined string.
c
c     This subroutine is called by:
c
c       XCON3/rd3d7.f
c       XCON6/rd6d7.f
c       XCON6/rd6d8.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       nfldmx = the maximum number of fields on line uline1 or
c                  line uline2
c       nfldtx = the expected number of fields on a line
c       nlchmx = the maximum number of characters on a line
c       ninpts = the unit number of the stripped input file
c       nttyo  = the unit number of the screen file
c       ulscr  = scratch character variable
c
c     Output:
c
c       nfldt  = the number of fields on line uline1 or line uline2
c       qrderr = logical flag, true if an error was encountered
c       ufield = the array of fields on a line
c       uline1 = the first line to be read
c       uline2 = the second line to be read
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nttyo
c
      integer nfldmx,nfldt,nfldtx,ninpts,nlchmx
c
      logical qrderr
c
      character*(*) ufield(nfldmx)
      character*(*) uline1,uline2,ulscr
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2
c
      integer ilnobl
c
c-----------------------------------------------------------------------
c
      qrderr = .false.
c
c     Read and parse the first line.
c
      read (ninpts,1000,err=990) uline1
 1000 format(a80)
c
      call parslj(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
c
c     Compare the number of fields found with that expected.
c
      if (nfldtx .gt. 0) then
        if (nfldt .ne. nfldtx) then
          j2 = ilnobl(uline1(1:70))
          write (nttyo,1010) nfldt,nfldtx,uline1(1:j2)
 1010     format(/' * Warning - (EQLIBU/rdd2l) Found ',i2,' fields',
     $    /7x,'where ',i2,' were expected on the line beginning with',
     $    /7x,'"',a,'".')
        endif
      endif
c
c     Read the second line.
c
      read (ninpts,1000,err=990) uline2
      if (uline2(2:9) .eq. '--------') go to 999
c
      qrderr = .true.
      j2 = ilnobl(uline2(1:70))
      write (nttyo,1030) uline2(1:j2)
 1030 format(/' * Error - (EQLIBU/rdd2l) Found the line beginning',
     $ ' with',/7x,'"',a,'"',
     $ /7x,'where a dashed separator line was expected.')
      go to 999
c
  990 qrderr = .true.
c
  999 continue
      end
