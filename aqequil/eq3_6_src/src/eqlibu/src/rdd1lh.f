      subroutine rdd1lh(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uheadr,uheadx,uline1,ulscr)
c
c     This subroutine reads one line containing an expected header
c     or tag string from an EQ3/6 input in menu-style ("D")
c     format. The line generally contains both the header and data.
c
c     EQLIBU/rdd2lh.f is like the present subroutine, but also reads
c     a second line, which is a separator containing a field of dashes
c     between two separators. EQLIBU/rdd1l.f is also very much like
c     the present subroutine. It differs in that there is no expected
c     header.
c
c     This subroutine is called by:
c
c       XCON6/rd6d7.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       nfldmx = the maximum number of fields on line uline1
c       nfldtx = the expected number of fields on line uline1
c       ninpts = the unit number of the stripped input file
c       nlchmx = the maximum number of characters on line uline1
c       nttyo  = the unit number of the screen file
c       uheadx = the expected header (the expected contents of the
c                  first of the two lines)
c       ulscr  = scratch character variable
c
c     Output:
c
c       nfldt  = the number of fields on line uline1
c       qrderr = logical flag, true if an error was encountered
c       ufield = the array of fields on line uline1
c       uheadr = the header (contents of the first of the two lines)
c       uline1 = the line to be read
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nfldmx,nlchmx
c
      integer nttyo
c
      integer nfldt,nfldtx,ninpts
c
      logical qrderr
c
      character*(*) ufield(nfldmx)
      character*(*) uline1,ulscr
      character*(*) uheadr,uheadx
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,j2,j3
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
      call parslj(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
c
c     Compare the number of fields found with that expected.
c
      if (nfldtx .gt. 0) then
        if (nfldt .ne. nfldtx) then
          j3 = ilnobl(uline1(1:70))
          write (nttyo,1010) nfldt,nfldtx,uline1(1:j3)
 1010     format(/' * Warning - (EQLIBU/rdd1lh) Found ',i2,' fields',
     $    /7x,'where ',i2,' were expected on the line beginning with',
     $    /7x,'"',a,'".')
        endif
      endif
c
c     Compare the header with that expected.
c
      uheadr = ufield(1)
      call locase(uheadr)
      call locase(uheadx)
c
      j2 = ilnobl(uheadx)
      j = index(uheadr,uheadx(1:j2))
      if (j .eq. 0) then
        qrderr = .true.
        j3 = ilnobl(uline1(1:70))
        j2 = ilnobl(uheadx(1:70))
        write (nttyo,1020) uheadx(1:j2),uline1(1:j3)
 1020   format(/' * Error - (EQLIBU/rdd1lh) Was expecting to find the',
     $  /7x,'header beginning with',/7x,'"',a,'"',
     $  /7x,'on the line beginning with',/7x,'"',a,'".')
        go to 999
      endif
c
      go to 999
c
  990 qrderr = .true.
c
  999 continue
      end
