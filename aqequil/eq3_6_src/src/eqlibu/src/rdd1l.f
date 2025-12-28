      subroutine rdd1l(nfldmx,nfldt,nfldtx,ninpts,nlchmx,
     $ nttyo,qrderr,ufield,uline1,ulscr)
c
c     This subroutine reads single line containing an unknown header
c     or tag string from an EQ3/6 input in menu-style ("D")
c     format. The line generally contains both the header and data.
c
c     EQLIBU/rdd2l.f is much like the present subroutine, but also reads
c     a second line, which is a separator containing a field of dashes
c     between two separators. EQLIBU/rdd1lh.f is also much like the
c     present subroutine. It differs in that it has an expected header.
c
c     This subroutine is called by:
c
c       XCON3/rd3d7.f
c       XCON6/rd6d7.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       ninpts = the unit number of the stripped input file
c       nfldmx = the maximum number of fields on line uline1
c       nfldtx = the expected number of fields on line uline1
c       nlchmx = the maximum number of characters on line uline1
c       nttyo  = the unit number of the screen file
c       ulscr  = scratch character variable
c
c
c     Output:
c
c       nfldt  = the number of fields on line uline1
c       qrderr = logical flag, true if an error was encountered
c       ufield = the array of fields on line uline1
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
      call parslj(nfldmx,nfldt,nlchmx,ufield,uline1,ulscr)
c
c     Compare the number of fields found with that expected.
c
      if (nfldtx .gt. 0) then
        if (nfldt .ne. nfldtx) then
          j2 = ilnobl(uline1(1:70))
          write (nttyo,1010) nfldt,nfldtx,uline1(1:j2)
 1010     format(/' * Warning - (EQLIBU/rdd1l) Found ',i2,' fields',
     $    /7x,'where ',i2,' were expected on the line beginning with',
     $    /7x,'"',a,'".')
        endif
      endif
c
      go to 999
c
  990 qrderr = .true.
c
  999 continue
      end
