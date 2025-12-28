      subroutine openin(noutpt,nttyo,ufiln,uform,ilu)
c
c     This subroutine opens an input type file. The file must already
c     exist. An unused logical unit number is obtained.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       noutpt = the unit number of the output file
c       nttyo  = the unit number of the screen file
c       ufiln  = the name of the file to open
c       uform  = the file format, 'formatted' or 'unformatted'
c
c     Output:
c
c       ilu     = the logical unit number of opened file
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
c     Note- the character length for the variables ufiln and uform
c     can not be exactly specified.
c
      integer ilu,noutpt,nttyo
c
      character*(*) ufiln,uform
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,nerr
c
      integer ilnobl
c
      logical qex
c
      character*8 uformo
c
c-----------------------------------------------------------------------
c
      j2 = ilnobl(ufiln)
c
c     Check to make sure the file exists.
c
      inquire(file=ufiln,exist=qex,formatted=uformo)
      if (.not.qex) then
        if (noutpt .gt. 0) write (noutpt,1000) ufiln(1:j2)
        write (nttyo,1000) ufiln(1:j2)
 1000   format(/' * Error - (EQLIBU/openin) The file "',a,'"',
     $  /7x,'does not exist.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the file format.
c
      if (uform.eq.'formatted' .and. uformo.eq.'no') then
        if (noutpt .gt. 0) write (noutpt,1010) ufiln(1:j2)
        write (nttyo,1010) ufiln(1:j2)
 1010   format(/' * Error - (EQLIBU/openin) The file "',a,'"',
     $  /7x,'should be formatted, but it is not.')
        stop
      endif
c
      if (uform.eq.'unformatted' .and .uformo.eq.'yes') then
c       if (noutpt .gt. 0) write (noutpt,1020) ufiln(1:j2)
c       write (nttyo,1020) ufiln(1:j2)
c1020   format(/' * Error - (EQLIBU/openin) The file "',a,'"',
c    $  /7x,'should be unformatted, but it is not.')
c     stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the next available logical unit number.
c
      call getlu(ilu,nerr)
      if (nerr .ne. 0) then
        if (noutpt .gt. 0) write (noutpt,1050) ufiln(1:j2)
        write (nttyo,1050) ufiln(1:j2)
 1050   format(/' * Error - (EQLIBU/openin) No logical unit number',
     $  /7x,'is available for the file "',a,'".')
      stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      open(ilu,file=ufiln,form=uform,status='old',err=10)
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
   10 continue
      if (noutpt .gt. 0) write (noutpt,1060) ufiln(1:j2)
      write (nttyo,1060) ufiln(1:j2)
 1060 format(/" * Error - (EQLIBU/openin) Can't open the file",
     $ /7x,'"',a,'".')
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
