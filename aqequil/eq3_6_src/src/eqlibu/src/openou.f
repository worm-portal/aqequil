      subroutine openou(noutpt,nttyo,ufiln,uform,nrecl,ilu)
c
c     This subroutine opens an output type file. If a file of the
c     same name already exists, it is first destroyed. An unused
c     logical unit number is obtained.
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
c       nrecl  = the record length (number of characters per line)
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
      integer noutpt,nttyo
c
      integer ilu,nrecl
c
c     Note- the character length for the variables ufiln and uform
c     can not be exactly specified.
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
      character*8 ustat
c
c-----------------------------------------------------------------------
c
      j2 = ilnobl(ufiln)
c
c     See if a file of the same name already exists. If so,
c     destroy it. This makes the logical unit number available.
c     If a file of the same name does not exist, get the next
c     available logical unit number.
c
      inquire(file=ufiln,exist=qex)
      if (qex) then
        ustat = 'old'
        call getlu(ilu,nerr)
        if (nerr .ne. 0) then
          if (noutpt .gt. 0) write (noutpt,1000) ustat,ufiln(1:j2)
          write (nttyo,1000) ustat,ufiln(1:j2)
 1000     format(/' * Error - (EQLIBU/openou) No logical unit number',
     $    /7x,'is available to open the ',a3,' file "',a,'".')
          stop
        endif
        open(ilu,file=ufiln,status=ustat,err=10)
        close(ilu,status='delete',err=15)
      else
        ustat = 'new'
        call getlu(ilu,nerr)
        if (nerr .ne. 0) then
          if (noutpt .gt. 0) write (noutpt,1000) ustat,ufiln(1:j2)
          write (nttyo,1000) ustat,ufiln(1:j2)
          stop
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Open the new file.
c
      if (nrecl .gt. 0) then
c
c       Use the specified record length.
c
        open(ilu,file=ufiln,form=uform,status='new',recl=nrecl,err=10)
      else
c
c       The record length is not specified. Open the file at the
c       default record length.
c
        open(ilu,file=ufiln,form=uform,status='new',err=10)
      endif
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
   10 if (noutpt .gt. 0) write (noutpt,1010) ustat,ufiln(1:j2)
      write (nttyo ,1010) ustat,ufiln(1:j2)
 1010 format(/" * Error - (EQLIBU/openou) Can't open the ",a3,' copy',
     $ /7x,'of the file "',a,'".')
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
   15 if (noutpt .gt. 0) write (noutpt,1020) ufiln(1:j2)
      write (nttyo,1020) ufiln(1:j2)
 1020 format(/" * Error - (EQLIBU/openou) Can't delete the old copy",
     $ /7x,'of the file "',a,'".')
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
