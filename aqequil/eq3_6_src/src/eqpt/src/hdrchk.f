      subroutine hdrchk(ndat0s,noutpt,nttyo)
c
c     This suboutine checks the first line of the DATA0 file to ensure
c     that the mandatory header ("data0" beginning in column 1) is
c     present.
c
c     This suboutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ndat0s = unit number of the stripped DATA0 file
c
c     Principal output:
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
      integer ndat0s,noutpt,nttyo
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2
c
      integer ilnobl
c
      character(len=80) uline
c
c-----------------------------------------------------------------------
c
c     Check the DATA0 file header. This insures that the file that
c     is supposed to be a DATA0 file really is one.
c
      read (ndat0s,1000,end=990,err=995) uline
 1000 format(a)
c
      if (uline(1:5).ne.'data0' .and. uline(1:5).ne.'Data0' .and.
     $  uline(1:5).ne.'DATA0') then
        j2 = ilnobl(uline)
        j2 = min(j2,50)
        write (noutpt,1010) uline(1:j2)
        write (nttyo,1010) uline(1:j2)
 1010   format(/' * Error - (EQPT/hdrchk) The DATA0 must have "data0"',
     $  ' beginning in',/7x,'column 1 of the first line. The first',
     $  ' line of this file begins',/7x,'instead with: "',a,'".')
        stop
      endif
c
      rewind(ndat0s)
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write a message for any read error.
c
  990 write (noutpt,2000)
      write (nttyo,2000)
 2000 format(/' * Error - (EQPT/hdrchk) The DATA0 file is an empty',
     $ ' file.')
      stop
c
  995 write (noutpt,2010)
      write (nttyo,2010)
 2010 format(/' * Error - (EQPT/hdrchk) Encountered a read format',
     $ ' error while',/7x,'checking for the mandatory header on the',
     $ ' first line of the',/7x,'DATA0 file.')
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
