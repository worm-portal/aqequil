      subroutine intchr(ivar,noutpt,nttyo,ustr)
c
c     This subroutine converts writes the integer variable ivar into the
c     character variable ustr, employing left justification and blank
c     fill. If right justification is desired, one should just write
c     the integer into the string.
c
c     The length of the string variable is unknown. A buffer variable
c     which is employed has a character length of ichpar. This limits
c     the size of the substring available for the input integer.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       ivar   = the input integer
c       noutpt = unit number of the output file
c       nttyo  = unit number of the screen file
c
c     Output:
c
c       ustr   = the string variable in which the integer is to
c                  be written
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt,nttyo,ivar
c
      character*(*) ustr
c
c-----------------------------------------------------------------------
c
c     Local parameter declarations.
c
      integer ichpar
c
      parameter (ichpar = 16)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,ilen
c
      integer ilnobl
c
      character*(ichpar) ux
c
c-----------------------------------------------------------------------
c
c     Write the integer into a character variable buffer.
c
      write (ux,1000,err = 100) ivar
c
c     The format below should specify an integer field which matches
c     the character length of the buffer variable (ichpar).
c
 1000 format(i16)
      go to 110
c
  100 write (noutpt,1010) ivar
      write (nttyo ,1010) ivar
 1010 format(/" * Error - (EQLIBU/intchr) Can't write the integer",
     $ /7x,'"',i16,'" into a character string due to formatted',
     $ /7x,'write error.')
      stop
c
c     Left justify.
c
  110 call lejust(ux)
c
c     Get the length of the integer string.
c
      j2 = ilnobl(ux)
c
c     Get the length of the string variable.
c
      ilen = len(ustr)
c
      if (j2 .gt. ilen) then
        write (noutpt,1020) ivar,j2,ilen
        write (nttyo ,1020) ivar,j2,ilen
 1020   format(/" * Error - (EQLIBU/intchr) Can't write the integer",
     $ /7x,'"',i16,'" into a string because the integer requires',i3,
     $ /7x,' characters, but the string variable has only ',i3,'.')
        stop
      endif
c
c     Load the string.
c
      ustr = ux(1:j2)
c
      end
