      subroutine realch(noutpt,nttyo,ustr,var)
c
c     This subroutine converts writes the real*8 variable var into the
c     string ustr, employing left justification and blank fill. If
c     right justification is desired, one should just write the
c     real*8 variable into the string.
c
c     The length of the string variable is unknown. A buffer variable
c     which is employed has a character length of ichpar. This limits
c     the size of the substring available for the input real*8 number.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       noutpt = unit number of the output file
c       nttyo  = unit number of the screen file
c       var   = the input real*8 number
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
      integer noutpt,nttyo
c
      character*(*) ustr
c
      real*8 var
c
c-----------------------------------------------------------------------
c
c     Local parameter declarations.
c
      integer ichpar
c
      parameter (ichpar = 24)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ilen,j2
c
      integer ilnobl
c
      character*(ichpar) ux
c
c-----------------------------------------------------------------------
c
c     Write the real*8 number into a character variable buffer.
c
      write (ux,1000,err = 100) var
c
c     The format below should specify a real*8 field which matches
c     the character length of the buffer variable (ichpar).
c
 1000 format(g24.4)
      go to 110
c
  100 write (noutpt,1010) var
      write (nttyo ,1010) var
 1010 format(/" * Error - (EQLIBU/realch) Can't write the real*8",
     $ /7x,'number ",g24.4," into a character string due to formatted',
     $ /7x,'write error.')
      stop
c
c     Left justify.
c
  110 call lejust(ux)
c
c     Get the length of the real*8 string.
c
      j2 = ilnobl(ux)
c
c     Get the length of the string variable.
c
      ilen = len(ustr)
c
      if (j2 .gt. ilen) then
        write (noutpt,1020) var,j2,ilen
        write (nttyo ,1020) var,j2,ilen
 1020   format(/" * Error - (EQLIBU/realch) Can't write the real*8",
     $ /7x,'number "',g24.4,'" into a string because the number',
     $ /7x,'requires',i3,' characters. The string variable',
     $ ' has only ',i3,' characters.')
        stop
      endif
c
c     Load the string.
c
      ustr = ux(1:j2)
c
      end
