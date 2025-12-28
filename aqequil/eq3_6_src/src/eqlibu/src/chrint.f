      subroutine chrint(ivar,nttyo,qrderr,ustr)
c
c     This subroutine reads the integer ivar from the character string
c     ustr. The string may contain non-blank characters other than the
c     integer. If there is more than one integer in the string, only the
c     first is read. If the string is completely blank, the integer is
c     returned with a value of zero. If no integer is present and at
c     least one blank is present, the integer is returned with a value
c     of zero.
c
c     The length of the string variable is unknown. A buffer variable
c     which is employed has a character length of ichpar. This limits
c     the size of the substring for the integer.
c
c     This subroutine is called by:
c
c       Any
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       nttyo  = unit number of the screen file
c       ustr   = the input string variable
c
c     Output:
c
c       ivar   = the output integer
c       qrderr = logical flag, .true. if a read format error occurred
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nttyo,ivar
c
      logical qrderr
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
      integer i,ifirst,ipos1,istr,ichars,j2,nchars
c
      integer ifnobl,ilnobl
c
      character*(ichpar) ux
      character*1 uc
c
      logical qsign
c
c-----------------------------------------------------------------------
c
      qrderr = .false.
c
c     Find the position of the first non-blank character in the
c     buffer variable.
c
      ipos1 = ifnobl(ustr)
c
      if (ipos1 .le. 0) then
        ivar = 0
        go to 999
      endif
c
c     Get the length of the string variable.
c
      nchars = len(ustr)
c
c     Find the first integer substring. This must begin with one of the
c     following characters: +-0123456789. The substring is considered
c     found when a number is found. If the substring begins with a
c     "+" or "-" sign, intervening blanks are permitted.
c
      ifirst = 0
      do i = ipos1,nchars
        uc(1:1) = ustr(i:i)
        if (uc.ge.'0' .and. uc.le.'9') then
          if (ifirst .eq. 0) ifirst = i
          go to 100
        elseif (uc .eq. '-') then
          ifirst = i
        elseif (uc .eq. '+') then
          ifirst = i
        elseif (uc .ne. ' ') then
          ifirst = 0
        endif
      enddo
c
c     No integer is present in the string. See if there is at
c     least one blank present.
c
      i = index(ustr,' ')
      if (i .gt. 0) then
        ivar = 0
        go to 999
      endif
c
      j2 = ilnobl(ustr)
      write (nttyo ,1000) ustr(1:j2)
 1000 format(/" * Error - (EQLIBU/chrint) Can't find an integer",
     $ /7x,'to read from the string "',a,'".')
      qrderr = .true.
      go to 999
c
c     Copy the integer field into the buffer string.
c
c       ichars = length of the substring containing the integer
c       qsign  = logical flag denoting the reading of a substring
c                containing an initial '-' or '+' sign followed
c                by zero or more blanks
c
  100 ichars = 0
      qsign = .false.
c
      do istr = ifirst,nchars
        uc(1:1) = ustr(istr:istr)
        if (uc.ge.'0' .and. uc.le.'9') then
c
c         Have encountered a number character in the integer substring.
c
          ichars = ichars + 1
          ux(ichars:ichars) = uc
          qsign = .false.
        elseif (uc.eq.'-' .or. uc.eq.'+') then
c
c         Have encountered a '-' or '+' sign which begins the
c         integer substring.
c
          ichars = ichars + 1
          ux(ichars:ichars) = uc
          qsign = .true.
        elseif (qsign .and. uc.eq.' ') then
c
c         Have encountered a blank between the '-' or '+' sign and
c         the first non-blank character and the integer string.
c         Allow it in the substring.
c
          ichars = ichars + 1
          ux(ichars:ichars) = uc
        else
          go to 110
        endif
      enddo
c
c     Load blank fill.
c
  110 do i = ichars + 1,ichpar
        ux(i:i) = ' '
      enddo
c
c     Read the integer from the string buffer.
c
      read (ux,1010,err=995) ivar
c
c     The format below should specify an integer field which matches
c     the character length of the buffer variable (ichpar).
c
 1010 format(i16)
      go to 999
c
  995 j2 = ilnobl(ux)
      write (nttyo ,1020) ux(1:j2)
 1020 format(/" * Error - (EQLIBU/chrint) Can't read an integer from",
     $ /7x,'the buffer string "',a,'" due to a formatted read error.')
      qrderr = .true.
c
  999 continue
      end
