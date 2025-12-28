      subroutine chreal(nttyo,qrderr,ustr,var)
c
c     This subroutine reads the real*8 var from the character string
c     ustr. The string may contain non-blank characters other than the
c     real*8 number. If there is more than one real*8 number in the
c     string, only the first is read. If the string is completely blank,
c     the real*8 number is returned with a value of zero. If no real*8
c     number is present and at least one blank is present, the
c     real*8 number is returned with a value of zero.
c
c     The length of the string variable is unknown. A buffer variable
c     which is employed has a character length of ichpar. This limits
c     the size of the substring for the real*8 number.
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
c       var    = the output real*8 number
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
      integer nttyo
c
      logical qrderr
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
      integer i,iexpf,ifirst,ilettr,inum1,ipoint,ipos1,istr,
     $ ichars,j2,nchars
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
        var = 0.
        go to 999
      endif
c
c     Get the length of the string variable.
c
      nchars = len(ustr)
c
c     Find the first real*8 substring. This must begin with one of the
c     following characters: .+-0123456789. The substring is considered
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
        elseif (uc .eq. '.') then
          if (ifirst .eq. 0) ifirst = i
        elseif (uc .ne. ' ') then
          ifirst = 0
        endif
      enddo
c
c     No real*8 number is present in the string. See if there is at
c     least one blank present.
c
      i = index(ustr,' ')
      if (i .gt. 0) then
        var = 0.
        go to 999
      endif
c
      j2 = ilnobl(ustr)
      write (nttyo,1000) ustr(1:j2)
 1000 format(/" * Error - (EQLIBU/chreal) Can't find a real*8 number",
     $ /7x,'to read from the string "',a,'".')
      qrderr = .true.
      go to 999
c
c     Copy the real*8 field into the buffer string.
c
c       ipoint = original string position containing '.'
c       iexpf  = original string position of the first character
c                  in the exponent field, normally 'e','E', 'd',
c                  or 'D', but possibly '-' or '+'
c       ilettr = original string position containing 'e','E',
c                  'd', or 'D'
c       inum1  = number of number characters in the non-exponent
c                  field
c       ichars = length of the substring containing the real*8
c                  number
c       qsign  = logical flag denoting the reading of a substring
c                containing an initial '-' or '+' sign followed
c                by zero or more blanks
c
  100 ipoint = 0
      iexpf = 0
      ilettr = 0
      inum1 = 0
      ichars = 0
      qsign = .false.
c
      do istr = ifirst,nchars
        uc(1:1) = ustr(istr:istr)
        if (uc.ge.'0' .and. uc.le.'9') then
c
c         Have encountered a number character in either the non-exponent
c         or exponent field.
c
          ichars = ichars + 1
          ux(ichars:ichars) = uc
          if (iexpf .eq. 0) inum1 = inum1 + 1
          qsign = .false.
        elseif (uc.eq.'-' .or. uc.eq.'+') then
          if (iexpf .eq. 0) then
            if (inum1 .eq. 0) then
c
c             Have encountered a '-' or '+' sign which begins the
c             non-exponent field.
c
              ichars = ichars + 1
              ux(ichars:ichars) = uc
              qsign = .true.
            else
c
c             Have encountered a '-' or '+' sign which starts
c             the exponent field.
c
              ichars = ichars + 1
              ux(ichars:ichars) = uc
              iexpf = istr
              qsign = .false.
            endif
          else
            if (istr .eq. (ilettr + 1)) then
c
c             Have a '-' or '+' sign following an 'e', 'E',
c             'd', or 'D' (in the exponent field).
c
              ichars = ichars + 1
              ux(ichars:ichars) = uc
              qsign = .false.
            else
c
c             Have encountered a second '-' or '+' sign in the
c             exponent field. Terminate the substring without it.
c
              go to 110
            endif
          endif
        elseif (uc.eq.' ') then
          if (qsign) then
c
c           Have encountered a blank between the initial '-' or
c           '+' sign and the first non-blank character in the
c           substring. Allow it in the substring.
c
            ichars = ichars + 1
            ux(ichars:ichars) = uc
          else
c
c           Have encountered an illegal blank. Terminate the
c           substring without it.
c
            go to 110
          endif
        elseif (uc .eq. '.') then
          if (iexpf .eq. 0) then
            if (ipoint .eq. 0) then
              ichars = ichars + 1
              ux(ichars:ichars) = uc
              ipoint = istr
              qsign = .false.
            else
c
c             Have encountered a second '.' in the non-exponent
c             field. Terminate the substring without it.
c
              go to 110
            endif
          else
c
c           Have encountered a '.' in the exponent field.
c           Terminate the substring without it.
c
            go to 110
          endif
        elseif (uc.eq.'e' .or. uc.eq.'E' .or. uc.eq.'d'
     $    .or. uc.eq.'D') then
          if (iexpf .eq. 0) then
c
c           Have encountered an exponent character.
c
            ichars = ichars + 1
            ux(ichars:ichars) = uc
            ilettr = istr
            iexpf = istr
            qsign = .false.
          else
c
c           Have encountered a second 'e', 'E', 'd', or 'D' in
c           the exponent field. Terminate the substring without it.
c
            go to 110
          endif
        else
c
c         Have encountered an illegal character. Terminate the
c         substring without it.
c
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
      if (inum1 .le. 0) then
        j2 = ilnobl(ux)
        write (nttyo,1010) ux(1:j2)
 1010   format(/" * Warning - (EQLIBU/chreal) The real*8 number in",
     $  /7x,'the buffer string "',a,'" has no numbers in the',
     $  /7x,'non-exponent part.')
      endif
c
c     Read the real*8 from the string buffer.
c
      read (ux,1030,err=995) var
c
c     The format below should specify a real*8 field which matches
c     the character length of the buffer variable (ichpar).
c
 1030 format(g24.0)
      go to 999
c
  995 j2 = ilnobl(ux)
      write (nttyo,1040) ux(1:j2)
 1040 format(/" * Error - (EQLIBU/chreal) Can't read a real*8 number",
     $ /7x,'from the buffer string "',a,'" due to a formatted read',
     $ ' error.')
      qrderr = .true.
c
  999 continue
      end
