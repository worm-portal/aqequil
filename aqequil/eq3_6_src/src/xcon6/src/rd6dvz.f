      subroutine rd6dvz(ninpt,ucode,urelno,ustage,ueqlrn,ueqlst)
c
c     This subroutine looks at the comments on the old EQ6 INPUT file
c     in menu-style ("D") format, and tries to recover the code and
c     version number data normally written on the bottom half of the
c     PICKUP file to identify the code and version which wrote this
c     part of that file.
c
c     This subroutine is called by:
c
c       XCON6/xcon6.f
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ninpt
c
      character*8 ucode,urelno,ustage,ueqlrn,ueqlst
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j
c
      character*80 uline
c
c-----------------------------------------------------------------------
c
      ucode  = ' '
      urelno = ' '
      ustage = ' '
      ueqlrn = ' '
      ueqlst = ' '
      rewind(ninpt)
c
      do 100 i = 1,10000
        read (ninpt,'(a)',end=999) uline
        if (uline(1:1).eq.'c' .or. uline(1:1).eq.'*') then
          j = index(uline,'pickup file written by')
          if (j .le. 0) j = index(uline,'Pickup file written by')
          if (j .gt. 0) go to 110
        endif
  100 continue
      go to 999
c
  110 j = index(uline,'eq3nr.')
      if (j .eq. 0) j = index(uline,'EQ3NR.')
      if (j .gt. 0) then
        ucode = 'eq3nr'
        urelno= uline(j + 6:j + 9)
        ustage= uline(j + 10:j+ 13)
        go to 120
      endif
      j = index(uline,'EQ3NR, version ')
      if (j .gt. 0) then
        ucode = 'EQ3NR'
        urelno= uline(j + 15:j + 18)
        ustage= uline(j + 21:j+ 28)
        go to 120
      endif
      j = index(uline,'eq6.')
      if (j .eq. 0) j = index(uline,'EQ6.')
      if (j .gt. 0) then
        ucode = 'eq6'
        urelno= uline(j + 4:j + 7)
        ustage= uline(j + 8:j + 11)
        go to 120
      endif
      j = index(uline,'EQ6, version ')
      if (j .gt. 0) then
        ucode = 'EQ6'
        urelno= uline(j + 13:j + 16)
        ustage= uline(j + 19:j + 26)
        go to 120
      endif
c
  120 read (ninpt,'(a)',end=999) uline
      j = index(uline,'supported by eqlib.')
      if (j .eq. 0) j = index(uline,'supported by EQLIB.')
      if (j .gt. 0) then
        ueqlrn = uline(j + 19:j + 22)
        ueqlst = uline(j + 23:j + 26)
        go to 999
      endif
      j = index(uline,'supported by EQLIB, version ')
      if (j .gt. 0) then
        ueqlrn = uline(j + 28:j + 31)
        ueqlst = uline(j + 34:j + 41)
        go to 999
      endif
      j = index(uline,'Supported by EQLIB, version ')
      if (j .gt. 0) then
        ueqlrn = uline(j + 28:j + 31)
        ueqlst = uline(j + 34:j + 41)
        go to 999
      endif
c
  999 continue
      end
