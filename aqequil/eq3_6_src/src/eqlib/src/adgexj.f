      subroutine adgexj(je,noutpt,nttyo,ugexjd)
c
c     This subroutine assigns a default name for the je-th site of a
c     generic exchange phase. The general model is "S(n)", where n
c     is the site number. However, if n is larger than 99999, the S
c     and both parentheses are all dropped.
c
c     This subroutine is called by:
c
c       EQLIB/intexi.f
c       EQ3NR/intge3.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       je     = site index
c
c     Principal output:
c
c       ugexjd = default name of the je-th site of a generic ion
c                  exchange phase
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
      integer je
c
      character*8 ugexjd
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c     None
c
c-----------------------------------------------------------------------
c
      ugexjd = ' '
      ugexjd(1:2) = 'S('
      if (je .le. 9) then
        write (ugexjd(3:3),'(i1)') je
        ugexjd(4:4) = ')'
      elseif (je .le. 99) then
        write (ugexjd(3:4),'(i2)') je
        ugexjd(5:5) = ')'
      elseif (je .le. 999) then
        write (ugexjd(3:5),'(i3)') je
        ugexjd(6:6) = ')'
      elseif (je .le. 9999) then
        write (ugexjd(3:6),'(i4)') je
        ugexjd(7:7) = ')'
      elseif (je .le. 99999) then
        write (ugexjd(3:7),'(i5)') je
        ugexjd(8:8) = ')'
      elseif (je .le. 999999) then
        ugexjd = 'S '
        write (ugexjd(2:7),'(i6)') je
      elseif (je .le. 9999999) then
        ugexjd = 'S '
        write (ugexjd(2:8),'(i7)') je
      elseif (je .le. 99999999) then
        ugexjd = ' '
        write (ugexjd(1:8),'(i8)') je
      else
        write (noutpt,1000) je
        write (nttyo,1000) je
 1000   format (/' * Error - (EQLIB/adgexj) Programming error trap:',
     $  ' The site',/7x,'number ',i12,' for a generic ion exchange',
     $  ' phase',/7x,'is too large to fit in the 8-character field',
     $  ' reserved for it',/7x,'in the default site name.')
        stop
      endif
c
      end
