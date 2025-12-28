      subroutine adgexp(ne,noutpt,nttyo,ugexpd)
c
c     This subroutine assigns a default name to the ne-th generic ion
c     exchange phase. The general model is "Exchanger(n)", where n
c     is the index number ne.
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
c       ne     = index number of the current generic ion exchange phase
c
c     Principal output:
c
c       ugexpd = default name of the ne-th generic ion exchange phase
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
      integer ne
c
      character*24 ugexpd
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c     None
c
c-----------------------------------------------------------------------
c
      ugexpd = ' '
      ugexpd(1:10) = 'Exchanger('
      if (ne .le. 9) then
        write (ugexpd(11:11),'(i1)') ne
        ugexpd(12:12) = ')'
      elseif (ne .le. 99) then
        write (ugexpd(11:12),'(i2)') ne
        ugexpd(13:13) = ')'
      elseif (ne .le. 999) then
        write (ugexpd(11:13),'(i3)') ne
        ugexpd(14:14) = ')'
      elseif (ne .le. 9999) then
        write (ugexpd(11:14),'(i4)') ne
        ugexpd(15:15) = ')'
      elseif (ne .le. 99999) then
        write (ugexpd(11:15),'(i5)') ne
        ugexpd(16:16) = ')'
      elseif (ne .le. 999999) then
        write (ugexpd(11:16),'(i6)') ne
        ugexpd(17:17) = ')'
      elseif (ne .le. 9999999) then
        write (ugexpd(11:17),'(i7)') ne
        ugexpd(18:18) = ')'
      elseif (ne .le. 99999999) then
        write (ugexpd(11:18),'(i8)') ne
        ugexpd(19:19) = ')'
      else
        write (noutpt,1000) ne
        write (nttyo,1000) ne
 1000   format (/' * Error - (EQLIB/adgexp) Programming error trap:',
     $  ' The index',/7x,'number ',i12,' of a generic ion exchange',
     $  ' phase',/7x,'is too large to fit in the 8-character field',
     $  ' reserved for it',/7x,'in the default site name.')
        stop
      endif
c
      end
