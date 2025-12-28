      subroutine g1dat(ier,noutpt,nttyo,udastr,var)
c
c     This subroutine reads a number (var) from a string (udastr)
c
c     This subroutine is called by:
c
c       EQPT/gnenb.f
c       EQPT/rdpca.f
c       EQPT/rdpth.f
c       EQPT/rdpni.f
c       EQPT/rdpn2.f
c       EQPT/rdpnn.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       udastr = a string containing one numerical data field
c
c     Principal output:
c
c       ier    = error flag
c       var    = the number contained in that field
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
      integer ier
c
      character(len=80) udastr
c
      real(8) var
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jj,jlen
c
      integer ilnobl
c
c-----------------------------------------------------------------------
c
      var = 0.
      ier = 0
c
      call lejust(udastr)
      jj = index(udastr,' ')
      if (jj .eq. 0) jj = 81
      if (jj .lt. 80) udastr(jj:80) = ' '
      jlen = ilnobl(udastr)
c
      if (jlen .gt. 25) then
        write (noutpt,1000) udastr(1:jlen)
        write (nttyo,1000) udastr(1:jlen)
 1000   format(/' * Error - (EQPT/g1dat) Have found a data',
     $  ' field in the string:',/7x,'"',a,'"',/7x,'that appears',
     $  ' to exceed the allowed 25 characters.')
        ier = 1
      else
        read (udastr,1010,err=995) var
 1010   format(e25.18)
      endif
      go to 999
c
  995 ier = 1
c
  999 continue
      end
