      subroutine intsbs(nb1,nb2,nbaspd,nbtd,nbtmax,noutpt,ns1,ns2,
     $ nsbsw,nstmax,nttyo,usbsw,uspeca)
c
c     This subroutine interprets the nsbsw-th directive read from the
c     input file for effecting special basis switching. The roles of
c     a strict basis species and an auxiliary basis species are
c     switched by this process. The name of the former species is
c     usbsw(1,nsbsw), that of the latter is usbsw(2,nsbsw). Here nb1
c     is the basis index of the strict basis species and nb2 is that
c     of the auxiliary basis species. Here also ns1 and ns2 are the
c     corresponding species indices.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbtmax,nstmax
c
      integer noutpt,nttyo
c
      integer nbaspd(nbtmax)
      integer nb1,nb2,nbtd,ns1,ns2,nsbsw
c
      character*48 usbsw(2,nbtmax),uspeca(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen,j2,nerr
c
      integer ilnobl
c
      character*48 ublk48,unam48
      character*56 uspn56
      character*8 ux8
c
c-----------------------------------------------------------------------
c
      data ublk48/'                                                '/
c
c-----------------------------------------------------------------------
c
      nerr = 0
c
      do nb1 = 1,nbtd
        ns1 = nbaspd(nb1)
        if (usbsw(1,nsbsw)(1:48) .eq. uspeca(ns1)(1:48)) go to 120
      enddo
c
      nerr = nerr + 1
c
      unam48 = usbsw(1,nsbsw)
      call fmspnx(jlen,unam48,uspn56)
      write (noutpt,1000) uspn56(1:jlen)
      write (nttyo,1000) uspn56(1:jlen)
 1000 format(/' * Error - (EQLIB/intsbs) Invalid special basis',
     $ ' directive on the',/7x,'input file: ',a,' is not in the basis',
     $ ' set and therefore',/7x,'can not be specified in a special',
     $ ' basis switch.')
c
  120 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (usbsw(2,nsbsw)(1:48) .eq. ublk48(1:48)) then
        nb2 = nb1
        ns2 = ns1
      else
        do nb2 = 1,nbtd
          ns2 = nbaspd(nb2)
          if (usbsw(2,nsbsw)(1:48) .eq. uspeca(ns2)(1:48)) go to 140
        enddo
c
        nerr = nerr + 1
c
        unam48 = usbsw(2,nsbsw)
        call fmspnx(jlen,unam48,uspn56)
        write (noutpt,1000) uspn56(1:jlen)
        write (nttyo,1000) uspn56(1:jlen)
c
  140   continue
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Stop if name match errors have been encountered.
c
      if (nerr .gt. 0) then
        write (ux8,'(i5)') nerr
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1020) ux8(1:j2)
        write (nttyo,1020) ux8(1:j2)
 1020   format(/' * Error - (EQLIB/intsbs) ',a,' errors were',
     $  /7x,'encountered in interpreting special basis switches',
     $  /7x,'specified on the input file.')
        stop
      endif
c
      end
