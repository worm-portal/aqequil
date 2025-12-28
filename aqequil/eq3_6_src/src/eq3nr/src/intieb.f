      subroutine intieb(iebal,iebal3,ier,jsflag,nbasp,nbt,nbtmax,noutpt,
     $ nstmax,nttyo,uebal,uspec,zchar)
c
c     This subroutine finds the basis index (iebal) of the species to be
c     adjusted for electrical balance, if any.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
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
      integer jsflag(nstmax),nbasp(nbtmax)
      integer iebal,iebal3,ier,nbt,noutpt,nttyo
c
      character*48 uspec(nstmax)
      character*24 uebal
c
      real(8) zchar(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,nb,nebal,ns
c
      integer ilnobl
c
c-----------------------------------------------------------------------
c
c     if iebal3 is 0 or iebal = 1 but uebal is 'none' or blank, do no
c     electrical balancing.
c
      iebal = 0
      if (iebal3 .le. 0) uebal(1:24) = 'None'
      if (uebal(1:8) .eq. '        ') uebal(1:24) = 'None'
      if (uebal(1:8) .eq. 'none    ') uebal(1:24) = 'None'
      if (uebal(1:8) .eq. 'NONE    ') uebal(1:24) = 'None'
      if (uebal(1:5) .eq. 'None ') go to 999
c
      j2 = ilnobl(uebal)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find the basis index of the species specified on the input file.
c
      do nb = 1,nbt
        ns = nbasp(nb)
        if (uebal(1:24).eq.uspec(ns)(1:24)) then
          iebal = nb
          nebal = ns
          go to 100
        endif
      enddo
c
      write (noutpt,1000) uebal(1:j2)
      write (nttyo,1000) uebal(1:j2)
 1000 format(/' * Error - (EQ3NR/intieb) The input file specifies',
     $ /7x,'that the concentration of ',a,' is to be adjusted',
     $ /7x,'to achieve electrical balance. However, this species',
     $ /7x,"isn't in the basis set.")
      ier = 1
  100 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Test the specified ion to see if it is okay.
c
c     Does the ion selected for electrical balancing have charge?
c
      if (zchar(nebal) .eq. 0.) then
        write (noutpt,1010) uebal(1:j2)
        write (nttyo,1010) uebal(1:j2)
 1010   format(/' * Warning - (EQ3NR/intieb) The input file specifies',
     $  /7x,'that the concentration of ',a,' is to be adjusted',
     $  /7x,'to achieve electrical balance. However, this species',
     $  /7x,"doesn't have an electrical charge. The success of this",
     $  /7x,'calculation will depend on the concentration adjustment',
     $  /7x,'influencing the concentration of one or more charged',
     $  /7x,'dependent species.')
      endif
c
c     Is the ion selected for electrical balancing in the model?
c
      if (jsflag(nebal) .gt. 0) then
        write (noutpt,1020) uebal(1:j2)
        write (nttyo,1020) uebal(1:j2)
 1020   format(/' * Error - (EQ3NR/intieb) The input file specifies',
     $  /7x,'that the concentration of ',a,' is to be adjusted',
     $  /7x,'to achieve electrical balance. However, this species',
     $  /7x,"is effectively suppressed and thus can't be present",
     $  /7x,'in the model.')
        ier = 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
