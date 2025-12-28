      subroutine intbsw(nbasp,nbaspx,nbt,nbtmax,nobswt,noutpt,nst,
     $ nstmax,nttyo,uobsw,uspec)
c
c     This subroutine interprets ordinary basis switching directives
c     specified on the input file.
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
c       nbasp  = the indices of the species in the basis set
c                  (modified by this subroutine)
c       nbt    = the number of species in the basis set
c       nobswt = the number of basis switches to make
c       uobsw  = the name pairs for the species to be switched
c       uspec  = array of species names
c
c     Principal output:
c
c       nbasp  = the indices of the species in the basis set
c                  (modified by this subroutine)
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
      integer nbasp(nbtmax),nbaspx(nbtmax)
      integer nobswt,nbt,nst
c
      character(len=48) uobsw(2,nbtmax),uspec(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen,n,nb,nerr,ns1,ns2
c
      character(len=56) uspn56
      character(len=48) unam48
c
c-----------------------------------------------------------------------
c
      nerr = 0
c
      do n = 1,nobswt
c
c       Find the first species.
c
        do nb = 1,nbt
          ns1 = nbaspx(nb)
          if (uspec(ns1)(1:48) .eq. uobsw(1,n)(1:48)) go to 110
        enddo
c
        nerr = nerr + 1
c
        unam48 = uobsw(1,n)
        call fmspnx(jlen,unam48,uspn56)
        write (noutpt,1000) uspn56(1:jlen)
        write (nttyo,1000) uspn56(1:jlen)
 1000   format(/' * Error - (EQLIB/intbsw) The species ',a,
     $  /7x,'is specified to be replaced in an ordinary basis switch,',
     $  /7x,'but it is not in the active basis set.')
        go to 150
c
  110   continue
c
c       Find the second species.
c
        do ns2 = 1,nst
          if (uspec(ns2)(1:48) .eq. uobsw(2,n)(1:48)) go to 130
        enddo
c
        nerr = nerr + 1
c
        unam48 = uobsw(2,n)
        call fmspnx(jlen,unam48,uspn56)
        write (noutpt,1010) uspn56(1:jlen)
        write (nttyo,1010) uspn56(1:jlen)
 1010   format(/' * Error - (EQLIB/intbsw) The species ',a,
     $  /7x,'is specified to be put into the basis set by an ordinary',
     $  /7x,'basis switch but it is not in the set of active',
     $  ' species.')
        go to 150
c
  130   continue
c
c       Mark to switch.
c
        nbasp(nb) = ns2
c
  150   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Stop if name match errors have been encountered.
c
      if (nerr .gt. 0) then
        write (noutpt,1020) nerr
        write (nttyo,1020) nerr
 1020   format(/' * Error - (EQLIB/intbsw) ',i4,' errors were',
     $  /7x,'encountered interpreting ordinary basis switches that',
     $  /7x,'were directed on the input file.')
        stop
      endif
c
      end
