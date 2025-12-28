      subroutine intmtb(mtb,mtbaq,mtbaqi,mtbi,nbasp,nbt,nbti,nbtmax,
     $ noutpt,nstmax,nttyo,ubmtbi,uspec)
c
c     This subroutine interprets the mass balance totals read from the
c     input file. It constructs the mtb and mtbaq arrays.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       mtbi   = total number of moles of basis species in the
c                  equilibrium system
c       mtbaqi = total number of moles of basis species in the
c                  aqueous phase
c       ubmtbi = array of names of the corresponding species
c
c     Principal output:
c
c       mtb    = total number of moles of basis species in the
c                  equilibrium system
c       mtbaq  = total number of moles of basis species in the
c                  aqueous phase
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
      integer nbasp(nbtmax)
c
      integer nbt,nbti
c
      character(len=48) ubmtbi(nbtmax),uspec(nstmax)
c
      real(8) mtb(nbtmax),mtbaq(nbtmax),mtbaqi(nbtmax),mtbi(nbtmax)
c
c----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen2,nb,nbi,nerr,ns
c
      character(len=56) uspn56
c
c----------------------------------------------------------------------
c
      nerr = 0
c
      do nbi = 1,nbti
        do nb = 1,nbt
          ns = nbasp(nb)
          if (ubmtbi(nbi)(1:48) .eq. uspec(ns)(1:48)) then
            mtb(nb) = mtbi(nbi)
            mtbaq(nb) = mtbaqi(nbi)
            go to 100
          endif
        enddo
c
        call fmspnm(jlen2,ubmtbi(nbi),uspn56)
        write (noutpt,1000) uspn56(1:jlen2)
        write (nttyo,1000) uspn56(1:jlen2)
 1000   format(/' * Error - (EQ6/intmtb) A mass balance is defined',
     $  ' on the input',/7x,'file for ',a,", but this species isn't",
     $  ' in the',/7x,"currently active basis set. Either it isn't",
     $  ' on the current data file',/7x,'or it has been suppressed',
     $  ' as by an nxmod or iopt(15) option.')
        nerr = nerr + 1
  100   continue
      enddo
c
      if (nerr .gt. 0) stop
c
      end
