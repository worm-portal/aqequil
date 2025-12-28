      subroutine grdxsp(nbasp,nbt,nbtmax,nct,ndrsr,noutpt,
     $ nrdxsp,nstmax,nttyo,uspec)
c
c     This subroutine finds the redox basis species (nrdxsp).
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
c       nbasp  = array of indices of basis species
c
c     Principal output:
c
c       nrdxsp = species index of the redox basis species (= 0 if
c                there is none)
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
      integer nbasp(nbtmax),ndrsr(2,nstmax)
      integer nbt,nct,nrdxsp
c
      character(len=48) uspec(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen,nb,nr1,nr2,ns,nt1
c
      character(len=56) uspn56
c
c-----------------------------------------------------------------------
c
c     The redox basis species, if it exists, folllows the other strict
c     basis species, each of which must associate one-to-one with a
c     chemical element. The strategy is to test the basis species after
c     those to see if it is a strict basis species (has no associated
c     real reaction). If that is the case, then it is the redox species
c     (species index nrdxsp). If the redox species doesn't exist, nrdxsp
c     is returned as zero.
c
      nrdxsp = 0
c
c     Get the candidate basis index (nb).
c
      nb = nct + 1
      if (nb .le. nbt) then
c
c       This index is in range.
c
        ns = nbasp(nb)
        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)
        nt1 = nr2 - nr1 + 1
        if (nt1 .lt. 2) then
c
c         The candidate species is a strict basis species.
c         It must be the redox species.
c
          nrdxsp = ns
        endif
      endif
c
      if (nrdxsp .gt. 0) then
c
c       Calling sequence substitutions:
c         uspec(nrdxsp) for unam48
c
        call fmspnx(jlen,uspec(nrdxsp),uspn56)
c
        write (noutpt,1000) uspn56(1:jlen)
        write (nttyo,1000) uspn56(1:jlen)
 1000   format(/' The redox basis species is ',a,'.')
      else
        write (noutpt,1010)
        write (nttyo,1010)
 1010   format(' No redox basis species was found.')
      endif
c
      end
