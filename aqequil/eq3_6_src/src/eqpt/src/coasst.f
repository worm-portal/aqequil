      subroutine coasst(jassan,jassca,jassne,nat,natmax,uaqsp,zaqsp)
c
c     Get the numbers of aqueous solute cations, anions, and neutral
c     species. These will be used to generate lists of species pairs
c     and triplets for use with Pitzer's equations.
c
c     The fictive redox species aqueous e- and aqueous O2(g), if
c     present, are excluded in computing these totals. So is solvent
c     water.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nat    = the number of aqueous species
c       uaqsp  = array of names of aqueous species
c       zaqsp  = array of electrical charages of aqueous species
c
c     Principal output:
c
c       jassca = number of aqueous cation species
c       jassan = number of aqueous anion species (excluding aqueous e-)
c       jassne = number of aqueous neutral species (excluding solvent
c                  water and aqeuous O2(g))
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer natmax
c
      integer jassan,jassca,jassne,nat
c
      character(len=24) uaqsp(natmax)
c
      real*8 zaqsp(natmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n
c
c-----------------------------------------------------------------------
c
c     Count the aqueous solute cations, anions, and neutrals. Include
c     in the counts only real solute species.
c
      jassne = 0
      jassca = 0
      jassan = 0
      do n = 2,nat
        if (zaqsp(n) .eq. 0.) then
          if (uaqsp(n)(1:6) .ne. 'O2(g) ') then
            jassne = jassne + 1
          endif
        elseif (zaqsp(n) .gt. 0.) then
          jassca = jassca + 1
        else
          if (uaqsp(n)(1:3) .ne. 'e- ') then
            jassan = jassan + 1
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
