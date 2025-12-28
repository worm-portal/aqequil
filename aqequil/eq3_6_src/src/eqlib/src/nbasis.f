      integer function nbasis(nbasp,nbt,nbtmax,ns)
c
c     This subroutine finds the position of the ns-th species in the
c     array nbasp, which defines the basis set. If the species
c     is not in the basis set, nbasis is returned as zero.
c
c     This subroutine is called by:
c
c       EQLIB/flgset.f
c       EQLIB/gcsts.f
c       EQLIB/switch.f
c       EQLIB/swtchk.f
c       EQ3NR/arrsim.f
c       EQ3NR/chkinx.f
c       EQ3NR/eq3nr.f
c       EQ3NR/intbs3.f
c       EQ6/combmb.f
c       EQ6/intbs6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nbasp  = array of indices of species in the basis set
c       nbt    = the number of species in the basis set
c       ns     = the index of a species
c
c     Principal output:
c
c       nbasis = the basis index of the ns-th species, if this species
c                  is in the basis set
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbtmax
c
      integer nbasp(nbtmax)
      integer nbt,ns
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nb,nss
c
c-----------------------------------------------------------------------
c
      nbasis = 0
      do nb = 1,nbt
        nss = nbasp(nb)
        if (ns .eq. nss) then
          nbasis = nb
          go to 999
        endif
      enddo
c
  999 continue
      end
