      subroutine csigm(conc,jcsort,narn1,narn2,nstmax,sigmmc)
c
c     This subroutine calculates the sum of the molalities of aqueous
c     solute species (sigmmc):
c
c       Sigma(i) m(i)
c
c     Note that a sorted summation is used.
c
c     This subroutine is called by:
c
c       EQLIB/ngcadv.f
c       EQLIB/ncmpex.f
c       EQ3NR/arrset.f
c       EQ3NR/scripx.f
c       EQ6/scripz.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       conc   = array of species concentrations
c       jcsort = array of species indices, in order of increasing
c                  concentration, but with sorting restricted to within
c                  phase ranges
c       narn1  = start of the range of aqueous species; this is
c                  also the index of the solvent, water
c       narn2  = end of the range of aqueous species
c
c     Principal output:
c
c       sigmmc = the sum of the molalities of the aqueous solute
c                  species
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nstmax
c
      integer jcsort(nstmax)
      integer narn1,narn2
c
      real*8 conc(nstmax)
      real*8 sigmmc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ileft,nss,nval
c
      real*8 cw,sx
c
c----------------------------------------------------------------------
c
c     Logically, this could be represented by:
c
c       sigmmc = 0.
c       do ns = narn1 + 1,narn2
c         sigmmc = sigmmc + conc(ns)
c       enddo
c
c     In doing a sorted summation, a complication arises in that
c     water (ns = narn1) must be included in the range of the
c     loop, because the jcsort array is organized to contain the
c     species indices sorted according to increasing concentration
c     within phase ranges. There is no available sorting among
c     aqueous solute species only. This complication is dealt with
c     by temporarily setting the concentration of water (conc(narn1))
c     to zero. Note the use of a local variable (sx) within the loop.
c     Note also that the loop is unrolled.
c
      sx = 0.
      cw = conc(narn1)
      conc(narn1) = 0.
      nval = narn2 - narn1 + 1
      ileft = (nval/8)*8 + narn1 - 1
c
      do nss = narn1,ileft,8
        sx = sx + conc(jcsort(nss))      + conc(jcsort(nss + 1))
     $          + conc(jcsort(nss + 2))  + conc(jcsort(nss + 3))
     $          + conc(jcsort(nss + 4))  + conc(jcsort(nss + 5))
     $          + conc(jcsort(nss + 6))  + conc(jcsort(nss + 7))
      enddo
c
      do nss = ileft + 1,narn2
        sx = sx + conc(jcsort(nss))
      enddo
c
      sigmmc = sx
      conc(narn1) = cw
c
      end
