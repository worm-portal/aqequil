      subroutine cfxist(ctb,fxistc,nbaspd,nbt,nbtmax,nstmax,zchsq2)
c
c     This subroutine calculates the stoichiometric ionic strength
c     (fxistc). Note that a negative value of total concentration
c     (ctb) for H+ is treated as a positive value for OH- (and vice
c     versa).
c
c     This subroutine is called by:
c
c       EQ3NR/scripx.f
c       EQ6/scripz.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ctb    = array of total concentrations of the basis species in
c                  the 'd' set
c       nbaspd = array containing the indices of the species in the 'd'
c                  basis set
c       nbt    = the number of active basis species
c       zchsq2 =  array of one-half the charge squared
c
c     Principal output:
c
c       fxistc = the stoichiometric ionic strength, defined in terms
c                  of the total molalities of the basis species in the
c                  'd' set
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
      integer nbaspd(nbtmax)
      integer nbt
c
      real*8 ctb(nbtmax),zchsq2(nstmax)
      real*8 fxistc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nb,ns
c
c-----------------------------------------------------------------------
c
c     Note: the relatively small likely value of nbt does not justify
c     the use of an unrolled loop.
c
      fxistc = 0.
      do nb = 1,nbt
        ns = nbaspd(nb)
        fxistc = fxistc + zchsq2(ns)*abs(ctb(nb))
      enddo
c
      end
