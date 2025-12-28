      subroutine zsrt(izmax,narn1,narn2,nstmax,zchar)
c
c     This subroutine finds the largest absolute value of the electrical
c     charge of any aqueous species (izmax). It is used as a limit in
c     calculating higher-order electrostatic terms in Pitzer's
c     equations.
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
c       zchar  = array of electrical charge numbers
c       narn1  = start of the range of aqueous species
c       narn2  = end of the range of aqueous species
c
c     Principal output:
c
c       izmax  = max norm of the electrical charges of the
c                aqueous species
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
      integer izmax,narn1,narn2
c
      real*8 zchar(nstmax)
c
c----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ileft,jzmax,jzmin,ns,nval
c
      real*8 zmax,zmin,zx,zx0,zx1,zx2,zx3,zx4,zx5,zx6,zx7
c
c----------------------------------------------------------------------
c
c     Find the extreme values of electrical charge among the aqueous
c     species. Note that the loop is unrolled.
c
      zmax = 0.
      zmin = 0.
      nval = narn2 - narn1 + 1
      ileft = (nval/8)*8 + narn1 - 1
c
      do ns = narn1,ileft,8
        zx0 = zchar(ns)
        zx1 = zchar(ns + 1)
        zx2 = zchar(ns + 2)
        zx3 = zchar(ns + 3)
        zx4 = zchar(ns + 4)
        zx5 = zchar(ns + 5)
        zx6 = zchar(ns + 6)
        zx7 = zchar(ns + 7)
        zmin = min(zx0,zx1,zx2,zx3,zx4,zx5,zx6,zx7,zmin)
        zmax = max(zx0,zx1,zx2,zx3,zx4,zx5,zx6,zx7,zmax)
      enddo
c
      do ns = ileft + 1,narn2
        zx = zchar(ns)
        zmin = min(zx,zmin)
        zmax = max(zx,zmax)
      enddo
c
c     Find the largest absolute value of electrical charge.
c
      jzmin = nint(zmin)
      jzmax = nint(zmax)
      izmax = max(-jzmin,jzmax)
c
      end
