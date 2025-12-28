      subroutine cfxi(conc,fxic,jcsort,narn1,narn2,nstmax,zchsq2)
c
c     This subroutine calculates the ionic strength. Note that a sorted
c     summation is used.
c
c     This subroutine is called by:
c
c       EQLIB/ngcadv.f
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
c       narn1  = start of the range of aqueous species
c       narn2  = end of the range of aqueous species
c       zchsq2 = array of one-half the electrical charge squared
c
c     Principal output:
c
c       fxic   = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
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
      real*8 conc(nstmax),zchsq2(nstmax)
      real*8 fxic
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ileft,ns,nss,ns0,ns1,ns2,ns3,ns4,ns5,ns6,ns7,nval
c
      real*8 xx
c
c-----------------------------------------------------------------------
c
c     Logically, this could be represented by:
c
c       fxic = 0.
c       do ns = narn1 + 1,narn2
c         fxic = fxic + conc(ns)*zchsq2(ns)
c       enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     In doing a sorted summation, a complication arises in that
c     water (ns = narn1) must be included in the range of the
c     loop, because the jcsort array is organized to contain the
c     species indices sorted according to increasing concentration
c     within phase ranges. There is no available sorting among
c     aqueous solute species only. Because the electrical charge
c     of water is zero, the calculation can be made looping from
c     narn1 to narn2. It is not necessary to temporarily set the
c     concentration of water (conc(narn1)) to zero. Note the use
c     of a local variable (xx) within the loop. Note also that the
c     loop is unrolled.
c
      xx = 0.
      nval = narn2 - narn1 + 1
      ileft = (nval/8)*8 + narn1 - 1
c
      do nss = narn1,ileft,8
        ns0 = jcsort(nss)
        ns1 = jcsort(nss + 1)
        ns2 = jcsort(nss + 2)
        ns3 = jcsort(nss + 3)
        ns4 = jcsort(nss + 4)
        ns5 = jcsort(nss + 5)
        ns6 = jcsort(nss + 6)
        ns7 = jcsort(nss + 7)
        xx = xx + conc(ns0)*zchsq2(ns0) + conc(ns1)*zchsq2(ns1)
     $          + conc(ns2)*zchsq2(ns2) + conc(ns3)*zchsq2(ns3)
     $          + conc(ns4)*zchsq2(ns4) + conc(ns5)*zchsq2(ns5)
     $          + conc(ns6)*zchsq2(ns6) + conc(ns7)*zchsq2(ns7)
      enddo
c
      do nss = ileft + 1,narn2
        ns = jcsort(nss)
        xx = xx + conc(ns)*zchsq2(ns)
      enddo
c
      fxic = xx
c
      end
