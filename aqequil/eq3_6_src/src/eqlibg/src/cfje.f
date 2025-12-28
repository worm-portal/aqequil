      subroutine cfje(conc,fjec,jcsort,narn1,narn2,nstmax,zchcu6)
c
c     This subroutine calculates the ionic asymmetry (the 3-rd order
c     electrostatic moment function J). This is defined as:
c
c       J = 1/6 SUM(i) m(i)z(i)**3
c
c     This is a higher order analogue of the ionic strength, I, which
c     is the 2nd-order electrostatic moment function. For comparison,
c     the ionic strength is defined as:
c
c       I = 1/2 SUM(i) m(i)z(i)**2
c
c     This J is not to be confused with the functions J0(x), J1(x), and
c     J2(x), which are also involved in higher-order electrostatic
c     contributions to activity coefficients (see EQLIBG/ghj0.f).
c     Note that a sorted summation is used.
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
c       zchcu6 = array of one-sixth the electrical charge cubed
c
c     Principal output:
c
c       fjec   = the ionic asymmetry (the 3rd-order electrostatic
c                  moment function J)
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
      real*8 conc(nstmax),zchcu6(nstmax)
      real*8 fjec
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
c       fjec = 0.
c       do ns = narn1 + 1,narn2
c         fjec = fjec + conc(ns)*zchcu6(ns)
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
        xx = xx + conc(ns0)*zchcu6(ns0) + conc(ns1)*zchcu6(ns1)
     $          + conc(ns2)*zchcu6(ns2) + conc(ns3)*zchcu6(ns3)
     $          + conc(ns4)*zchcu6(ns4) + conc(ns5)*zchcu6(ns5)
     $          + conc(ns6)*zchcu6(ns6) + conc(ns7)*zchcu6(ns7)
      enddo
c
      do nss = ileft + 1,narn2
        ns = jcsort(nss)
        xx = xx + conc(ns)*zchcu6(ns)
      enddo
c
      fjec = xx
c
      end
