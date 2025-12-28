      subroutine cabar(abar,azero,conc,jcsort,fxi,narn1,narn2,natmax,
     $ nstmax,zchsq2)
c
c     This subroutine calculates the average hard core diameter of
c     aqueous ionic species (abar). This average is defined
c     using "ionic strength weighting"; that is, the weighting factor
c     for each species is 1/2 the square of its electrical charge.
c     This thus excludes contributions from electrically neutral
c     solute species. This "abar" is used in empirical ion size
c     mixing models of the first-order Debye-Huckel contribution to
c     aqueous species activity coefficients.
c
c     This subroutine is called by:
c
c       EQLIBG/gcoeff.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       azero  = array of hard core diameters
c       conc   = array of species concentrations
c       jcsort = array of species indices, in order of increasing
c                  concentration, but with sorting restricted to within
c                  phase ranges
c       fxi    = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
c       narn1  = start of the range of aqueous species
c       narn2  = end of the range of aqueous species
c       zchsq2 = array of z(i)**2/2 values
c
c     Principal output:
c
c       abar  = average ionic hard core diameter
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer natmax,nstmax
c
      integer jcsort(nstmax)
c
      integer narn1,narn2
c
      real*8 azero(natmax),conc(nstmax),zchsq2(nstmax)
      real*8 abar,fxi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ileft,nnn,ns,nss,ns0,ns1,ns2,ns3,ns4,ns5,ns6,ns7,nval
c
      real*8 sx
c
c-----------------------------------------------------------------------
c
      abar = 0.
      if (fxi .le. 0.) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Logically, this could be represented by:
c
c       sx = 0.
c       do ns = narn1 + 1,narn2
c         na = ns - narn1 + 1
c         sx = sx + conc(ns)*zchsq2(ns)*azero(na)
c       enddo
c       abar = sx/fxi
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
c     of a local variable (sx) within the loop. Note also that the
c     loop is unrolled.
c
      sx = 0.
      nnn = narn1 - 1
      nval = narn2 - nnn
      ileft = (nval/8)*8 + nnn
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
        sx = sx + conc(ns0)*zchsq2(ns0)*azero(ns0 - nnn)
     $          + conc(ns1)*zchsq2(ns1)*azero(ns1 - nnn)
     $          + conc(ns2)*zchsq2(ns2)*azero(ns2 - nnn)
     $          + conc(ns3)*zchsq2(ns3)*azero(ns3 - nnn)
     $          + conc(ns4)*zchsq2(ns4)*azero(ns4 - nnn)
     $          + conc(ns5)*zchsq2(ns5)*azero(ns5 - nnn)
     $          + conc(ns6)*zchsq2(ns6)*azero(ns6 - nnn)
     $          + conc(ns7)*zchsq2(ns7)*azero(ns7 - nnn)
      enddo
c
      do nss = ileft + 1,narn2
        ns = jcsort(nss)
        sx = sx + conc(ns)*zchsq2(ns)*azero(ns - nnn)
      enddo
c
      abar = sx/fxi
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
