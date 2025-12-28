      subroutine ctds(jcsort,mosp,mwtsp,narn1,narn2,nstmax,wotds)
c
c     This subroutine calculates the total dissolved solute mass
c     (wotds, g). Note that a sorted summation is used.
c
c     This subroutine is called by:
c
c       EQLIB/gwdenp.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       mosp   = array of numbers of moles of species
c       jcsort = array of species indices, in order of increasing
c                  concentration/number of moles, but with sorting
c                  restricted to within phase ranges
c       mwtsp  = array of molecular weights of species
c       narn1  = start of the range of aqueous species
c       narn2  = end of the range of aqueous species
c
c     Principal output:
c
c       wotds  = the total dissolved solute mass (g)
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
      real*8 mosp(nstmax),mwtsp(nstmax)
      real*8 wotds
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ileft,ns,nss,ns0,ns1,ns2,ns3,ns4,ns5,ns6,ns7,nval
c
      real*8 mxw,wx
c
c-----------------------------------------------------------------------
c
c     Logically, this could be represented by:
c
c       wotds = 0.
c       do ns = narn1 + 1,narn2
c         wotds = wotds + mosp(ns)*mwtsp(ns)
c       enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     In doing a sorted summation, a complication arises in that
c     water (ns = narn1) must be included in the range of the
c     loop, because the jcsort array is organized to contain the
c     species indices sorted according to increasing concentration
c     within phase ranges. There is no available sorting among
c     aqueous solute species only. Thus, the concentration of water
c     (conc(nanrn1)) is temporarily set to zero. Note the use
c     of a local variable (wx) within the loop. Note also that the
c     loop is unrolled.
c
      wx = 0.
      mxw = mosp(narn1)
      mosp(narn1) = 0.
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
        wx = wx + mosp(ns0)*mwtsp(ns0) + mosp(ns1)*mwtsp(ns1)
     $          + mosp(ns2)*mwtsp(ns2) + mosp(ns3)*mwtsp(ns3)
     $          + mosp(ns4)*mwtsp(ns4) + mosp(ns5)*mwtsp(ns5)
     $          + mosp(ns6)*mwtsp(ns6) + mosp(ns7)*mwtsp(ns7)
      enddo
c
      do nss = ileft + 1,narn2
        ns = jcsort(nss)
        wx = wx + mosp(ns)*mwtsp(ns)
      enddo
c
      wotds = wx
      mosp(narn1) = mxw
c
      end
