      subroutine ca3bar(azero,a3bar,a3bars,conc,jcsort,narn1,
     $ narn2,natmax,nstmax,sigmam)
c
c     This subroutine calculates the characteristic average cubed
c     distance of closest approach for each aqueous solute species
c     (a3bars) and the average cubed distance of closest approach for
c     all solute species (a3bar). These quantities appear in the
c     first order term representing hard core repulsion in models
c     of aqueous species activity coefficients. In the calculation
c     of these quantities, the weighting factor corresponding to a
c     species is its molality.
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
c       sigmam = sum of the molalities of the aqueous solute species
c
c     Prinicpal output:
c
c       a3bars = characteristic average cubed distance of closest
c                  approach array
c       a3bar  = average cubed distance of closest approach
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
      integer narn1,narn2
c
      real*8 azero(natmax),a3bars(natmax),conc(nstmax)
      real*8 a3bar,sigmam
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ileft,na,nnn,ns,nsj,nss,ns0,ns1,ns2,ns3,ns4,ns5,ns6,
     $ ns7,nval
c
      real*8 aij,aij0,aij1,aij2,aij3,aij4,aij5,aij6,aij7,azi,cw,sx
c
c-----------------------------------------------------------------------
c
      do na = 1,narn2 - narn1 + 1
        a3bars(na) = 0
      enddo
      a3bar = 0
      if (sigmam .le. 0.) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the average value of the cube of the distance of closest
c     approach for each aqueous solute species. Logically, this could be
c     represented by:
c
c       do ns = narn1 + 1,narn2
c         na = ns - narn1 + 1
c         sx = 0.
c         do nsj = narn1 + 1,narn2
c           naj = nsj - narn1 + 1
c           aij = 0.5*(azero(na) + azero(naj))
c           sx = sx + conc(nsj)*(aij**3)
c         enddo
c         a3bars(na) = sx/sigmam
c       enddo
c
c     In doing a sorted summation in the inner loop, a complication
c     arises in that water (ns = narn1) must be included in the
c     range of the loop, because the jcsort array is organized to
c     contain the species indices sorted according to increasing
c     concentration within phase ranges. There is no available
c     sorting among aqueous solute species only. This complication
c     is dealt with by temporarily setting the concentration of water
c     (conc(narn1)) to zero. Note the use of a local variable (sx)
c     within the loop. Note also that the loop is unrolled.
c
      cw = conc(narn1)
      conc(narn1) = 0.
      nnn = narn1 - 1
      nval = narn2 - nnn
      ileft = (nval/8)*8 + nnn
c
      do ns = narn1 + 1,narn2
        na = ns - narn1 + 1
        azi = azero(na)
        sx = 0.
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
          aij0 = 0.5*(azi + azero(ns0 - nnn))
          aij1 = 0.5*(azi + azero(ns1 - nnn))
          aij2 = 0.5*(azi + azero(ns2 - nnn))
          aij3 = 0.5*(azi + azero(ns3 - nnn))
          aij4 = 0.5*(azi + azero(ns4 - nnn))
          aij5 = 0.5*(azi + azero(ns5 - nnn))
          aij6 = 0.5*(azi + azero(ns6 - nnn))
          aij7 = 0.5*(azi + azero(ns7 - nnn))
          sx = sx + conc(ns0)*(aij0**3) + conc(ns1)*(aij1**3)
     $            + conc(ns2)*(aij2**3) + conc(ns3)*(aij3**3)
     $            + conc(ns4)*(aij4**3) + conc(ns5)*(aij5**3)
     $            + conc(ns6)*(aij6**3) + conc(ns7)*(aij7**3)
        enddo
c
        do nss = ileft + 1,narn2
          nsj = jcsort(nss)
          aij = 0.5*(azi + azero(nsj - nnn))
          sx = sx + conc(ns)*(aij**3)
        enddo
c
        a3bars(na) = sx/sigmam
      enddo
c
      conc(narn1) = cw
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the average value of the cube of the distance of closest
c     approach for all aqeuous solute species. Logically, this could be
c     represented by:
c
c       sx = 0.
c       do ns = narn1 + 1,narn2
c         na = ns - narn1 + 1
c         sx = sx + conc(ns)*a3bars(na)
c       enddo
c       a3bar = sx/sigmam
c
c     Because a3bars(narn1) is defined to be zero, a sorted summation
c     can be safely made by looping from narn1 to narn2. Note the use
c     of a local variable (sx) within the loop. Note also that the loop
c     is unrolled.
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
        sx = sx + conc(ns0)*a3bars(ns0 - nnn)
     $          + conc(ns1)*a3bars(ns1 - nnn)
     $          + conc(ns2)*a3bars(ns2 - nnn)
     $          + conc(ns3)*a3bars(ns3 - nnn)
     $          + conc(ns4)*a3bars(ns4 - nnn)
     $          + conc(ns5)*a3bars(ns5 - nnn)
     $          + conc(ns6)*a3bars(ns6 - nnn)
     $          + conc(ns7)*a3bars(ns7 - nnn)
      enddo
c
      do nss = ileft + 1,narn2
        ns = jcsort(nss)
        sx = sx + conc(ns)*a3bars(ns - nnn)
      enddo
c
      a3bar = sx/sigmam
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
