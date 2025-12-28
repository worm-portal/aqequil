      subroutine gszm(conc,jcsort,narn1,narn2,nstmax,sigza,sigzc,
     $ sigzi,sigzm,zchar)
c
c     This subroutine calculates the sums of equivalent concentrations
c     and the charge imbalance. Note that a sorted summation is
c     used.
c
c     This subroutine is called by:
c
c       EQLIB/betas.f
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
c       narn1  = start of the range of aqeuous species; also
c                  the species index of solvent water
c       narn2  = end of the range of aqeuous species
c       zchar  = array of electrical charge numbers
c
c     Principal output:
c
c       sigza  = the sum of equivalent concentrations of aqueous
c                  anions, SUM(i) abs(z(i))m(i), for z(i) < 0
c       sigzc  = the sum of equivalent concentrations of aqueous
c                  cations, SUM(i) z(i)m(i), for z(i) > 0
c       sigzi  = the calculated charge imbalance, sigzc - sigza
c       sigzm  = the sum of equivalent concentrations of aqueous
c                  ions, sigzc + sigza
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
      real*8 conc(nstmax),zchar(nstmax)
      real*8 sigza,sigzc,sigzi,sigzm
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ileft,nnn,ns,nss,ns0,ns1,ns2,ns3,ns4,ns5,ns6,ns7,nval
c
      real*8 ec,sxa,sxc,sxi,sxm
c
c-----------------------------------------------------------------------
c
c     Logically, this could be represented by:
c
c       sigzc = 0.
c       sigza = 0.
c       do ns = narn1,narn2
c         if (zchar(ns) .gt. 0.) then
c           sigzc = sigzc + conc(ns)*zchar(ns)
c         elseif (zchar(ns) .lt. 0.) then
c           sigza = sigza + abs(conc(ns)*zchar(ns))
c         endif
c       enddo
c       sigzi = sigzc - sigza
c       sigzm = sigzc + sigza
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      sigzc = 0.
      sigza = 0.
      sigzm = 0.
      sigzi = 0.
c
c     Note the use of local variables (sxc, sxa, sxi, and sxm) within
c     the loop. Note also that the loop is unrolled.
c
      sxc = 0.
      sxa = 0.
      sxi = 0.
      sxm = 0.
      nnn = narn1 - 1
      nval = narn2 - nnn
      ileft = (nval/8)*8 + nnn
c
      do nss = narn1,ileft,8
        ns0 = jcsort(nss)
        if (zchar(ns0) .gt. 0.) then
          ec = zchar(ns0)*conc(ns0)
          sxc = sxc + ec
          sxi = sxi + ec
          sxm = sxm + ec
        elseif (zchar(ns0) .lt. 0.) then
          ec = - zchar(ns0)*conc(ns0)
          sxa = sxa + ec
          sxi = sxi - ec
          sxm = sxm + ec
        endif
        ns1 = jcsort(nss + 1)
        if (zchar(ns1) .gt. 0.) then
          ec = zchar(ns1)*conc(ns1)
          sxc = sxc + ec
          sxi = sxi + ec
          sxm = sxm + ec
        elseif (zchar(ns1) .lt. 0.) then
          ec = - zchar(ns1)*conc(ns1)
          sxa = sxa + ec
          sxi = sxi - ec
          sxm = sxm + ec
        endif
        ns2 = jcsort(nss + 2)
        if (zchar(ns2) .gt. 0.) then
          ec = zchar(ns2)*conc(ns2)
          sxc = sxc + ec
          sxi = sxi + ec
          sxm = sxm + ec
        elseif (zchar(ns2) .lt. 0.) then
          ec = - zchar(ns2)*conc(ns2)
          sxa = sxa + ec
          sxi = sxi - ec
          sxm = sxm + ec
        endif
        ns3 = jcsort(nss + 3)
        if (zchar(ns3) .gt. 0.) then
          ec = zchar(ns3)*conc(ns3)
          sxc = sxc + ec
          sxi = sxi + ec
          sxm = sxm + ec
        elseif (zchar(ns3) .lt. 0.) then
          ec = - zchar(ns3)*conc(ns3)
          sxa = sxa + ec
          sxi = sxi - ec
          sxm = sxm + ec
        endif
        ns4 = jcsort(nss + 4)
        if (zchar(ns4) .gt. 0.) then
          ec = zchar(ns4)*conc(ns4)
          sxc = sxc + ec
          sxi = sxi + ec
          sxm = sxm + ec
        elseif (zchar(ns4) .lt. 0.) then
          ec = - zchar(ns4)*conc(ns4)
          sxa = sxa + ec
          sxi = sxi - ec
          sxm = sxm + ec
        endif
        ns5 = jcsort(nss + 5)
        if (zchar(ns5) .gt. 0.) then
          ec = zchar(ns5)*conc(ns5)
          sxc = sxc + ec
          sxi = sxi + ec
          sxm = sxm + ec
        elseif (zchar(ns5) .lt. 0.) then
          ec = - zchar(ns5)*conc(ns5)
          sxa = sxa + ec
          sxi = sxi - ec
          sxm = sxm + ec
        endif
        ns6 = jcsort(nss + 6)
        if (zchar(ns6) .gt. 0.) then
          ec = zchar(ns6)*conc(ns6)
          sxc = sxc + ec
          sxi = sxi + ec
          sxm = sxm + ec
        elseif (zchar(ns6) .lt. 0.) then
          ec = - zchar(ns6)*conc(ns6)
          sxa = sxa + ec
          sxi = sxi - ec
          sxm = sxm + ec
        endif
        ns7 = jcsort(nss + 7)
        if (zchar(ns7) .gt. 0.) then
          ec = zchar(ns7)*conc(ns7)
          sxc = sxc + ec
          sxi = sxi + ec
          sxm = sxm + ec
        elseif (zchar(ns7) .lt. 0.) then
          ec = - zchar(ns7)*conc(ns7)
          sxa = sxa + ec
          sxi = sxi - ec
          sxm = sxm + ec
        endif
      enddo
c
      do nss = ileft + 1,narn2
        ns = jcsort(nss)
        if (zchar(ns) .gt. 0.) then
          ec = zchar(ns)*conc(ns)
          sxc = sxc + ec
          sxi = sxi + ec
          sxm = sxm + ec
        elseif (zchar(ns) .lt. 0.) then
          ec = - zchar(ns)*conc(ns)
          sxa = sxa + ec
          sxi = sxi - ec
          sxm = sxm + ec
        endif
      enddo
c
      sigzc = sxc
      sigza = sxa
      sigzi = sxi
      sigzm = sxm
c
      end
