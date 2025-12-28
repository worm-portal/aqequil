      subroutine fdomsp(jssort,mosp,nsi,nst,nstmax,weight,wsi)
c
c     This subroutine finds the species that dominates a mass balance.
c     The primary purpose of this subroutine is to support the
c     pre-Newton-Raphson optimization algorithm. This subroutine is not
c     intended to support optimizing the basis set.
c
c     Note - to save time, it is assumed that the ratio of the
c     largest value of weight to the smallest non-zero value is
c     no greater than 100.
c
c     This subroutine is called by:
c
c       EQLIB/cfracf.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       jssort = array of species indices in order of increasing
c                  mass
c       weight = array of stoichiometric weighting factors
c       mosp   = array of moles of species
c       nst    = number of species
c
c     Principal output:
c
c       nsi    = index of the dominant species
c       wsi    = weight(nsi)
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
      integer jssort(nstmax)
      integer nsi,nst
c
      real*8 mosp(nstmax),weight(nstmax)
      real*8 wsi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nsc,nss,ns
c
      real*8 ap0,ap1,m0,m1,w0,w1,p1,rx
c
c-----------------------------------------------------------------------
c
      m0 = 0.
      w0 = 1.
      ap0 = 0.
      nsi = 0
c
      do nsc = 1,nst
        nss = nst + 1 - nsc
        ns = jssort(nss)
        w1 = weight(ns)
        m1 = mosp(ns)
        if (m1 .le. 0.) go to 20
        if (w1 .ne. 0.) then
          p1 = w1*m1
          ap1 = abs(p1)
          if (ap1 .gt. ap0) then
            m0 = m1
            w0 = w1
            ap0 = ap1
            nsi = ns
            go to 15
          endif
        endif
        rx = m0/m1
        if (rx .gt. 100.) go to 20
   15   continue
      enddo
c
   20 continue
      wsi = w0
      end
