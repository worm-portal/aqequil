      subroutine betacf(acflg,acflgo,bacfmx,nst,nstmax,ubacmx,uspec)
c
c     This subroutine finds the activity coefficient residual with
c     the largest magnitude (bacfmx). The range of activity coefficient
c     residuals covers all species.
c
c     This subroutine is called by:
c
c       EQLIB/ngcadv.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       uspec  = name array for aqueous species
c       acflg  = array of current values of the log activity
c                  coefficients
c       acflgo = array of the previous values of the log activity
c                  coefficients
c
c     Principal output:
c
c       bacfmx = the largest activity coefficient residual
c       ubacmx = the name of the species with the largest activity
c                  coefficient residual is bacfmx
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
      integer nst
c
      character*48 uspec(nstmax),ubacmx
c
      real*8 acflg(nstmax),acflgo(nstmax)
      real*8 bacfmx
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ileft,ns
c
      character*48 ux48
c
      real*8 adacfl,bx
c
c-----------------------------------------------------------------------
c
c     Note the use of local variables (bx and ux48) within the loop.
c     Note also that the loop is unrolled.
c
      bacfmx = 0.
      ubacmx = 'None'
      bx = 0.
      ux48 = 'None'
      ileft = (nst/8)*8
c
      do ns = 1,ileft,8
        adacfl = abs(acflg(ns) - acflgo(ns))
        if (adacfl .gt. bx) then
          bx = adacfl
          ux48 = uspec(ns)
        endif
        adacfl = abs(acflg(ns + 1) - acflgo(ns + 1))
        if (adacfl .gt. bx) then
          bx = adacfl
          ux48 = uspec(ns + 1)
        endif
        adacfl = abs(acflg(ns + 2) - acflgo(ns + 2))
        if (adacfl .gt. bx) then
          bx = adacfl
          ux48 = uspec(ns + 2)
        endif
        adacfl = abs(acflg(ns + 3) - acflgo(ns + 3))
        if (adacfl .gt. bx) then
          bx = adacfl
          ux48 = uspec(ns + 3)
        endif
        adacfl = abs(acflg(ns + 4) - acflgo(ns + 4))
        if (adacfl .gt. bx) then
          bx = adacfl
          ux48 = uspec(ns + 4)
        endif
        adacfl = abs(acflg(ns + 5) - acflgo(ns + 5))
        if (adacfl .gt. bx) then
          bx = adacfl
          ux48 = uspec(ns + 5)
        endif
        adacfl = abs(acflg(ns + 6) - acflgo(ns + 6))
        if (adacfl .gt. bx) then
          bx = adacfl
          ux48 = uspec(ns + 6)
        endif
        adacfl = abs(acflg(ns + 7) - acflgo(ns + 7))
        if (adacfl .gt. bx) then
          bx = adacfl
          ux48 = uspec(ns + 7)
        endif
      enddo
c
      do ns = ileft + 1,nst
        adacfl = abs(acflg(ns) - acflgo(ns))
        if (adacfl .gt. bx) then
          bx = adacfl
          ux48 = uspec(ns)
        endif
      enddo
c
      bacfmx = bx
      ubacmx = ux48
c
      end
