      subroutine elmdd(aphi,el,elp,elpp,fxi,ijz,qpit75)
c
c     This subroutine calculates the E-lambda function (el) and its
c     first two derivatives (elp and elpp) with respect to ionic
c     strength (fxi) for the charge pair product ijz. The E-lambda
c     function and its derivatives are used to evaluate Pitzer's
c     equations. Here aphi is the Debye-Huckel A(phi) parameter.
c
c     Note: el = elp = elpp = 0 if ijz is less than or equal to 0.
c
c     This subroutine is called by:
c
c       EQLIBG/gelam.f.
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       aphi   = Debye-Huckel A(phi) parameter
c       fxi    = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
c       ijz    = input charge product
c
c     Principal output:
c
c       el     = array of values of E-lambda functions
c       elp    = array of first ionic strength derivatives of
c                  E-lambda functions
c       elpp   = array of second ionic strength derivatives of
c                  E-lambda functions
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ijz
c
      logical qpit75
c
      real(8) aphi,el,elp,elpp,fxi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      real(8) ck,cx,dhj0,d2hj0,hj0,hj1,hj2,ri,sjp,sjpp,x,fxisqt,xp,xpp
c
c-----------------------------------------------------------------------
c
      el = 0.
      elp = 0.
      elpp = 0.
c
      if (ijz .gt. 0.) then
c
        fxisqt = sqrt(fxi)
        cx = 3.*ijz*aphi
        x = 2.*cx*fxisqt
c
c       Compute x' and x''.
c
        xp = cx/fxisqt
        xpp = -(0.5*xp)/fxi
c
c       Get J0(x) and related functions.
c
        if (.not.qpit75) then
c
c         Use the Chebyshev polynomial approximation
c         of Harvie (1981). This is what should normally
c         be used. Note that hj1 and hj2 are outputs
c         that are not used.
c
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
c
c         Use the less accurate approximation of
c         Pitzer (1975).
c
c         Calling sequence substitutions:
c           dhj0 for dpj0
c           d2hj0 for d2pj0
c           hj0 for pj0
c
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
c
c       Convert these to derivatives with respect to I.
c
        sjp = dhj0*xp
        sjpp = dhj0*xpp + d2hj0*xp*xp
c
        ck = 0.25*ijz
        ri = 1./fxi
c
        el = ri*ck*hj0
        elp = ri*(ck*sjp - el)
        elpp = ri*(ck*sjpp - 2.*elp)
      endif
c
      end
