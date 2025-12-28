      subroutine gpj0(dpj0,d2pj0,pj0,x)
c
c     This subroutine evaluates the function J0(x) (pj0) and its
c     first and second derivatives with respect to x (dpj0 and
c     d2pj0). The evaluation here uses an approximation given by
c     Pitzer (1975, eq. 47). Most modern usage of Pitzer's equations
c     is consistent with the later approximation of Harvie (1981).
c
c     These functions are used to compute E-lambda and E-theta
c     coefficients and their derivatives with respect to ionic
c     strength. They describe the higher-order electrostatic effects
c     in Pitzer's equations.
c
c     It is important to note that the tabulation of J(x) and J'(x)
c     in Table II of Pitzer (1975) is not for the approximation
c     represented by eq. 47. Pitzer (1975) does not provide such a
c     table for the eq. 47 approximation.
c
c                            References
c
c     Harvie, C. E. 1981. Theoretical Investigations in Geochemistry
c       and Atom Surface Scattering. Ph.D. dissertation, University
c       of California, San Diego (Available as #8203026 from
c       University Microfilms International, Ann Arbor, Michigan).
c
c     Pitzer, K.S. 1975. Thermodynamics of electrolytes. V. Effects
c       of higher-order electrostatic terms. Journal of Solution
c       Chemistry, v. 4, p. 249-265.
c
c     This subroutine is called by:
c
c       EQLIBG/elmdd.f
c       EQLIBG/cwrpjt.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       x      = independent variable
c
c     Output:
c
c       dpj0   = the derivative dJ0(x)/dx
c       d2pj0  = the derivative d2 J0(x)/dx2
c       pj0    = the function J0(x)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      real(8) dpj0,d2pj0,pj0,x
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      real(8) alx,bk2,bk3,c1,c2,c3,c4,xc2,c3xc4,ec3xc4,t,xc2p
c
c-----------------------------------------------------------------------
c
c     These parameters are those for eq. 47 given in Table III of
c     Pitzer (1975).
c
      data c1,c2,c3,c4 /4.581,0.7237,0.0120,0.528/
c
c-----------------------------------------------------------------------
c
      pj0 = 0.
      dpj0 = 0.
      d2pj0 = 0.
c
      if (x .le. 0.) go to 999
c
      alx = log(x)
c
      xc2 = exp(-c2 * alx)
      c3xc4 = c3*exp(c4 * alx)
      ec3xc4 = exp(-c3xc4)
c
      pj0 = x / (4. + c1*xc2*ec3xc4)
c
      t = c4*c3xc4
      xc2p = exp((-c2 - 1.)*alx)
      bk2 = c1*xc2p *(c2 + t)*ec3xc4
c
      dpj0 = (pj0 / x)*(1. + pj0*bk2)
c
      bk3 = (t/x)*((c4 / (c2 + t)) - (c2 + 1.) / t)
      d2pj0 = ((dpj0 / pj0) - 1. / x)*(2.*dpj0 + pj0*bk3)
c
  999 continue
      end
