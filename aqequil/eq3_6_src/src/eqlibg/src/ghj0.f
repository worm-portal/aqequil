      subroutine ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
c
c     This subroutine evaluates the function J0(x) (hj0), its first
c     two derivatives with respect to x (dhj0 and d2hj0), and
c     the related functions J1(x) (hj1) and J2(x) (hj2). The
c     algorithm used here is based on Chebyshev polynomials and is
c     taken from the Appendix B of Harvie (1981). The treatment of
c     the second derivative and of J2(x) is an extension of Harvie's
c     work. The Harvie (1981) approximation is considered to be
c     more accurate than the earlier approximation given by Pitzer
c     (1975) and is used in almost all modern work involving
c     Pitzer's equations.
c
c     These functions are used to compute E-lambda and E-theta
c     coefficients and their derivatives with respect to ionic
c     strength. They describe the higher-order electrostatic effects
c     in Pitzer's equations.
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
c       dhj0   = the derivative dJ0(x)/dx
c       d2hj0  = the derivative d2 J0(x)/dx2
c       hj0    = the function J0(x)
c       hj1    = the function J1(x)
c       hj2    = the function J2(x)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      real*8 dhj0,d2hj0,hj0,hj1,hj2,x
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer k,ki
c
      real*8 ai(0:20),aii(0:20),b(0:22),d(0:22),e(0:22)
      real*8 bdif,ddif,edif,dz,d2z,z
c
c-----------------------------------------------------------------------
c
      data b(21) /0./,b(22) /0./
      data d(21) /0./,d(22) /0./
      data e(21) /0./,e(22) /0./
c
      data (ai(k),k = 0,9)
     $         /1.92515 40148 14667,
     $          -.06007 64777 53119,
     $          -.02977 90774 56514,
     $          -.00729 94996 90937,
     $           .00038 82606 36404,
     $           .00063 68745 99598,
     $           .00003 65836 01823,
     $          -.00004 50369 75204,
     $          -.00000 45378 95710,
     $           .00000 29377 06971/
c
      data (ai(k), k = 10,20)
     $          /.00000 03965 66462,
     $          -.00000 02020 99617,
     $          -.00000 00252 67769,
     $           .00000 00135 22610,
     $           .00000 00012 29405,
     $          -.00000 00008 21969,
     $          -.00000 00000 50847,
     $           .00000 00000 46333,
     $           .00000 00000 01943,
     $          -.00000 00000 02563,
     $          -.00000 00000 10991/
c
      data (aii(k), k = 0,9)
     $         / .62802 33205 20852,
     $           .46276 29853 38493,
     $           .15004 46371 87895,
     $          -.02879 60576 04906,
     $          -.03655 27459 10311,
     $          -.00166 80879 45272,
     $           .00651 98403 98744,
     $           .00113 03780 79086,
     $          -.00088 71713 10131,
     $          -.00024 21076 41309/
c
      data (aii(k), k = 10,20)
     $          /.00008 72944 51594,
     $           .00003 46821 22751,
     $          -.00000 45837 68938,
     $          -.00000 35486 84306,
     $          -.00000 02504 53880,
     $           .00000 02169 91779,
     $           .00000 00807 79570,
     $           .00000 00045 58555,
     $          -.00000 00069 44757,
     $          -.00000 00028 49257,
     $           .00000 00002 37816/
c
c-----------------------------------------------------------------------
c
      hj0 = 0.
      dhj0 = 0.
      d2hj0 = 0.
c
      if (x .le. 0.) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (x .le. 1.) then
c
c       Case 1.
c
        z = 4.*(x**0.2) - 2.
        dz = 0.8*(x**(-0.8))
        d2z = 0.64*(x**(-1.8))
c
        do ki = 0,20
          k = 20 - ki
          b(k) = z*b(k + 1) - b(k + 2) + ai(k)
          d(k) = b(k + 1) + z*d(k + 1) - d(k + 2)
          e(k) = 2.*d(k + 1) + z*e(k + 1) - e(k + 2)
      enddo
c
      else
c
c       Case 2.
c
        z = 4.44444444444444*(x**(-0.1)) - 2.44444444444444
        dz = -0.444444444444444*(x**(-1.1))
        d2z = 0.488888888888889*(x**(-2.1))
c
        do ki = 0,20
          k = 20 - ki
          b(k) = z*b(k + 1) - b(k + 2) + aii(k)
          d(k) = b(k + 1) + z*d(k + 1) - d(k + 2)
          e(k) = 2.*d(k + 1) + z*e(k + 1) - e(k + 2)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute (b0 - b2), (d0 - d2), and (e0 - e2).
c
      bdif = b(0) - b(2)
      ddif = d(0) - d(2)
      edif = e(0) - e(2)
c
      hj0 = 0.25*x - 1.0 + bdif/2.
      dhj0 = 0.25 + 0.5*dz*ddif
      d2hj0 = 0.5*(edif*(dz**2) + ddif*d2z)
c
      hj1 = x*dhj0
      hj2 = x*d2hj0
c
  999 continue
      end
