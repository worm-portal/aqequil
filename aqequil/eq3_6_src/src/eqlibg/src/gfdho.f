      subroutine gfdho(aphi,bt,f,fp,fpp,fxi)
c
c     This subroutine computes the Debye-Huckel function f and its
c     derivatives with respect to ionic strength for the case of
c     the Debye-Huckel-osmotic (DHO) model. This is the Debye-Huckel
c     model in Pitzer's (1973, 1975) equations.
c
c     This subroutine is called by:
c
c       EQLIBG/gcoeff.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       aphi   = the Debye-Huckel A(phi) parameter
c       bt     = product of Debye-Huckel B and the hard core diameter,
c                usually taken as a constant (1.2)
c       fxi    = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
c
c     Output:
c
c       f      = Debye-Huckel function f
c       fp     = f' (df/dI)
c       fpp    = f'' (d2f/dI2)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      real*8 aphi,bt,f,fp,fpp,fxi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      real*8 alvxp1,fc,fxisqt,vx,vxp1,vxp1sq
c
c-----------------------------------------------------------------------
c
      fxisqt = sqrt(fxi)
      vx = bt*fxisqt
      vxp1 = vx + 1.
      vxp1sq = vxp1*vxp1
      alvxp1 = log(vxp1)
c
      fc = -4.*aphi/bt
      f = fc*fxi*alvxp1
      fp = fc*((vx / (2.*vxp1)) + alvxp1)
      fpp = (-aphi/fxisqt)*((1. + 2.*vxp1) / vxp1sq)
c
      end
