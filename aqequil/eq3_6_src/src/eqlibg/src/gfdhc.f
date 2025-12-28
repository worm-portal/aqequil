      subroutine gfdhc(adh,bt,fdhc,fdhcp,fdhcpp,fxi)
c
c     This subroutine computes the Debye-Huckel function f and its
c     derivatives with respect to ionic strength for the case of the
c     Debye-Huckel-charging (DHC) model.
c
c     This subroutine is called by:
c
c       EQLIBG/gcoeff.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       adh    = the Debye-Huckel A(gamma) parameter
c       bt     = product of Debye-Huckel B and the hard core diameter
c       fxi    = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
c
c     Output:
c
c       fdhc   = Debye-Huckel function f (DHC model)
c       fdhcp  = f' (df/dI)
c       fdhcpp = f'' (d2f/dI2)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      real*8 adh,bt,fdhc,fdhcp,fdhcpp,fxi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      real*8 fxisqt,vx,vxp1
c
c-----------------------------------------------------------------------
c
      fxisqt = sqrt(fxi)
      vx = bt*fxisqt
      vxp1 = vx + 1.
c
      fdhc = -6.*adh*(vx*(vx-2.) + 2.*log(vxp1))/(bt**3)
      fdhcp = -6.*adh*fxisqt/vxp1
      fdhcpp = -3.*adh/(fxisqt*vxp1*vxp1)
c
      end
