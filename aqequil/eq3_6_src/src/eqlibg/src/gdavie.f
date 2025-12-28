      subroutine gdavie(acflgc,actwlc,adh,al10,fxi,narn1,narn2,nstmax,
     $ omega,sigmam,xbrwlc,zchsq2)
c
c     This subroutine computes activity coefficients of aqueous species
c     using the Davies (1961) equation. The activity of water is
c     computed from an expression that was derived from the Davies
c     equation using thermodynamic consistency relations.
c
c     This subroutine is called by:
c
c        EQLIBG/gcoeff.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       adh    = Debye-Huckel A(gamma) parameter
c       omega  = water constant; ~55.51.
c       sigmam = sum of solute molalities
c       xbrwlc = log mole fraction of water
c       fxi    = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
c       zchsq2 = one-half the charge squared array
c
c     Principal output:
c
c       acflgc = array of log activity coefficients
c       actwlc = log activity of water
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
      integer narn1,narn2
c
      real*8 acflgc(nstmax),zchsq2(nstmax)
      real*8 actwlc,adh,al10,fxi,omega,sigmam,xbrwlc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ns
c
      real*8 factor,fxisqt,sga,sgx,xxp1
c
c-----------------------------------------------------------------------
c
      fxisqt = sqrt(fxi)
      xxp1 = 1. + fxisqt
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute log lambda(w).
c
      sgx = (3./(fxisqt**3))*( xxp1 - (1./xxp1) - 2.*log(xxp1) )
      sga = sigmam/al10
      actwlc = ( -sga + 2.*adh*(fxi**1.5)*sgx/3.
     $ - 0.2*adh*fxi*fxi )/omega
      acflgc(narn1) = actwlc - xbrwlc
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute log gamma(i).
c
      factor = (fxisqt/xxp1) - 0.2*fxi
      factor = -2*adh*factor
      do ns = narn1 + 1,narn2
        acflgc(ns) = zchsq2(ns)*factor
      enddo
c
      end
