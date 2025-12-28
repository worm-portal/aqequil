      subroutine gbdot(acflgc,actwlc,adh,al10,azero,bdh,bdot,cco2,
     $ fxi,insgf,narn1,narn2,natmax,nstmax,omega,sigmam,tempk,
     $ xbrwlc,zchar,zchsq2)
c
c     This subroutine computes activity coefficients of aqueous species
c     using the B-dot equation and a set of related approximations.
c     This model is suitable only for relatively dilute solutions
c     (ionic strength no greater than one molal, though demonstrable
c     inaccuracy appears at much lower values).
c
c     The B-dot equation is an extended Debye-Huckel equation (Helgeson,
c     1969). It is applied only to electrically charged species. The
c     activity coefficients of electrically neutral solute species are
c     treated according to the value of the insgf flag. If it is -1,
c     the activity coefficient is computed from the Drummond (1981)
c     equation. This is a fit of the log activity coefficient of
c     CO2(aq) as a function of temperature and ionic strength (ignoring
c     ion pairing) in pure aqueous sodium chloride solutions. This is
c     an appropriate treatment for neutral species that are nonpolar.
c     If insgf is 0, the log activity coefficient is set to zero.
c     The activity of water is computed from an equation which was
c     derived from the B-dot equation using thermodynamic consistency
c     relations, though it is actually consistent with that equation
c     only for the case of all ions having the same size and an
c     electrical charge of +1 or -1.
c
c     This subroutine is called by:
c
c       EQLIBG/gcoeff.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       adh    = Debye-Huckel A(gamma) parameter
c       azero  = hard core diameter array
c       bdh    = Debye-Huckel B parameter
c       bdot   = B-dot parameter
c       cco2   = coefficients of the Drummond (1981) equation
c       insgf  = flag array for electrically neutral solute species
c       omega  = water constant; ~55.51.
c       sigmam = sum of solute molalities
c       xbrwlc = log mole fraction of water
c       fxi    = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
c       zchar  = electrical charge array
c       zchsq2 = one-half the charge squared array
c
c      Principal output:
c
c       acflgc = log activity coefficient array
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
      integer natmax,nstmax
c
      integer insgf(natmax)
      integer narn1,narn2
c
      real*8 acflgc(nstmax),azero(natmax),cco2(5),zchar(nstmax),
     $ zchsq2(nstmax)
      real*8 actwlc,adh,al10,bdh,bdot,fxi,omega,sigmam,tempk,xbrwlc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer na,ns
c
      real*8 acfco2,art,bdtfxi,brt,fxisqt,sga,sgx,xx,xxp1
c
c-----------------------------------------------------------------------
c
      fxisqt = sqrt(fxi)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute log lambda(w).
c
      xx = 4.0*bdh*fxisqt
      xxp1 = xx + 1.
      sgx = (3./(xx**3)) * ( xxp1 - (1./xxp1) - 2.*log(xxp1) )
      sga = sigmam/al10
      actwlc = ( -sga + 2.*adh*(fxi**1.5)*sgx/3. - bdot*fxi*fxi )/omega
      acflgc(narn1) = actwlc - xbrwlc
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute log gamma(i).
c
c     Compute the log activity coefficient of CO2(aq), using the
c     equation proposed by Drummond (1981, thesis, Penn State
c     University).
c
      acfco2 = ( cco2(1) + cco2(2)*tempk + cco2(3)/tempk )*fxi
     $ - ( cco2(4) + cco2(5)*tempk )*(fxi/(fxi + 1.))
      acfco2 = acfco2/al10
c
      art = 2*adh*fxisqt
      brt = bdh*fxisqt
      bdtfxi = bdot*fxi
c
      do ns = narn1 + 1,narn2
        na = ns - narn1 + 1
        if (zchar(ns) .ne. 0.) then
c
c         Have an ion. Use the B-dot equation.
c
          acflgc(ns)  = -( (art*zchsq2(ns)) / (1.0 + brt*azero(na)) )
     $    + bdtfxi
        else
c
c         Have a neutral species.
c
          if (insgf(na) .le. -1) then
c
c           For nonpolar species, using the result for CO2(aq) from
c           the Drummond (1981) equation.
c
            acflgc(ns) = acfco2
          else
c
c           For polar species, assign a log gamma value of zero.
c
            acflgc(ns) = 0.
          endif
        endif
      enddo
c
      end
