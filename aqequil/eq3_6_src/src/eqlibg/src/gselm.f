      subroutine gselm(conc,delam,dselm,elam,izmax,narn1,narn2,nazmmx,
     $ nazpmx,nstmax,selm,zchar)
c
c     This subroutine computes the following first order sums used
c     in Pitzer's equations:
c
c       selm(i):  SUM(j) E-lambda(ij)*m(j)
c       dselm(1,i): SUM(j) E-lambda'(ij)*m(j)
c
c     This subroutine is called by:
c
c       EQLIBG/gcoeff.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       conc   = array of concentrations
c       elam   = array of E-lambda functions
c       delam  = array of ionic strength derivatives of E-lambda
c                  functions
c       narn1  = start of species range for aqueous solution
c       narn2  = end of species range for aqueous solution
c       zchar  = array of charges
c
c     Principal output:
c
c       selm   = array of sums: 2 SUM(j) E-lambda(ij)*m(j)
c       dselm  = array of the corresponding derivatives with respect
c                  to ionic strength
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nazmmx,nazpmx,nstmax
c
      integer izmax,narn1,narn2
c
      real*8 delam(2,nazpmx,nazpmx),dselm(2,nazmmx:nazpmx),
     $ conc(nstmax),elam(nazpmx,nazpmx),selm(nazmmx:nazpmx),
     $ zchar(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer iajz,iz,j,jz
c
      real*8 sum1,sum1p,sum2,sum2p,zj
c
c-----------------------------------------------------------------------
c
c     Note: nazmmx = -nazpmx.
c
      do iz = nazmmx,nazpmx
        selm(iz) = 0.
        dselm(1,iz) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Looping over iz from 1 to izmax, get the quantities for charges
c     iz and -iz simultaneously.
c
      do iz = 1,izmax
        sum1 = 0.
        sum1p = 0.
        sum2 = 0.
        sum2p = 0.
c
        do j = narn1 + 1,narn2
          zj = zchar(j)
          jz = nint(zj)
          iajz = abs(jz)
          if (jz .gt. 0) then
            sum1 = sum1 + elam(iz,iajz)*conc(j)
            sum1p = sum1p + delam(1,iz,iajz)*conc(j)
          elseif (jz .lt. 0) then
            sum2 = sum2 + elam(iz,iajz)*conc(j)
            sum2p = sum2p + delam(1,iz,iajz)*conc(j)
          endif
        enddo
c
        selm(iz) = sum1
        dselm(1,iz) = sum1p
        selm(-iz) = sum2
        dselm(1,-iz) = sum2p
      enddo
c
      end
