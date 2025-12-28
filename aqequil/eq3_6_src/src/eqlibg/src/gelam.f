      subroutine gelam(aphi,delam,dpelm,elam,fxi,izmax,nazpmx,
     $ pelm,qpit75)
c
c     This subroutine calculates the E-lambda function (elam(i,j)) and
c     its first two ionic strength derivatives (delam(1,i,j) and
c     delam(2,i,j)). Here i, j refers to the charge pair zi, zj or
c     -zi, -zj, where zi and zj are both positive. The E-lambda
c     function and its derivatives are used in Pitzer's equations
c     to represent higher-order electrical interactions. Here izmax
c     is the max norm of the electrical charges of the aqueous species
c     and aphi is the Debye-Huckel A(phi) parameter. The array
c     pelm(i,j) contains a set of "primitive" E-lambda functions
c     corresponding to those in elam(i,j). The arrays dpelm(1,i,j)
c     and dpelm(2,i,j) similarly parallel delam(1,i,j) and delam(2,i,j).
c
c     This subroutine is called by:
c
c       EQLIBG/gcoeff.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       aphi   = the Debye-Huckel A(phi) parameter
c       izmax  = the max norm of the electrical charge numbers of
c                  the aqueous species
c       fxi    = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
c
c     Principal output:
c
c       elam   = array of values of the E-lambda functions
c       delam  = array of values of the ionic strength derivatives
c                  of the E-lambda functions
c
c     Work space:
c
c       pelm   = array of values of primitive E-lambda functions
c       dpelm  = array of values of the ionic strength derivatives
c                  of the primitive E-lambda functions
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nazpmx
c
      integer izmax
c
      logical qpit75
c
      real*8 delam(2,nazpmx,nazpmx),dpelm(2,nazpmx,nazpmx),
     $ elam(nazpmx,nazpmx),pelm(nazpmx,nazpmx)
      real*8 aphi,fxi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ijz,iz,jz
c
      real*8 el,elp,elpp,wi,wj
c
c-----------------------------------------------------------------------
c
c     First calculate "primitive" E-lambdas (pelm) and their
c     derivatives (dpelm) for the various charge pairs. These depend
c     only on the product of the electrical charges.
c
      do jz = 1,izmax
        do iz = jz,izmax
          ijz = iz*jz
c
c         Get the primitive E-lambda and its derivatives for
c         this current charge pair.
c
          call elmdd(aphi,el,elp,elpp,fxi,ijz,qpit75)
          pelm(iz,jz) = el
          dpelm(1,iz,jz) = elp
          dpelm(2,iz,jz) = elpp
          if (jz .ne. iz) then
            pelm(jz,iz) = el
            dpelm(1,jz,iz) = elp
            dpelm(2,jz,iz) = elpp
          endif
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Now calculate conventional E-lambdas (elam) and their derivatives
c     (delam) for the various charge pairs.
c
      do jz = 1,izmax
        do iz = 1,izmax
          wj = 0.5*iz/jz
          wi = 0.5*jz/iz
          elam(iz,jz) = pelm(iz,jz) - wj*pelm(jz,jz) - wi*pelm(iz,iz)
          delam(1,iz,jz) = dpelm(1,iz,jz) - wj*dpelm(1,jz,jz)
     $    - wi*dpelm(1,iz,iz)
          delam(2,iz,jz) = dpelm(2,iz,jz) - wj*dpelm(2,jz,jz)
     $    - wi*dpelm(2,iz,iz)
        enddo
      enddo
c
c     Set the diagonal elements to zero.
c
      do i = 1,izmax
        elam(i,i) = 0.
        delam(1,i,i) = 0.
        delam(2,i,i) = 0.
      enddo
c
      end
