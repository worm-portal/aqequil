      subroutine betars(alphar,betar,btrfnc,btrmax,btrmxo,ibtrmx,
     $ nrct,nrct1,nrctmx,rirec1,rirecp,rrelr1,rrelrp)
c
c     This subroutine calculates residual functions for the ODE
c     integrator.
c
c     This subroutine is called by:
c
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c       alphar = array of alpha[r] residuals
c       betar  = array of beta[r] residuals
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nrctmx
c
      integer ibtrmx,nrct,nrct1
c
      real(8) alphar(nrct1),betar(nrct1),rrelr1(nrctmx),rrelrp(nrctmx)
c
      real(8) btrfnc,btrmax,btrmxo,rirec1,rirecp
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nrc
c
      integer iarmxn
c
c-----------------------------------------------------------------------
c
c     Calculate absolute (alphar) and normalized (betar) residuals.
c
      do nrc = 1,nrct
        alphar(nrc) = rrelr1(nrc) - rrelrp(nrc)
        betar(nrc) = 0.
        if (rrelrp(nrc) .ne. 0.) betar(nrc) = alphar(nrc)/rrelrp(nrc)
      enddo
c
      alphar(nrct1) = rirec1 - rirecp
      if (rirecp .ne. 0.) betar(nrct1) = alphar(nrct1)/rirecp
c
c     Find the max norm of the normalized residual vector (betar).
c
      ibtrmx = iarmxn(betar,nrct1)
      btrmax = 0.
      if (ibtrmx .gt. 0) btrmax = abs(betar(ibtrmx))
c
c     Calculate the betar improvement function (btrfnc).
c
      btrfnc = 0.
      if (btrmxo .gt. 0.) btrfnc = (btrmxo -btrmax)/btrmxo
c
      btrmxo = btrmax
c
      end
