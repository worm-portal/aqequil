      subroutine evratc(eact,hact,iact,imchmx,imech,nrct,nrctmx,nrk,
     $ rk,rkb,rtcnst,tempk,trkb)
c
c     This subroutine evaluates rate constants as functions of
c     temperature. There are two alternative treatments. One assumes
c     a constant activation energy. The other assumes a constant
c     activation enthalpy.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
c       EQ6/tpadv.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer imchmx,nrctmx
c
      integer iact(imchmx,2,nrctmx),imech(2,nrctmx),nrk(2,nrctmx)
c
      integer nrct
c
      real*8 eact(imchmx,2,nrctmx),hact(imchmx,2,nrctmx),
     $ rk(imchmx,2,nrctmx),rkb(imchmx,2,nrctmx),trkb(imchmx,2,nrctmx)
c
      real*8 rtcnst,tempk
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,nrc
c
      real*8 tfx,tk1
c
c-----------------------------------------------------------------------
c
      do nrc = 1,nrct
c
c       Forward rate laws (forward equates to dissolution,
c       dissociation, or other means of destruction).
c
        if (nrk(1,nrc) .ge. 1) then
          do i = 1,imech(1,nrc)
            if (iact(i,1,nrc) .eq. 0) then
c
c             No temperature dependence.
c
              rk(i,1,nrc) = rkb(i,1,nrc)
            elseif (iact(i,1,nrc) .eq. 1) then
c
c             Constant activation energy.
c
              tk1 = trkb(i,1,nrc) + 273.15
              tfx = (eact(i,1,nrc)/rtcnst) * ((tempk/tk1) - 1.)
              rk(i,1,nrc) = rkb(i,1,nrc)*exp(tfx)
            elseif (iact(i,1,nrc) .eq. 2) then
c
c             Constant activation enthalpy.
c
              tk1 = trkb(i,1,nrc) + 273.15
              tfx = (hact(i,1,nrc)/rtcnst) * ((tempk/tk1) - 1.)
              rk(i,1,nrc) = rkb(i,1,nrc)*(tempk/tk1)*exp(tfx)
            endif
          enddo
        endif
c
c       Backward rate laws (backward equates to precipitation,
c       association, or other means of formation).
c
        if (nrk(2,nrc) .ge. 1) then
          do i = 1,imech(2,nrc)
            if (iact(i,2,nrc) .eq. 0) then
c
c             No temperature dependence.
c
              rk(i,2,nrc) = rkb(i,2,nrc)
            elseif (iact(i,2,nrc) .eq. 1) then
c
c             Constant activation energy.
c
              tk1 = trkb(i,2,nrc) + 273.15
              tfx = (eact(i,2,nrc)/rtcnst) * ((tempk/tk1) - 1.)
              rk(i,2,nrc) = rkb(i,2,nrc)*exp(tfx)
            elseif (iact(i,2,nrc) .eq. 2) then
c
c             Constant activation enthalpy.
c
              tk1 = trkb(i,2,nrc) + 273.15
              tfx = (hact(i,2,nrc)/rtcnst) * ((tempk/tk1) - 1.)
              rk(i,2,nrc) = rkb(i,2,nrc)*(tempk/tk1)*exp(tfx)
            endif
          enddo
        endif
      enddo
c
      end
