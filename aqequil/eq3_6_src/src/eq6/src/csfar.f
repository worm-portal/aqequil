      subroutine csfar(afrc1,morr,morr0,mwtrc,noutpt,nrc,nrctmx,
     $ nsk,nttyo,prcinf,sfcar,sfcar0,ssfcar,ureac)
c
c     This subroutine calculates the surface area required to calculate
c     the rate for the nrc-th irreversible reaction.
c
c     This subroutine is called by:
c
c       EQ6/rtcalc.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
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
      integer noutpt,nttyo
c
      integer nsk(nrctmx)
c
      integer nrc
c
      character(len=24) ureac(nrctmx)
c
      real(8) afrc1(nrctmx),morr(nrctmx),morr0(nrctmx),mwtrc(nrctmx),
     $ sfcar(nrctmx),sfcar0(nrctmx),ssfcar(nrctmx)
c
      real(8) prcinf
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2
c
      integer ilnobl
c
      real(8) gx
c
c-----------------------------------------------------------------------
c
c     Compute surface area.
c
      gx = morr(nrc)*mwtrc(nrc)
      if (nsk(nrc) .eq. 0) then
c
c       Constant surface area (cm2).
c       Compute the specific surface area (cm2/g).
c
        ssfcar(nrc) = 0.
        if (gx .gt. 0.) ssfcar(nrc) = sfcar(nrc)/gx
      elseif (nsk(nrc) .eq. 1) then
c
c       Constant specific surface area (cm2/g).
c       Compute the surface area (cm2).
c
        if (morr(nrc) .gt. 0.) then
          sfcar(nrc) = ssfcar(nrc)*gx
        elseif (afrc1(nrc) .le. 0.) then
          sfcar(nrc) = 1.e+5
          j2 = ilnobl(ureac(nrc))
          write (noutpt,1000) ureac(nrc)(1:j2),sfcar(nrc)
          write (nttyo,1000) ureac(nrc)(1:j2),sfcar(nrc)
 1000     format(/' * Note - (EQ6/csfar) The surface area of',
     $    ' reactant ',a,/7x,'has been temporarily set to ',1pe11.4,
     $    ' cm2 to permit this phase',/7x,'to begin precipitating.',
     $    ' This is necessary because the surface area',/7x,
     $    ' model requires the phase to have a mass, which it',
     $    /7x,'presently lacks.')
        else
          sfcar(nrc) = 0.
        endif
      elseif (nsk(nrc) .eq. 2) then
c
c       Constant particle number surface area growth law.
c       Compute the surface area (cm2).
c       Compute the specific surface area (cm2/g).
c
        if (morr0(nrc) .gt. 0.) then
          sfcar(nrc) = sfcar0(nrc)*(morr(nrc)/morr0(nrc))**(2./3.)
          ssfcar(nrc) = 0.
          if (gx .gt. 0) ssfcar(nrc) = sfcar(nrc)/gx
        elseif (afrc1(nrc) .le. 0.) then
          sfcar(nrc) = 1.e+5
          ssfcar(nrc) = prcinf
          j2 = ilnobl(ureac(nrc))
          write (noutpt,1010) ureac(nrc)(1:j2),sfcar(nrc)
          write (nttyo,1010) ureac(nrc)(1:j2),sfcar(nrc)
 1010     format(/' * Note - (EQ6/csfar) The surface area of',
     $    ' reactant ',a,/7x,'has been temporarily set to ',1pe11.4,
     $    ' cm2 to permit this phase',/7x,'to begin precipitating.',
     $    ' This is necessary because the surface area',/7x,
     $    ' model requires the phase to have a mass at the',
     $    /7x,'previous point, which it lacks.')
        else
          sfcar(nrc) = 0.
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
