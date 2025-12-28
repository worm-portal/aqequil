      subroutine crrate(act,afrc1,cdac,csigma,eps100,fkrc,idirec,
     $ imchmx,imech,iodb,jreac,morr,ndac,ndact,ndctmx,nodbmx,noutpt,
     $ nrc,nrctmx,nrk,nstmax,nttyo,rk,rreac1,rrelr1,rrxfi1,rtcnst,
     $ sfcar,udac,ureac)
c
c     This subroutine calculates the rate (relative or absolute) from
c     the specified rate law for the nrc-th irreversible reaction.
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
      integer imchmx,ndctmx,nodbmx,nrctmx,nstmax
c
      integer noutpt,nttyo
c
      integer nrc
c
      integer idirec(nrctmx),imech(2,nrctmx),iodb(nodbmx),jreac(nrctmx),
     $ ndac(ndctmx,imchmx,2,nrctmx),ndact(imchmx,2,nrctmx),nrk(2,nrctmx)
c
      character*24 udac(ndctmx,imchmx,2,nrctmx),ureac(nrctmx)
c
      real(8) act(nstmax),afrc1(nrctmx),cdac(ndctmx,imchmx,2,nrctmx),
     $ csigma(imchmx,2,nrctmx),fkrc(nrctmx),morr(nrctmx),
     $ rk(imchmx,2,nrctmx),rreac1(nrctmx),rrelr1(nrctmx),
     $ rrxfi1(imchmx,nrctmx),sfcar(nrctmx)
c
      real(8) eps100,rtcnst
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j2,j3,n,ns
c
      integer ilnobl
c
      logical qovstp
c
      character(len=8) ux8
c
      real(8) affr,aprod,efx,fs,rx
c
c-----------------------------------------------------------------------
c
c     Calculate the net rate (relative or absolute) for each
c     kinetically-governed reaction. Use either a forward rate constant
c     plus affinity form or a backward rate constant plus affinity form.
c     This subroutine does not use a form involving a forward rate
c     constant plus backward rate constant.
c
      rreac1(nrc) = 0.
      rrelr1(nrc) = 0.
      do i = 1,imchmx
        rrxfi1(i,nrc) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for exhausted or saturated reactant.
c
      if (jreac(nrc) .gt. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the sign of the affinity. For minerals, the affinity
c     here is the affinity to dissolve, so:
c
c       Negative value -> precipitation
c       Positive value -> dissolution
c
c     Some reactants (special reactants, aqeuous species, and
c     gases) currently do not have actual affinities associated
c     with them. For these, the affinity is currently taken as
c     +9999999.
c
      affr = afrc1(nrc)
      qovstp = affr.lt.0. .and. nrk(2,nrc).eq.0 .and. jreac(nrc).le.0
c
      if ( iodb(2) .ge. 2 ) then
        j2 = ilnobl(ureac(nrc))
        if (affr.gt.0. .or. qovstp) then
          write (noutpt,1010) ureac(nrc)(1:j2)
 1010     format(/' --- Rate calculations for destruction of ',a,
     $    ' ---',/)
        elseif (affr .lt. 0.) then
          write (noutpt,1020) ureac(nrc)(1:j2)
 1020     format(/' --- Rate calculations for formation of ',a,
     $    ' ---',/)
        else
          write (noutpt,1030) sfcar(nrc)
 1030     format(/5x,'Surface area= ',e12.5,' cm2.',/)
        endif
      endif
c
      fs = fkrc(nrc)*sfcar(nrc)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (idirec(nrc) .eq. 1) then
c
c       Calculate the net rate using the forward (dissolution) form
c       of rate law.
c
        if (nrk(1,nrc) .eq. 0) then
          j2 = ilnobl(ureac(nrc))
          write (noutpt,1100) ureac(nrc)(1:j2)
          write (nttyo,1100) ureac(nrc)(1:j2)
 1100     format(/' * Error - (EQ6/crrate) Programming error trap:',
     $    /7x,'The forward rate law code (nrk(1,nrc)) is zero for',
     $    /7x,'reactant ',a,'.')
          stop
        elseif (nrk(1,nrc) .eq. 1) then
c
c         Specified relative rate (arbitrary kinetics).
c
          rrelr1(nrc) = rk(1,1,nrc)
          rrxfi1(1,nrc) = rrelr1(nrc)
c
        elseif (nrk(1,nrc) .eq. 2) then
c
c         Transition state theory (more properly, the TST-like form).
c         The net rate is equal to an absolute forward rate minus an
c         absolute backward rate. However, that is obvious only in the
c         equivalent two rate constant form.
c
          rreac1(nrc) = 0.
          do i = 1,imech(1,nrc)
            efx = affr/(csigma(i,1,nrc)*rtcnst)
            aprod = 1.
            do n = 1,ndact(i,1,nrc)
              ns = ndac(n,i,1,nrc)
              aprod = aprod*act(ns)**cdac(n,i,1,nrc)
            enddo
            rx = rk(i,1,nrc)*fs*aprod*(1.0 - exp(-efx))
            rreac1(nrc) = rreac1(nrc) + rx
            rrxfi1(i,nrc) = rx
          enddo
c
        elseif (nrk(1,nrc) .eq. 3) then
c
c         Specified rate (no affinity dependence).
c
          rreac1(nrc) = fs*rk(1,1,nrc)
          rrxfi1(1,nrc) = rreac1(nrc)
c
        else
c
c         Unrecognized rate law code.
c
          j2 = ilnobl(ureac(nrc))
          write (ux8,'(i5)') nrk(1,nrc)
          call lejust(ux8)
          j3 = ilnobl(ux8)
          write (noutpt,1110) ureac(nrc)(1:j2),ux8(1:j3)
          write (nttyo,1110) ureac(nrc)(1:j2),ux8(1:j3)
 1110     format(/' * Error - (EQ6/crrate) The reaction for the',
     $    /7x,'destruction of reactant ',a,' has an unrecognized',
     $    /7x,'forward rate law code of ',a,'.')
          stop
        endif
      else
c
c       Calculate the net rate using the backward (formation) form
c       of rate law.
c
c       Note: All rates are expressed in this code as the net of
c       forward minus backward rates, so a postive formation rate
c       is expressed as a negative number.
c
        if (nrk(2,nrc). eq. 0) then
c
c         Instantaneous equilibrium. Escape for case of a
c         supersaturated reactant that is not in the matrix.
c
        elseif (nrk(2,nrc) .eq. 1) then
c
c         Specified relative rate (aribitrary kinetics).
c
          rrelr1(nrc) = -rk(1,2,nrc)
          rrxfi1(1,nrc) = rrelr1(nrc)
c
        elseif (nrk(2,nrc) .eq. 2) then
c
c         Transition state theory (more properly, the TST-like form).
c
          rreac1(nrc) = 0.
          do i = 1,imech(2,nrc)
c
c           Note: if affr < 0, then efx < 0, affac < 0, and rreac1 < 0.
c
            efx = affr/(csigma(i,2,nrc)*rtcnst)
            aprod = 1.
            do n = 1,ndact(i,2,nrc)
              ns = ndac(n,i,2,nrc)
              aprod = aprod*act(ns)**cdac(n,i,2,nrc)
            enddo
            rx = rk(i,2,nrc)*fs*aprod*(1.0 - exp(-efx))
            rreac1(nrc) = rreac1(nrc) + rx
            rrxfi1(i,nrc) = rx
          enddo
c
        elseif (nrk(2,nrc) .eq. 3) then
c
c         Linear rate law.
c
          rreac1(nrc) = -fs*rk(1,2,nrc)
          rrxfi1(1,nrc) = rreac1(nrc)
c
        else
c
          j2 = ilnobl(ureac(nrc))
          write (ux8,'(i5)') nrk(2,nrc)
          call lejust(ux8)
          j3 = ilnobl(ux8)
          write (noutpt,250) ureac(nrc)(1:j2),ux8(1:j3)
          write (nttyo,250) ureac(nrc)(1:j2),ux8(1:j3)
  250     format(/' * Error - (EQ6/crrate) The reaction for the',
     $    /7x,'destruction of reactant ',a,' has an unrecognized',
     $    /7x,'backward rate law code of ',a,'.')
          stop
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
