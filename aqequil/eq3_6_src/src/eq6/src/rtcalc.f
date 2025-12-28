      subroutine rtcalc(act,afrc1,cdac,csigma,eps100,fkrc,idirec,
     $ imchmx,imech,iodb,iopt,jcode,jreac,morr,morr0,mwtrc,ndac,
     $ ndact,ndctmx,nodbmx,noptmx,nord,noutpt,nrk,nrct,nrctmx,nsk,
     $ nstmax,nttyo,prcinf,prminf,qriinf,rirec1,rk,rreac1,rrelr1,
     $ rtcnst,rrxfi1,sfcar,sfcar0,ssfcar,udac,ureac)
c
c     This subroutine calculates the relative and absolute rates of
c     the nrc-th irreversible reaction. This rate is computed from the
c     specified rate expression.
c
c     Under iopt(2) = 1, xi1 is now defined as the sum of the absolute
c     values of the progress variables for all kinetically governed
c     (nrk(1,nrc).gt.1 or nrk(2,nrc).gt.1) reactions. The inverse rate
c     is the inverse of the sum of the absolute values of the rates of
c     these reactions.
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
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer imchmx,ndctmx,nodbmx,noptmx,nrctmx,nstmax
c
      integer idirec(nrctmx),imech(2,nrctmx),iodb(nodbmx),iopt(noptmx),
     $ jcode(nrctmx),jreac(nrctmx),ndac(ndctmx,imchmx,2,nrctmx),
     $ ndact(imchmx,2,nrctmx),nrk(2,nrctmx),nsk(nrctmx)
c
      integer noutpt,nttyo
c
      integer nord,nrct
c
      logical qriinf
c
      character*24 udac(ndctmx,imchmx,2,nrctmx),ureac(nrctmx)
c
      real(8) act(nstmax),afrc1(nrctmx),cdac(ndctmx,imchmx,2,nrctmx),
     $ csigma(imchmx,2,nrctmx),fkrc(nrctmx),morr(nrctmx),morr0(nrctmx),
     $ mwtrc(nrctmx),rk(imchmx,2,nrctmx),rreac1(nrctmx),rrelr1(nrctmx),
     $ rrxfi1(imchmx,nrctmx),sfcar(nrctmx),sfcar0(nrctmx),
     $ ssfcar(nrctmx)
c
      real*8 eps100,prcinf,prminf,rirec1,rtcnst
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer id,j2,nrc
c
      integer ilnobl
c
      logical qovstp
c
      real*8 affr,trate
c
c-----------------------------------------------------------------------
c
c     Loop over all kinetically-governed reactions. Calculate the
c     relevant surface areas for surface-area-controlled reactions.
c
      do nrc = 1,nrct
        call csfar(afrc1,morr,morr0,mwtrc,noutpt,nrc,nrctmx,
     $  nsk,nttyo,prcinf,sfcar,sfcar0,ssfcar,ureac)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Determine the rate law form (forward or backward) to use for
c     each reactant.
c
      do nrc = 1,nrct
        affr = afrc1(nrc)
        qovstp = affr.lt.0. .and. nrk(2,nrc).eq.0 .and. jreac(nrc).le.0
        if (affr .ge. 0.) then
c
c         The affinity favors the forward direction.
c
          idirec(nrc) = 1
          if (nrk(1,nrc) .eq. -1) idirec(nrc) = 2
        else
c
c         The affinity favors the backward direction.
c
          idirec(nrc) = 2
          if (nrk(2,nrc).eq.-1 .or. qovstp) idirec(nrc) = 1
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Loop over all kinetically-governed reactions. Calculate the
c     net rate (relative or absolute) for each one.
c
      do nrc = 1,nrct
        call crrate(act,afrc1,cdac,csigma,eps100,fkrc,idirec,
     $  imchmx,imech,iodb,jreac,morr,ndac,ndact,ndctmx,nodbmx,noutpt,
     $  nrc,nrctmx,nrk,nstmax,nttyo,rk,rreac1,rrelr1,rrxfi1,rtcnst,
     $  sfcar,udac,ureac)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     If in time mode, calculate the inverse rate. Calculate relative
c     rates from the corresponding absolute rates. If any primary
c     calculated rates were relative rates, calculate the corresponding
c     absolute rates.
c
      if (iopt(2) .gt. 0) then
c
c       Compute the total rate.
c
        trate = 0.
        do nrc = 1,nrct
          id = idirec(nrc)
          if (nrk(id,nrc) .gt. 1) trate = trate + abs(rreac1(nrc))
        enddo
c
c       Compute absolute rates for reactants constrained by
c       relative rates.
c
        do nrc = 1,nrct
          id = idirec(nrc)
          if (nrk(id,nrc) .eq. 1) rreac1(nrc) = rrelr1(nrc)*trate
        enddo
c
c       Compute the inverse rate.
c
        qriinf = .false.
        if (trate .le. prminf) then
          qriinf = .true.
          nord = 0
          rirec1 = prcinf
        else
          rirec1 = 1./trate
        endif
c
c       Compute relative rates for reactants constrained by
c       absolute rates.
c
        do nrc = 1,nrct
          id = idirec(nrc)
          if (nrk(id,nrc) .gt. 1) rrelr1(nrc) = rreac1(nrc)*rirec1
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check if the rate and the affinity are compatible.
c
      do nrc = 1,nrct
        if (jreac(nrc).le.0 .and. jcode(nrc) .le. 1) then
          affr = afrc1(nrc)
          if (affr.gt.0. .and. rrelr1(nrc).lt.0.) then
            j2 = ilnobl(ureac(nrc))
            write (noutpt,1000) ureac(nrc)(1:j2)
            write (nttyo,1000) ureac(nrc)(1:j2)
 1000       format(/' * Warning - (EQ6/rtcalc) The calculated rate',
     $      ' for reactant',/7x,a,' is less than zero, but the',
     $      ' affinity',/7x,'is greater than zero:')
            write (noutpt,1010) rreac1(nrc),rrelr1(nrc),affr
            write (nttyo,1010) rreac1(nrc),rrelr1(nrc),affr
 1010       format(/9x,'Rate = ',1pe12.5,' mol/d',
     $      /9x,'Relative rate= ',e12.5,' mol/mol',
     $      /9x,'Affinity= ',e12.5,' kcal')
            write (noutpt,1020)
            write (nttyo,1020)
 1020       format(/7x,'Check the signs of the rate constants.',/)
          endif
c
          qovstp = affr.lt.0. .and. nrk(2,nrc).eq.0 .and.
     $    jreac(nrc).le.0
          if (affr.lt.0. .and. .not.qovstp .and.
     $      rrelr1(nrc).gt.0.) then
            j2 = ilnobl(ureac(nrc))
            write (noutpt,1050) ureac(nrc)(1:j2)
            write (nttyo,1050) ureac(nrc)(1:j2)
 1050       format(/' * Warning - (EQ6/rtcalc) The calculated rate',
     $      ' for reactant',/7x,a,' is greater than zero, but the',
     $      ' affinity',/7x,'is less than zero:')
            write (noutpt,1010) rreac1(nrc),rrelr1(nrc),affr
            write (nttyo,1010) rreac1(nrc),rrelr1(nrc),affr
            write (noutpt,1020)
            write (nttyo,1020)
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(2) .gt. 1) then
c
c       Write a summary of calculated rate information on the
c       output file.
c
        do nrc = 1,nrct
          if (jreac(nrc) .le. 0) then
            j2 = ilnobl(ureac(nrc))
            write (noutpt,1070) ureac(nrc)(1:j2),rreac1(nrc),
     $      rrelr1(nrc),affr
 1070       format(/5x,'Calculated Rate Summary',
     $      /5x,a,/7x,'Rate = ',1pe12.5,' mol/d',
     $      /7x,'Relative rate= ',e12.5,' mol/mol',
     $      /7x,'Affinity= ',e12.5,' kcal')
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
