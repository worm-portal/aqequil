      subroutine rdpz2(abeta,alpha,acphi,ipbtmx,jpfcmx,nat,natmax,
     $ ndat0s,nerr,noutpt,npx2mx,npx2t,nttyo,nwarn,uaqsp,upair,
     $ zaqsp)
c
c     This subroutine reads coefficients for computing the following
c     Pitzer interaction parameters from the DATA0 file:
c
c       beta(MX)(0), beta(MX)(1), beta(MX)(2), and their
c         temperature derivatives
c       Cphi(MX) and its temperature derivatives
c
c     (Here M = a cation and X = an anion). The alpha(1) and alpha(2)
c     coefficients which go with the beta parameters are also read.
c
c     The beta parameters map to conventional lambda parameters for
c     cation-anion pairs. Conventional lambda coefficients for MN
c     and NX pairs, where N = a neutral species may be entered in
c     place of beta parameters. The same is true for the observable
c     lambda coefficients of NN and NN' pairs (here N' is a second
c     neutral species).
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ndat0s = unit number of the stripped DATA0 file
c       ndat1f = unit number of the DATA1F file
c
c     Principal output:
c
c       alpha  = array of Pitzer alpha constants
c       acphi  = array of coefficients for computing Cphi parameters
c       abeta  = array of coefficients for computing beta parameters
c       nerr   = cumulative error counter
c       npx2t  = the number of pairs of species for which Pitzer
c                  interaction parameters are defined
c       nwarn  = cumulative warning counter
c       uaqsp  = array of names of aqueous species
c       upair  = array of species names for species pairs
c       uaqsp  = array of electrical charge numbers of aqueous species
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipbtmx,jpfcmx,natmax,npx2mx
c
      integer ndat0s,noutpt,nttyo
c
      integer nat,nerr,npx2t,nwarn
c
      character(len=24) uaqsp(natmax),upair(2,npx2mx)
c
      real(8) abeta(jpfcmx,0:ipbtmx,npx2mx),alpha(ipbtmx,npx2mx),
     $ acphi(jpfcmx,npx2mx),zaqsp(natmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ier,iz1,iz1x,iz2,iz2x,j,j2,j3,j4,j5,j6,j7,npx2,n1,n2
c
      integer ilnobl
c
      logical qzepfc
c
      character(len=80) uline
      character(len=24) unam1,unam2,ux24a,ux24b
      character(len=8) umixt,usingl,uterm,ux8,ux8a
c
      real(8) z1,z2
c
c-----------------------------------------------------------------------
c
      data uterm  /'+-------'/
      data usingl /'single-s'/
      data umixt  / 'mixture'/
c
c-----------------------------------------------------------------------
c
      npx2 = 0
c
c     Advance to the start of the superblock for these parameters.
c
  100 read (ndat0s,1000,end=990,err=995) uline
 1000 format(a80)
      if (uline(1:8) .ne. usingl(1:8)) go to 100
c
c     Skip the terminator.
c
      read (ndat0s,1000,end=990,err=995) uline
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the names of the species comprising a pair.
c
  110 read (ndat0s,1010,end=990,err=995) unam1,unam2
 1010 format(a24,2x,a24)
c
      if (unam1(1:8) .eq. umixt(1:8)) then
c
c       Found the header for the mixture parameters superblock.
c
        npx2t = npx2
        go to 999
      endif
c
      npx2 = npx2 + 1
c
c     Check for blank names.
c
      j2 = ilnobl(unam1)
      j3 = ilnobl(unam2)
      if (j2.le.0 .or. j3.le.0) then
        write (ux8,'(i5)') npx2
        call lejust(ux8)
        j5 = ilnobl(ux8)
        if (j2 .le. 0) then
          unam1 = '<blank>'
          j2 = ilnobl(unam1)
        endif
        if (j3 .le. 0) then
          unam2 = '<blank>'
          j3 = ilnobl(unam2)
        endif
        write (noutpt,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
        write (nttyo,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
 1020   format(/' * Error - (EQPT/rdpz2) Have encountered a blank',
     $  ' species name',/7x,'in the species pair for block ',a,
     $  ' of the superblock for',/7x,'Pitzer parameter data',
     $  ' corresponding to species pairs.',/7x,'The offending',
     $  ' pair shows as ',a,', ',a,'.')
        if (npx2 .gt. 1) then
          ux24a = upair(1,npx2 - 1)
          ux24b = upair(2,npx2 - 1)
          if (ux24a(1:7).ne.'<blank>' .or.
     $      ux24b(1:7).ne.'<blank>') then
            j6 = ilnobl(ux24a)
            j7 = ilnobl(ux24b)
            write (noutpt,1030) ux24a(1:j6),ux24b(1:j7)
            write (nttyo,1030) ux24a(1:j6),ux24b(1:j7)
 1030       format(7x,'This block follows the one for ',a,
     $      ', ',a,'.')
          endif
        endif
        nerr = nerr + 1
      endif
c
c     Read the electrical charges.
c
      read (ndat0s,1050,end=990,err=995) z1,z2
 1050 format(f3.0,t15,f3.0)
c
      n1 = 0
      n2 = 0
      iz1 = 0
      iz2 = 0
c
c     Get the indices of the two species (n1, n2). These
c     indices point to the data read from the previously
c     read species blocks.
c
c     Calling sequence substitutions:
c       n1 for na
c       unam1 for unams
c
      call gspidx(ier,n1,nat,natmax,uaqsp,unam1)
c
      if (ier .gt. 0) then
        if (unam1(1:7).ne.'<blank>' .and.
     $    unam2(1:7).ne.'<blank>') then
          j2 = ilnobl(unam1)
          j3 = ilnobl(unam2)
          write (noutpt,1100) unam1(1:j2),unam1(1:j2),unam2(1:j3)
          write (nttyo,1100) unam1(1:j2),unam1(1:j2),unam2(1:j3)
 1100     format(/' * Error - (EQPT/rdpz2) The aqueous species ',
     $    a,' appearing in',/7x,'the Pitzer parameter data block',
     $    ' for the species pair',/7x,a,', ',a,' does not appear',
     $    ' in an aqueous species',/7x,'data block.')
          nerr = nerr + 1
        endif
      else
c
        iz1x = nint(z1)
        z1 = zaqsp(n1)
        iz1 = nint(z1)
c
        if (iz1 .ne. iz1x) then
          j2 = ilnobl(unam1)
          j3 = ilnobl(unam2)
          write (ux8a,'(i5)') iz1x
          call lejust(ux8a)
          j4 = ilnobl(ux8a)
          write (ux8,'(i5)') iz1
          call lejust(ux8)
          j5 = ilnobl(ux8)
          write (noutpt,1110) ux8a(1:j4),unam1(1:j2),
     $    unam1(1:j2),unam2(1:j3),ux8(1:j5)
          write (nttyo,1110) ux8a(1:j4),unam1(1:j2),
     $    unam1(1:j2),unam2(1:j3),ux8(1:j5)
 1110     format(/' * Note - (EQPT/rdpz2) A charge of ',a,
     $    ' was incorrectly',/7x,'specifed for ',a,' in the'
     $    ' data block for the pair',/7x,a,', ',a,'. The',
     $    ' correct charge is ',a,'.',/7x,'The incorrect',
     $    ' charge will be ignored.')
        endif
      endif
c
c     Calling sequence substitutions:
c       n2 for na
c       unam2 for unams
c
      call gspidx(ier,n2,nat,natmax,uaqsp,unam2)
c
      if (ier .gt. 0) then
        if (unam1(1:7).ne.'<blank>' .and.
     $    unam2(1:7).ne.'<blank>') then
          j2 = ilnobl(unam1)
          j3 = ilnobl(unam2)
          write (noutpt,1100) unam2(1:j3),unam1(1:j2),unam2(1:j3)
          write (nttyo,1100) unam2(1:j3),unam1(1:j2),unam2(1:j3)
          nerr = nerr + 1
        endif
      else
c
        iz2x = nint(z2)
        z2 = zaqsp(n2)
        iz2 = nint(z2)
c
        if (iz2 .ne. iz2x) then
          j2 = ilnobl(unam1)
          j3 = ilnobl(unam2)
          write (ux8a,'(i5)') iz2x
          call lejust(ux8a)
          j4 = ilnobl(ux8a)
          write (ux8,'(i5)') iz2
          call lejust(ux8)
          j5 = ilnobl(ux8)
          write (noutpt,1110) ux8a(1:j4),unam2(1:j3),
     $    unam1(1:j2),unam2(1:j3),ux8(1:j5)
          write (nttyo,1110) ux8a(1:j4),unam2(1:j3),
     $    unam1(1:j2),unam2(1:j3),ux8(1:j5)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for illegal appearance of solvent water.
c
      if (unam1(1:4).eq.'H2O ' .or. unam2(1:4).eq.'H2O ' .or.
     $  n1.eq.1 .or. n2.eq.1) then
        j2 = ilnobl(unam1)
        j3 = ilnobl(unam2)
        write (noutpt,1200) unam1(1:j2),unam2(1:j3)
        write (nttyo,1200) unam1(1:j2),unam2(1:j3)
 1200   format(/' * Error - (EQPT/rdpz2) Have found an illegal',
     $  ' data block for the',/7x,'species pair ',a,', ',a,
     $  ' among the blocks for "single-salt',/7x,'parameters".',
     $  ' Solvent water may not appear in such a pair in the',
     $  /7x,'normal Pitzer treatment of electrolyte solutions.')
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (n1.gt.0 .and. n2.gt.0) then
c
c       Check for illegal combinations (cc, aa, cc', and aa').
c       Note that nc, na, and nn' combinations are allowed here
c       under restricted circumstances (see further below), even
c       though they technically correspond to mixtures.
c
        if ((iz1*iz2) .gt. 0) then
          j2 = ilnobl(unam1)
          j3 = ilnobl(unam2)
          if (n1 .ne. n2) then
            ux8 = 'cations'
            if (iz1 .lt. 0) ux8 = 'anions'
            j5 = ilnobl(ux8)
            write (noutpt,1220) unam1(1:j2),unam2(1:j3),ux8(1:j5)
            write (nttyo,1220) unam1(1:j2),unam2(1:j3),ux8(1:j5)
 1220       format(/' * Error - (EQPT/rdpz2) Have found an',
     $      ' illegal data block for the',/7x,'species pair ',
     $      a,', ',a,' among the blocks for "single-salt',
     $      /7x,'parameters". Both species are ',a,'. This is',
     $      ' not a valid',/7x,'combination for the present kind',
     $      ' of block. Enter these',/7x,'data in a block for',
     $      '"mixture parameters".')
          else
            write (noutpt,1230) unam1(1:j2),unam2(1:j3)
            write (nttyo,1230) unam1(1:j2),unam2(1:j3)
 1230       format(/' * Error - (EQPT/rdpz2) Have found an',
     $      ' illegal data block for the',/7x,'species pair ',
     $      a,', ',a,' among the blocks for "single-salt',
     $      /7x,'parameters". The same ion appears twice.',
     $      ' This is not',/7x,'a valid combination, as the',
     $      ' corresponding lambda and mu',/7x,'data are',
     $      ' defined to be zero by convention in the normal',
     $      ' Pitzer',/7x,'treatment of electrolyte solutions.')
          endif
          nerr = nerr + 1
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Species pair storage rules:
c
c       neutral < cation < anion
c       if both species are of the same charge type,
c       store alphabetically.
c
      if (iz1.eq.0 .and. iz2.ne.0) then
c
c       The first species is a neutral and the second species
c       isn't.
c
        upair(1,npx2) = unam1
        upair(2,npx2) = unam2
      elseif (iz2.eq.0 .and. iz1.ne.0) then
c
c       The second species is a neutral and the first species
c       isn't.
c
        upair(1,npx2) = unam2
        upair(2,npx2) = unam1
      elseif (iz1.gt.0 .and. iz2.lt.0) then
c
c       The first species is a cation and the second species
c       is an anion.
c
        upair(1,npx2) = unam1
        upair(2,npx2) = unam2
      elseif (iz2.gt.0 .and. iz1.lt.0) then
c
c       The second species is a cation and the first species
c       is an anion.
c
        upair(1,npx2) = unam2
        upair(2,npx2) = unam1
      else
c
c       Both species have the same charge type.
c
        if (unam1 .le. unam2) then
          upair(1,npx2) = unam1
          upair(2,npx2) = unam2
        else
          upair(1,npx2) = unam2
          upair(2,npx2) = unam1
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the beta (lambda) coefficients.
c
      read (ndat0s,1300,end=990,err=995) abeta(1,0,npx2),
     $ abeta(1,1,npx2),abeta(1,2,npx2)
 1300 format(13x,f9.5,11x,f9.5,11x,f9.5)
c
c     Read the alpha parameters.
c
      read (ndat0s,1310,end=990,err=995) alpha(1,npx2),alpha(2,npx2)
 1310 format(18x,2(16x,f5.1))
c
c     Read the Cphi parameter.
c
      read (ndat0s,1320,end=990,err=995) acphi(1,npx2)
 1320 format(13x,f9.5,12x,f5.1)
c
c     Read source 1.
c
      read (ndat0s,1330,end=990,err=995) uline
 1330 format(13x,a24)
c
c     Read the beta (lambda) derivatives.
c
      read (ndat0s,1340,end=990,err=995) abeta(2,0,npx2),
     $ abeta(3,0,npx2)
      read (ndat0s,1340,end=990,err=995) abeta(2,1,npx2),
     $ abeta(3,1,npx2)
      read (ndat0s,1340,end=990,err=995) abeta(2,2,npx2),
     $ abeta(3,2,npx2)
 1340 format(13x,e10.3,13x,e10.3)
c
c     Read the Cphi derivatives.
c
      read (ndat0s,1340,end=990,err=995) acphi(2,npx2),
     $ acphi(3,npx2)
c
c     Read source 2.
c
      read (ndat0s,1330,end=990,err=995) uline
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (n1.gt.0 .and. n2.gt.0) then
c
        if ((iz1.eq.0 .and. iz2.eq.0) .and. (n1 .ne. n2)) then
c
c         Have an nn' combination. This is really a mixture
c         mixture parameter combination. Only lambda coefficients
c         (as "beta0") are permitted in a "single-salt parameters"
c         block. The coresponding two mu coefficients (as "psi")
c         must be entered in separate blocks under "mixture
c         parameters". The lambda data must be entered in a
c         block under "single-salt parameters" in order to
c         preserve current block formatting and to avoid the
c         problem that a given lambda of this type corresponds to
c         two distinct mu values. If such a lambda were to be
c         entered in the "theta" field of the blocks for the
c         corresponding mu values, then the same parameter would
c         be entered in two places.
c
          j2 = ilnobl(unam1)
          j3 = ilnobl(unam2)
          qzepfc = .true.
          do j = 1,jpfcmx
            if (acphi(j,npx2) .ne. 0.) qzepfc = .false.
          enddo
          if (.not. qzepfc) then
            write (noutpt,1400) unam1(1:j2),unam2(1:j3)
            write (nttyo,1400) unam1(1:j2),unam2(1:j3)
 1400       format(/' * Error - (EQPT/rdpz2) Have found an',
     $      ' illegal data block for the',/7x,'species pair ',
     $      a,', ',a,' among the blocks for "single-salt',
     $      /7x,'parameters". Have non-zero input for mu',
     $      ' data in the "Cphi" field.',/7x,'Only lambda data',
     $      ' may be specified (in the "beta0" field) in',
     $      /7x,'a block of this type. To specify mu data,',
     $      ' enter each of the',/7x,'two corresponding mu',
     $      " values [mu(nnn') and mu(n'n'n)] in separate",
     $      /7x,'"mixture parameters" blocks, using the "psi"',
     $      ' fields. The',/7x,'corresponding lambda data must',
     $      ' remain in a "single-salt parameters"',/7x,'block',
     $      ' (in the "beta0" field).')
            nerr = nerr + 1
          endif
        endif
c
        if ((iz1.ne.0 .and. iz2.eq.0) .or.
     $    (iz1.eq.0 .and. iz2.ne.0)) then
c
c         Have an nc or na combination. This is also really a
c         mixture parameter combination. Only lambda coefficients
c         (as "beta0") are permitted in a "single-salt parameters"
c         block. The coresponding zeta coefficient (as "psi")
c         must be entered in a block under "mixture parameters".
c         Note that there is a lambda(nc) and a lambda(na) for
c         each zeta(nca). The lambda data must be entered in
c         separate blocks under "single-salt parameters" because
c         a given lambda of this type is shared by various nca-
c         type mixtures and duplicate values would occur
c         if the lambdas were put into blocks under "mixture
c         parameters". For example, lambda(Na+, CO2(aq)) pertains
c         to mixtures involving various anions.
c
          j2 = ilnobl(unam1)
          j3 = ilnobl(unam2)
          qzepfc = .true.
          do j = 1,jpfcmx
            if (acphi(j,npx2) .ne. 0.) qzepfc = .false.
          enddo
          if (.not. qzepfc) then
            write (noutpt,1410) unam1(1:j2),unam2(1:j3)
            write (nttyo,1410) unam1(1:j2),unam2(1:j3)
 1410       format(/' * Error - (EQPT/rdpz2) Have found an',
     $      ' illegal data block for the',/7x,'species pair ',
     $      a,', ',a,' among the blocks for "single-salt',
     $      /7x,'parameters". Have non-zero input for zeta',
     $      ' data in the "Cphi"',/7x,'field. Only lambda data',
     $      ' may be specified in a block of this type,',
     $      /7x,'(in the "beta0" field). To specify zeta data,',
     $      ' enter it in',/7x,'a "mixture parameters"',
     $      ' block in the "psi" field. The',/7x,'corresponding',
     $      ' lambda data must remain in a "single-salt',
     $      /7x,'parameters" block (in the "beta0" field).')
            nerr = nerr + 1
          endif
        endif
c
        if (iz1*iz2 .ge. 0) then
c
c         Have a combination other than a ca combination.
c         The "beta1" and "beta2" parameters are ordinarily used
c         only for ca combinations. Issue a warning if they
c         are used for any other combinations, such as nc, na,
c         nn, or nn'.
c
          j2 = ilnobl(unam1)
          j3 = ilnobl(unam2)
          qzepfc = .true.
          do i = 1,ipbtmx
            do j = 1,jpfcmx
              if (abeta(j,i,npx2) .ne. 0.) qzepfc = .false.
            enddo
          enddo
          if (.not. qzepfc) then
            write (noutpt,1420) unam1(1:j2),unam2(1:j3)
            write (nttyo,1420) unam1(1:j2),unam2(1:j3)
 1420       format(/' * Warning - (EQPT/rdpz2) Have found a',
     $      ' data block for the',/7x,'species pair ',
     $      a,', ',a,' among the blocks for "single-salt',
     $      /7x,'parameters" that has non-zero input in the',
     $      ' "beta1" and/or "beta2"',/7x,'fields. The "beta1"',
     $      ' and "beta2" parameters are ordinarily used',
     $      /7x,'only for cation-anion pairs.')
            nwarn = nwarn + 1
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Skip past the block delimiter line.
c
  120 read (ndat0s,1500,end=990,err=995) uline
 1500 format(a8)
      if (uline(1:8) .eq. uterm(1:8)) then
c
c       Process the next data block.
c
        go to 110
      endif
      go to 120
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 write (noutpt,2000)
      write (nttyo,2000)
 2000 format(/' * Error - (EQPT/rdpz2) Unexpectedly encountered',
     $ /7x,'end-of-file while reading the DATA0 file.')
      stop
c
  995 write (noutpt,2010)
      write (nttyo,2010)
 2010 format(/' * Error - (EQPT/rdpz2) Encountered a read format',
     $ /7x,'error while reading the DATA0 file.')
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
