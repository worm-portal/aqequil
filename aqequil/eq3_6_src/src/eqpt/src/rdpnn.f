      subroutine rdpnn(abeta,ipbtmx,jpfcmx,nat,natmax,ndat0s,nerr,
     $ noutpt,npxnn,npx2mx,npx2t,nttyo,nwarn,uaqsp,upair,zaqsp)
c
c     This subroutine reads from the DATA1 file the coefficients
c     required to compute those Pitzer interaction parameters
c     associated with neutral-different neutral (nn) pairs. This
c     particular set of interaction coefficients includes the
c     following:
c
c       1. lambda(nn')
c
c     Although mu(nnn') and mu(nn'n') are closely related to
c     lambda(nn'), they are not read from the current superblock
c     by this subroutine, but from a triplets superblock read by
c     another subroutine, rdpn2n.f. There are several reasons for
c     this treatment. First, the relevant type of mu here is more
c     analogous to a psi parameter than a Cphi parameter. Note that
c     mu(nnn) is included with lambda(nn) in part because mu(nnn)
c     is analogous to a Cphi parameter. Second, data for the current
c     type of mu is both scarce and not very likely to be needed
c     in most applications.
c
c     Coefficients required to compute Pitzer interaction coefficients
c     associated with other pair types or with triplets are not read
c     by this subroutine. This subroutine is one of several that
c     replace rdpz2.f and rdpz3.f.
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
c       abeta  = array of coefficients for computing beta parameters
c       nerr   = cumulative error counter
c       npxnn  = the number of nn pairs for which Pitzer parameters
c                  were read from the data file
c       npx2t  = the number of pairs of all species types for which
c                  Pitzer parameters were read from the data file
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
      integer nat,nerr,npxnn,npx2t,nwarn
c
      character(len=24) uaqsp(natmax),upair(2,npx2mx)
c
      real(8) abeta(jpfcmx,0:ipbtmx,npx2mx),zaqsp(natmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ier,iz1,iz2,j,j2,j3,j4,j5,j6,npx2,npx2te,n1,n2
c
      integer ilnobl
c
      character(len=80) udastr,uline,ux80
      character(len=24) unam1,unam2,ux24a,ux24b
      character(len=24) uhdnn,uhdnxt
      character(len=16) ustr16
      character(len=8) uterm,ux8
c
      real(8) var,z1,z2
c
c-----------------------------------------------------------------------
c
      data uterm  /'+-------'/
      data uhdnn  /"nn' combinations: lambda"/
      data uhdnxt /"cc'a and aa'c combinatio"/
c
c-----------------------------------------------------------------------
c
      npx2te = npx2t
      npx2 = npx2te
      j2 = 0
      j3 = 0
c
c     Advance to the start of the superblock for these parameters.
c
  100 read (ndat0s,1000,end=990,err=995) uline
 1000 format(a80)
      if (uline(1:24) .ne. uhdnn(1:24)) go to 100
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
      if (unam1(1:24) .eq. uhdnxt(1:24)) then
c
c       Found the header for the next superblock of Pitzer parameters.
c       Back up one line and terminate reading the current superblock.
c
        npxnn = npx2 - npx2te
        npx2t = npx2
        backspace(ndat0s)
        go to 999
      endif
c
      npx2 = npx2 + 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Setup up some string data for error messages.
c
      j2 = ilnobl(unam1)
      if (j2 .le. 0) then
        unam1 = '<blank>'
        j2 = ilnobl(unam1)
      endif
c
      j3 = ilnobl(unam2)
      if (j3 .le. 0) then
        unam2 = '<blank>'
        j3 = ilnobl(unam2)
      endif
c
      write (ux8,'(i5)') npx2
      call lejust(ux8)
      j4 = ilnobl(ux8)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for blank names.
c
      if (j2.le.0 .or. j3.le.0) then
        write (noutpt,1020) ux8(1:j4),unam1(1:j2),unam2(1:j3)
        write (nttyo,1020) ux8(1:j4),unam1(1:j2),unam2(1:j3)
 1020   format(/' * Error - (EQPT/rdpnn) Have encountered a blank',
     $  ' species name',/7x,'in the species pair for block ',a,
     $  ' of the superblock for',/7x,"Pitzer nn' parameter data.",
     $  ' The offending pair shows as ',a,', ',a,'.')
c
        if (npx2 .gt. (npx2te + 1)) then
          ux24a = upair(1,npx2 - 1)
          ux24b = upair(2,npx2 - 1)
          if (ux24a(1:7).ne.'<blank>' .or.
     $      ux24b(1:7).ne.'<blank>') then
            j5 = ilnobl(ux24a)
            j6 = ilnobl(ux24b)
            write (noutpt,1030) ux24a(1:j5),ux24b(1:j6)
            write (nttyo,1030) ux24a(1:j5),ux24b(1:j6)
 1030       format(7x,'This block follows the one for ',a,
     $      ', ',a,'.')
          endif
        endif
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
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
          write (noutpt,1100) unam1(1:j2),unam1(1:j2),unam2(1:j3)
          write (nttyo,1100) unam1(1:j2),unam1(1:j2),unam2(1:j3)
 1100     format(/' * Error - (EQPT/rdpnn) The aqueous species ',
     $    a,' appearing in',/7x,"the Pitzer nn' parameter data block",
     $    ' for the species pair',/7x,a,', ',a,' does not appear',
     $    ' in an aqueous species',/7x,'data block.')
          nerr = nerr + 1
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
          write (noutpt,1100) unam2(1:j3),unam1(1:j2),unam2(1:j3)
          write (nttyo,1100) unam2(1:j3),unam1(1:j2),unam2(1:j3)
          nerr = nerr + 1
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for illegal appearance of solvent water.
c
      if (unam1(1:4).eq.'H2O ' .or. unam2(1:4).eq.'H2O ' .or.
     $  n1.eq.1 .or. n2.eq.1) then
        write (noutpt,1200) unam1(1:j2),unam2(1:j3)
        write (nttyo,1200) unam1(1:j2),unam2(1:j3)
 1200   format(/' * Error - (EQPT/rdpnn) Have found an illegal',
     $  ' data block for the',/7x,'species pair ',a,', ',a,
     $  ' in the superblock of Pitzer data',/7x,"for nn' pairs.",
     $  ' Solvent water may not appear in such a pair.')
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the electrical charges of the species.
c
      z1 = 0.
      z2 = 0.
      iz1 = 0
      iz2 = 0
      if (n1 .gt. 0) then
        z1 = zaqsp(n1)
        iz1 = nint(z1)
      endif
      if (n2 .gt. 0) then
        z2 = zaqsp(n2)
        iz2 = nint(z2)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (n1.gt.0 .and. n2.gt.0) then
c
c       Make sure that the species pair is appropriate for the
c       current superblock. It must be of type nn '. Other
c       pair types such as ca, cc', aa', nc, na, and nn are not
c       permitted here.
c
        if (iz1.ne.0 .or. iz2.ne.0 .or. n1.eq.n2) then
c
c         The current data block is in the wrong superblock.
c
          write (noutpt,1230) unam1(1:j2),unam2(1:j3)
          write (nttyo,1230) unam1(1:j2),unam2(1:j3)
 1230     format(/' * Error - (EQPT/rdpnn) Have found an',
     $    ' illegal data block for the',/7x,'species pair ',
     $    a,', ',a,' in the',/7x,'superblock of Pitzer data for',
     $    " nn' pairs.",/7x,'The present pair type may not',
     $    ' appear in the current superblock.')
          nerr = nerr + 1
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Species pair storage rules:
c
c       1. neutral < cation < anion
c
c       2. If both species are of the same charge type,
c          store alphabetically.
c
c     Here: alphabetical.
c
      if (lle(unam1,unam2)) then
c
c       The first species is first alphabetically.
c
        upair(1,npx2) = unam1
        upair(2,npx2) = unam2
      else
c
c       The second species is first alphabetically.
c
        upair(1,npx2) = unam2
        upair(2,npx2) = unam1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the coefficients for the lambda(nc) and lambda(na)
c     parameters. Store the lambda data as beta(0) data.
c
c     Check the required parameter header.
c
      read (ndat0s,1000,end=990,err=995) uline
      j4 = ilnobl(uline)
      j4 = min(j4,70)
      ux80 = uline
      call lejust(ux80)
      ustr16 = 'lambda'
      j5 = 6
      if (ux80(1:j5) .ne. ustr16(1:j5)) then
        write (noutpt,1330) uline(1:j4),ustr16(1:j5),unam1(1:j2),
     $  unam2(1:j3)
        write (nttyo,1330) uline(1:j4),ustr16(1:j5),unam1(1:j2),
     $  unam2(1:j3)
 1330   format(/' * Error - (EQPT/rdpnn) Have found a line',
     $  ' starting with:',/7x,'"',a,'"',/7x,'where "',a,'" was',
     $  ' expected in the block for the species',/7x,'pair ',
     $  a,', ',a,' in the superblock of Pitzer data',
     $  /7x,"for nn' pairs.")
        nerr = nerr + 1
        go to 999
      endif
c
c     Read the coefficients for the associated temperature function.
c
      do j = 1,jpfcmx
        read (ndat0s,1000,end=990,err=995) uline
        j4 = ilnobl(uline)
        j4 = min(j4,70)
        ux80 = uline
        call lejust(ux80)
        ustr16 = 'a  ='
        j5 = 4
        write (ustr16(2:2),'(i1)') j
        if (ux80(1:j5) .eq. ustr16(1:j5)) then
          udastr = ux80
          udastr(1:j5) = ' '
          call g1dat(ier,noutpt,nttyo,udastr,var)
          if (ier .gt. 0) then
            write (noutpt,1320) uline(1:j4),unam1(1:j2),unam2(1:j3)
            write (nttyo,1320) uline(1:j4),unam1(1:j2),unam2(1:j3)
 1320       format(/' * Error - (EQPT/rdpnn) Have found a line',
     $      ' starting with:',/7x,'"',a,'"',/7x,'containing an',
     $      ' expected numerical field that could not be read.',
     $      /7x,'This occurred in the block for the species pair',
     $      /7x,a,', ',a,' in the superblock of Pitzer data',
     $      /7x,"for nn' pairs.")
            nerr = nerr + 1
            go to 999
          endif
          abeta(j,0,npx2) = var
        else
          write (noutpt,1330) uline(1:j4),ustr16(1:j5),unam1(1:j2),
     $    unam2(1:j3)
          write (nttyo,1330) uline(1:j4),ustr16(1:j5),unam1(1:j2),
     $    unam2(1:j3)
          nerr = nerr + 1
          go to 999
        endif
      enddo
c
c     Skip past the block delimiter line.
c
      read (ndat0s,1000,end=990,err=995) uline
c
      if (uline(1:8) .ne. uterm(1:8)) then
        j4 = ilnobl(uline)
        j4 = min(j4,70)
        write (noutpt,1470) uline(1:j4),unam1(1:j2),unam2(1:j3)
        write (nttyo,1470) uline(1:j4),unam1(1:j2),unam2(1:j3)
 1470   format(/' * Error - (EQPT/rdpnn) Have found a line starting',
     $  ' with:',/7x,'"',a,'"',/7x,' in the data block for the species',
     $  ' pair ',a,', ',a,' in the',/7x,'superblock of Pitzer data for',
     $  " nn' pairs. Should have",/7x,'found the delimiter line',
     $  ' marking the end of the block.')
        nerr = nerr + 1
        go to 999
      endif
c
c     Process the next data block.
c
      go to 110
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 write (noutpt,2010)
      write (nttyo,2010)
 2010 format(/' * Error - (EQPT/rdpnn) Unexpectedly encountered',
     $ /7x,"an end-of-file error while reading the nn' (lambda)",
     $ /7x,'superblock of the DATA0 file.')
      stop
c
  995 write (noutpt,2020)
      write (nttyo,2020)
 2020 format(/' * Error - (EQPT/rdpnn) Encountered a read format',
     $ /7x,"error while reading the nn' (lambda) superblock of",
     $ /7x,'the DATA0 file.')
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
