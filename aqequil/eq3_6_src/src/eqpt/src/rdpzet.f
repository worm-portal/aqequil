      subroutine rdpzet(apsi,jpfcmx,nat,natmax,ndat0s,nerr,noutpt,
     $ npxzet,npx3mx,npx3t,nttyo,nwarn,uaqsp,utripl,zaqsp)
c
c     This subroutine reads from the DATA1 file the coefficients
c     required to compute those Pitzer interaction parameters
c     associated with neutral-cation-anion (nca) triplets. This
c     particular set of interaction coefficients includes the
c     following:
c
c       1. zeta(nca)
c
c     Coefficients required to compute Pitzer interaction coefficients
c     associated with other triplet types or with pairs are not read
c     by this subroutine. This subroutine is one of several that
c     replace rdpz2.f and rdpz3.f.
c
c     This suboutine is called by:
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
c       apsi   = array of coefficients for computing psi parameters
c       nerr   = cumulative error counter
c       npxzet = the number of nca triplets for which Pitzer parameters
c                 were read from the data file
c       npx3t  = the number of triplets of all species types for which
c                  Pitzer parameters were read from the data file.
c       nwarn  = cumulative warning counter
c       uaqsp  = array of names of aqueous species
c       utripl = array of species names for species triplets
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
      integer jpfcmx,natmax,npx3mx
c
      integer ndat0s,noutpt,nttyo
c
      integer nat,nerr,npxzet,npx3t,nwarn
c
      character(len=24) uaqsp(natmax),utripl(3,npx3mx)
c
      real(8) apsi(jpfcmx,npx3mx),zaqsp(natmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ier,iz1,iz2,iz3,j,j2,j3,j4,j5,j6,j7,j8,j9,
     $ na,nc,nn,npx3,npx3te,n1,n2,n3
c
      integer ilnobl
c
      logical qdupan,qdup12,qdup13,qdup23
c
      character(len=80) udastr,uline,ux80
      character(len=24) unam1,unam2,unam3,ux24a,ux24b,ux24c,u1,u2,u3
      character(len=24) uhdzet,uhdnxt
      character(len=16) ustr16
      character(len=8) uterm,ux8
c
      real(8) var,z1,z2,z3
c
c-----------------------------------------------------------------------
c
      data uterm  /'+-------'/
      data uhdzet /"nca combinations: zeta(n"/
      data uhdnxt /"nnn' combinations: mu(nn"/
c
c-----------------------------------------------------------------------
c
      npx3te = npx3t
      npx3 = npx3te
      j2 = 0
      j3 = 0
      j4 = 0
c
c     Advance to the start of the superblock for these parameters.
c
  100 read (ndat0s,1000,end=990,err=995) uline
 1000 format(a80)
      if (uline(1:24) .ne. uhdzet(1:24)) go to 100
c
c     Skip the terminator.
c
      read (ndat0s,1000,end=990,err=995) uline
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the names of the species comprising a triplet.
c
  110 read (ndat0s,1010,end=990,err=995) unam1,unam2,unam3
 1010 format(a24,2x,a24,2x,a24)
c
      if (unam1(1:24) .eq. uhdnxt(1:24)) then
c
c       Found the header for the next superblock of Pitzer parameters.
c       Back up one line and terminate reading the current superblock.
c
        npxzet = npx3 - npx3te
        npx3t = npx3
        backspace(ndat0s)
        go to 999
      endif
c
      npx3 = npx3 + 1
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
      j4 = ilnobl(unam3)
      if (j4 .le. 0) then
        unam3 = '<blank>'
        j4 = ilnobl(unam3)
      endif
c
      write (ux8,'(i5)') npx3
      call lejust(ux8)
      j5 = ilnobl(ux8)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for blank names.
c
      if (j2.le.0 .or. j3.le.0 .or. j4.le.0) then
        write (noutpt,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3),
     $  unam3(1:j3)
        write (nttyo,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3),
     $  unam3(1:j3)
 1020   format(/' * Error - (EQPT/rdpzet) Have encountered a blank',
     $  ' species name',/7x,'in the species triplet for block ',a,
     $  ' of the superblock for',/7x,"Pitzer nca parameter data.",
     $  ' The offending triplet shows as',/7x,a,', ',a,', ',a,'.')
c
        if (npx3 .gt. (npx3te + 1)) then
          ux24a = utripl(1,npx3 - 1)
          ux24b = utripl(2,npx3 - 1)
          ux24c = utripl(3,npx3 - 1)
          if (ux24a(1:7).ne.'<blank>' .or.
     $      ux24b(1:7).ne.'<blank>' .or. ux24c(1:7).ne.'<blank>') then
            j7 = ilnobl(ux24a)
            j8 = ilnobl(ux24b)
            j9 = ilnobl(ux24c)
            write (noutpt,1030) ux24a(1:j7),ux24b(1:j8),ux24c(1:j9)
            write (nttyo,1030) ux24a(1:j7),ux24b(1:j8),ux24c(1:j9)
 1030       format(7x,'This block follows the one for ',a,
     $      ', ',a,', ',a,'.')
          endif
        endif
        nerr = nerr + 1
      endif
c
cXXXXXXXXX
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the indices of the three species (n1, n2, n3). These
c     indices point to the data read from the previously read
c     species blocks.
c
c     Calling sequence substitutions:
c       n1 for na
c       unam1 for unams
c
      call gspidx(ier,n1,nat,natmax,uaqsp,unam1)
c
      if (ier .gt. 0) then
        if (unam1(1:7).ne.'<blank>' .and.
     $    unam2(1:7).ne.'<blank>' .and. unam3(1:7).ne.'<blank>') then
          write (noutpt,1100) unam1(1:j2),unam1(1:j2),unam2(1:j3),
     $    unam3(1:j4)
          write (nttyo,1100) unam1(1:j2),unam1(1:j2),unam2(1:j3),
     $    unam3(1:j4)
 1100     format(/' * Error - (EQPT/rdpzet) The aqueous species ',
     $    a,' appearing in',/7x,"the Pitzer nca parameter data",
     $    ' block for the species triplet',/7x,a,', ',a,', ',a,
     $    ' does not appear in an',/7x,'aqueous species data block.')
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
     $    unam2(1:7).ne.'<blank>' .and. unam3(1:7).ne.'<blank>') then
          write (noutpt,1100) unam2(1:j3),unam1(1:j2),unam2(1:j3),
     $    unam3(1:j4)
          write (nttyo,1100) unam2(1:j3),unam1(1:j2),unam2(1:j3),
     $    unam3(1:j4)
          nerr = nerr + 1
        endif
      endif
c
c     Calling sequence substitutions:
c       n3 for na
c       unam3 for unams
c
      call gspidx(ier,n3,nat,natmax,uaqsp,unam3)
c
      if (ier .gt. 0) then
        if (unam1(1:7).ne.'<blank>' .and.
     $    unam2(1:7).ne.'<blank>' .and. unam3(1:7).ne.'<blank>') then
          write (noutpt,1100) unam3(1:j3),unam1(1:j2),unam2(1:j3),
     $    unam3(1:j4)
          write (nttyo,1100) unam3(1:j3),unam1(1:j2),unam2(1:j3),
     $    unam3(1:j4)
          nerr = nerr + 1
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for illegal appearance of solvent water.
c
      if (unam1(1:4) .eq.'H2O ' .or. unam2(1:4) .eq.'H2O ' .or.
     $  unam3(1:4).eq.'H2O ' .or. n1.eq.1 .or. n2.eq.1 .or.
     $  n3.eq.1) then
        write (noutpt,1200) unam1(1:j2),unam2(1:j3),unam3(1:j4)
        write (nttyo,1200) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1200   format(/' * Error - (EQPT/rdpzet) Have found an illegal data',
     $  ' block for the',/7x,'species triplet ',a,', ',a,', ',a,
     $  ' in the superblock',/7x,"of Pitzer data for nca triplets.",
     $  ' Solvent water may not',/7x,'appear in such a triplet.')
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the electrical charges of the species.
c
      z1 = 0.
      z2 = 0.
      z3 = 0.
      iz1 = 0
      iz2 = 0
      iz3 = 0
      if (n1 .gt. 0) then
        z1 = zaqsp(n1)
        iz1 = nint(z1)
      endif
      if (n2 .gt. 0) then
        z2 = zaqsp(n2)
        iz2 = nint(z2)
      endif
      if (n3 .gt. 0) then
        z3 = zaqsp(n3)
        iz3 = nint(z3)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (n1.ne.0 .and. n2.ne.0 .and. n3.ne.0) then
c
c       Make sure that the species triplet is appropriate for the
c       current superblock. It must be of type nca. Other triplet
c       types such as cc'a, aa'c, nnn, nnn', and nn'n' are not
c       permitted here.
c
c       Count the neutral species.
c
        nn = 0
        if (iz1 .eq. 0) nn = 1
        if (iz2 .eq. 0) nn = nn + 1
        if (iz3 .eq. 0) nn = nn + 1
c
c       Count the cations.
c
        nc = 0
        if (iz1 .gt. 0) nc = 1
        if (iz2 .gt. 0) nc = nc + 1
        if (iz3 .gt. 0) nc = nc + 1
c
c       Count the anions.
c
        na = 0
        if (iz1 .lt. 0) na = 1
        if (iz2 .lt. 0) na = na + 1
        if (iz3 .lt. 0) na = na + 1
c
c       Check for duplications.
c
        qdup12 = unam1(1:24) .eq. unam2(1:24)
        qdup13 = unam1(1:24) .eq. unam3(1:24)
        qdup23 = unam2(1:24) .eq. unam3(1:24)
        qdupan = qdup12 .or. qdup23 .or. qdup13
c
        if (nc.eq.3 .or. na.eq.3
     $    .or. (nc.eq.2 .and. nn.eq.1)
     $    .or. (na.eq.2 .and. nn.eq.1)
     $    .or. (nn.eq.2 .and. nc.eq.1)
     $    .or. (nn.eq.2 .and. na.eq.1)
     $    .or. (nc.eq.2 .and. na.eq.1 .and. qdupan)
     $    .or. (na.eq.2 .and. nc.eq.1 .and. qdupan)) then
c
c         The current data block refers to a triplet combination
c         that is not valid in the Pitzer framework.
c
          write (noutpt,1220) unam1(1:j2),unam2(1:j3),unam3(1:j4)
          write (nttyo,1220) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1220     format(/' * Error - (EQPT/rdpzet Have found an',
     $    ' illegal data block for the',/7x,'species triplet ',
     $    a,', ',a,', ',a,' in the',/7x,'superblock of Pitzer',
     $    " data for nca triplets.",/7x,'The present triplet',
     $    ' type is not valid in the Pitzer framework.')
          nerr = nerr + 1
        elseif (.not.(nn.eq.1 .and. nc.eq.1 .and. na.eq.1)) then
c
c         The current data block is in the wrong superblock.
c
          write (noutpt,1230) unam1(1:j2),unam2(1:j3),unam3(1:j4)
          write (nttyo,1230) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1230     format(/' * Error - (EQPT/rdpzet Have found an',
     $    ' illegal data block for the',/7x,'species triplet ',
     $    a,', ',a,', ',a,' in the',/7x,'superblock of Pitzer',
     $    " data nca triplets. The present",/7x,'triplet type',
     $    ' may not appear in the current superblock.')
          nerr = nerr + 1
        endif
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Species triplet storage rules:
c
c       nnn : neutral, neutral, neutral
c       nnn': neutral 1, neutral1, neutral 2
c       nca : neutral, cation, anion
c       cc'a: cation 1, cation 2, anion
c       aa'c: anion 1, anion 2, cation
c
c       A duplicated species appears first.
c
c       Otherwise, if two species are of the same charge type,
c       store alphabetically.
c
c     Here: neutral, cation, anion.
c
c     Copy the names before rearranging for storage.
c
      u1 = unam1
      u2 = unam2
      u3 = unam3
c
c     Rearrange according to the storage rules. Note that
c     n1, n2, n3, z1, z2, z3, iz1, iz2, iz3 are all changed
c     in addtion to u1, u2, u3. Note that unam1, unam2, and
c     unam3 are not changed.
c
      call artrip(iz1,iz2,iz3,na,nc,nn,n1,n2,n3,
     $ u1,u2,u3,z1,z2,z3)
c
c     Store the triplet names in the order required by the
c     storage rules. The triplet in the original order is
c     preserved in unam1, unam2, unam3 for use in error messages
c     and such.
c
      utripl(1,npx3) = u1
      utripl(2,npx3) = u2
      utripl(3,npx3) = u3
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the coefficients for the zeta(nca) parameters. Store the
c     zeta data as psi data.
c
c     Check the required parameter header.
c
      read (ndat0s,1000,end=990,err=995) uline
      j5 = ilnobl(uline)
      j5 = min(j5,70)
      ux80 = uline
      call lejust(ux80)
      ustr16 = 'zeta'
      j6 = 4
      if (ux80(1:j6) .ne. ustr16(1:j6)) then
        write (noutpt,1330) uline(1:j5),ustr16(1:j6),unam1(1:j2),
     $  unam2(1:j3),unam3(1:j4)
        write (nttyo,1330) uline(1:j5),ustr16(1:j6),unam1(1:j2),
     $  unam2(1:j3),unam3(1:j4)
 1330   format(/' * Error - (EQPT/rdpzet) Have found a line',
     $  ' starting with:',/7x,'"',a,'"',/7x,'where "',a,'" was',
     $  ' expected in the block for the species',/7x,'triplet ',
     $  a,', ',a,', ',a,' in the superblock of Pitzer data',
     $  /7x,"for nca triplets.")
        nerr = nerr + 1
        go to 999
      endif
c
c     Read the coefficients for the associated temperature function.
c
      do j = 1,jpfcmx
        read (ndat0s,1000,end=990,err=995) uline
        j5 = ilnobl(uline)
        j5 = min(j5,70)
        ux80 = uline
        call lejust(ux80)
        ustr16 = 'a  ='
        j6 = 4
        write (ustr16(2:2),'(i1)') j
        if (ux80(1:j6) .eq. ustr16(1:j6)) then
          udastr = ux80
          udastr(1:j6) = ' '
          call g1dat(ier,noutpt,nttyo,udastr,var)
          if (ier .gt. 0) then
            write (noutpt,1320) uline(1:j5),unam1(1:j2),unam2(1:j3),
     $      unam3(1:j4)
            write (nttyo,1320) uline(1:j5),unam1(1:j2),unam2(1:j3),
     $      unam3(1:j4)
 1320       format(/' * Error - (EQPT/rdpzet) Have found a line',
     $      ' starting with:',/7x,'"',a,'"',/7x,'containing an',
     $      ' expected numerical field that could not be read.',
     $      /7x,'This occurred in the block for the species triplet',
     $      /7x,a,', ',a,', ',a,' in the superblock of Pitzer data',
     $      /7x,"for nca triplets.")
            nerr = nerr + 1
            go to 999
          endif
          apsi(j,npx3) = var
        else
          write (noutpt,1330) uline(1:j5),ustr16(1:j6),unam1(1:j2),
     $    unam2(1:j3),unam3(1:j4)
          write (nttyo,1330) uline(1:j5),ustr16(1:j6),unam1(1:j2),
     $    unam2(1:j3),unam3(1:j4)
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
        j5 = ilnobl(uline)
        j5 = min(j5,70)
        write (noutpt,1470) uline(1:j5),unam1(1:j2),unam2(1:j3),
     $  unam3(1:j4)
        write (nttyo,1470) uline(1:j5),unam1(1:j2),unam2(1:j3),
     $  unam3(1:j4)
 1470   format(/' * Error - (EQPT/rdpzet) Have found a line starting',
     $  ' with:',/7x,'"',a,'"',/7x,' in the data block for the species',
     $  ' triplet ',a,', ',a,', ',a,/7x,'in the superblock of',
     $  " Pitzer data for nca triplets. Should have",
     $  /7x,'found the delimiter line marking the end of the block.')
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
 2010 format(/' * Error - (EQPT/rdpzet) Unexpectedly encountered',
     $ /7x,"an end-of-file error while reading the nca (zeta)",
     $ /7x,'superblock of the DATA0 file.')
      stop
c
  995 write (noutpt,2020)
      write (nttyo,2020)
 2020 format(/' * Error - (EQPT/rdpzet) Encountered a read format',
     $ /7x,"error while reading the nca (zeta) superblock of",
     $ /7x,'the DATA0 file.')
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
c
      end
