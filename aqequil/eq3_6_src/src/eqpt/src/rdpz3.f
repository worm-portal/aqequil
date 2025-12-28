      subroutine rdpz3(apsi,atheta,jpfcmx,nat,natmax,ndat0s,nerr,
     $ noutpt,npx3mx,npx3t,nthdt,nttyo,nwarn,uaqsp,uethfl,uthdtr,
     $ utripl,zaqsp)
c
c     This suboutine reads the following Pitzer interaction parameter
c     data from the DATA0 file:
c
c       S-theta(MM'X) (and S-theta (MXX')) and their temperature
c         derivatives
c       psi(MM'X) (and psi (MXX')) and their temperature derivatives
c
c     (Here M = a cation, M' = a second cation, X = an anion, and
c     X' = a second anion).
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
c
c     Principal output:
c
c       apsi   = array of coefficients for computing psi parameters
c       atheta = array of coefficients for computing S-theta parameters
c       nerr   = cumulative error counter
c       nwarn  = cumulative warning counter
c       nthdt  = the number of distinct S-theta parameter entries
c                  (i.e., not counting duplicates)
c       npx3t  = the number of triplets of species for which Pitzer
c                  interaction parameters are defined
c       utripl = array of names in a species triplet
c       uthdtr = array of names of cc'a and aa'c triplets which
c                  provide the values used for coefficients for
c                  computing S-theta parameters for the cc' and aa'
c                  pairs
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
      integer nat,nerr,npx3t,nthdt,nwarn
c
      character*24 uaqsp(natmax),uthdtr(3,npx3mx),utripl(3,npx3mx)
      character*8 uethfl
c
      real*8 atheta(jpfcmx,npx3mx),apsi(jpfcmx,npx3mx),zaqsp(natmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ier,ifound,iz1,iz2,iz3,j,j2,j3,j4,j5,j6,j7,j8,
     $ na,nc,nn,npx3,nthd,n1,n2,n3
c
      integer ilnobl
c
      logical qdup12,qdup13,qdup23
c
      character*80 uline
      character*24 unam1,unam2,unam3,uthd1,uthd2,uthd3,ux24a,ux24b,
     $ ux24c,u1,u2,u3
      character*8 uelemt,uterm,ux8
c
      real(8), dimension(:), allocatable :: athetx
c
      real(8) z1,z2,z3
c
c-----------------------------------------------------------------------
c
      data uelemt / 'elements'/
      data uterm  / '+-------'/
c
c-----------------------------------------------------------------------
c
c     Allocate work array.
c
      ALLOCATE(athetx(jpfcmx))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      nthdt = 0
      npx3 = 0
c
c     Read the E-theta flag.
c
      uethfl = ' '
      read (ndat0s,1000,end=990,err=995) uline
 1000 format(a80)
      i = index(uline,'E-theta flag')
      if (i .eq. 0) then
        write (noutpt,1010)
        write (nttyo,1010)
 1010   format(/' * Error - (EQPT/rdpz3) The E-theta flag line is',
     $  ' missing from',
     $  /7x,'the DATA0 file. This line is specific to DATA0 files',
     $  /7x,"tied to Pitzer's equations. It should appear after",
     $  /7x,'the line containing "mixture parameters". The E-theta',
     $  /7x,'flag line should contain "E-theta flag = on" or',
     $  /7x,'"E-theta flag = off". A value of "on" causes inclusion',
     $  /7x,'of the higher-order electrostatic (E-theta) terms in',
     $  /7x,"Pitzer's",' equations. A value of "off" causes the',
     $  /7x,"exclusion of these terms. Most models based on Pitzer's",
     $  /7x,'equations are consistent with the inclusion of these',
     $  /7x,'terms, but some are not. The values of the mixture',
     $  /7x,'coefficients theta and psi are consistent with either',
     $  /7x,'the inclusion or the exclusion of these terms.')
        stop
      endif
      i = index(uline,'= on')
      if (i .gt. 0) then
        uethfl = uline(i + 2:i + 3)
      else
        i = index(uline,'= off')
        if (i .gt. 0) then
          uethfl = uline(i + 2:i + 4)
        else
          write (noutpt,1020) uline
          write (nttyo,1020) uline
 1020     format(/' * Error - (EQPT/rdpz3) The value of the E-theta',
     $    ' flag on the line:',
     $    //7x,'"',a32,'"',
     $    //7x,"isn't ",'"on" or "off", as is required.')
          stop
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     LOOP HERE for each species block
c
  100 continue
c
c     Skip past the block terminator line.
c
      read (ndat0s,1000,end=990,err=995) uline
      if (uline(1:8) .ne. uterm(1:8)) go to 100
c
c     Check for the end of the present superblock.
c
      read (ndat0s,1030,end=990,err=995) unam1,unam2,unam3
 1030 format(a24,2x,a24,2x,a24)
      if (unam1(1:8) .eq. uelemt(1:8)) then
c
c       Found the header for the elements block.
c
        npx3t = npx3
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Have another species triplet.
c
      npx3 = npx3 + 1
c
      j2 = ilnobl(unam1)
      j3 = ilnobl(unam2)
      j4 = ilnobl(unam3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for blank names.
c
      if (j2.le.0 .or. j3.le.0 .or. j4.le.0) then
        write (ux8,'(i5)') npx3
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
        if (j4 .le. 0) then
          unam3 = '<blank>'
          j4 = ilnobl(unam3)
        endif
        write (noutpt,1040) ux8(1:j5),unam1(1:j2),unam2(1:j3),
     $  unam3(1:j4)
        write (nttyo,1040) ux8(1:j5),unam1(1:j2),unam2(1:j3),
     $  unam3(1:j4)
 1040   format(/' * Error - (EQPT/rdpz3) Have encountered a blank',
     $  ' species name',/7x,'in the species triplet for block ',a,
     $  ' of the superblock for',/7x,'Pitzer parameter data',
     $  ' corresponding to species triplets.',/7x,'The offending',
     $  ' triplet shows as ',a,', ',a,', ',a,'.')
        if (npx3 .gt. 1) then
          ux24a = utripl(1,npx3 - 1)
          ux24b = utripl(2,npx3 - 1)
          ux24c = utripl(3,npx3 - 1)
          if (ux24a(1:7).ne.'<blank>' .or.
     $      ux24b(1:7).ne.'<blank>' .or.
     $      ux24c(1:7).ne.'<blank>') then
            j6 = ilnobl(ux24a)
            j7 = ilnobl(ux24b)
            j8 = ilnobl(ux24c)
            write (noutpt,1050) ux24a(1:j6),ux24b(1:j7),ux24c(1:j8)
            write (nttyo,1050) ux24a(1:j6),ux24b(1:j7),ux24c(1:j8)
 1050       format(7x,'This block follows the one for ',a,
     $      ', ',a,', ',a,'.')
          endif
        endif
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the indices of the three species (n1, n2, n3). These
c     indices point to the data read from the previously read
c     species blocks.
c
c     Get the index for the first species name.
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
 1100     format(/' * Error - (EQPT/rdpz3) The aqueous species ',
     $    a,' appearing in',/7x,'the Pitzer parameter data block',
     $    ' for the species triplet',/7x,a,', ',a,', ',a,' does',
     $    ' not appear in an aqueous',/7x,'species data block.')
          nerr = nerr + 1
        endif
      endif
c
c     Get the index for the second species name.
c
      qdup12 = .false.
      if (unam2(1:24) .eq. unam1(1:24)) then
c
        qdup12 = .true.
        n2 = n1
      else
c
c       Calling sequence substitutions:
c         n2 for na
c         unam2 for unams
c
        call gspidx(ier,n2,nat,natmax,uaqsp,unam2)
c
        if (ier .gt. 0) then
          if (unam1(1:7).ne.'<blank>' .and.
     $      unam2(1:7).ne.'<blank>' .and. unam3(1:7).ne.'<blank>') then
            write (noutpt,1100) unam2(1:j3),unam1(1:j2),unam2(1:j3),
     $      unam3(1:j4)
            write (nttyo,1100) unam2(1:j3),unam1(1:j2),unam2(1:j3),
     $      unam3(1:j4)
            nerr = nerr + 1
          endif
        endif
      endif
c
c     Get the index for the third species name.
c
      qdup13 = .false.
      qdup23 = .false.
      if (unam3(1:24) .eq. unam1(1:24)) then
c
        qdup13 = .true.
        if (qdup12) qdup23 = .true.
        n3 = n1
      elseif (unam3(1:24) .eq. unam2(1:24)) then
c
        qdup23 = .true.
        if (qdup12) qdup13 = .true.
        n3 = n2
      else
c
c       Calling sequence substitutions:
c         n3 for na
c         unam3 for unams
c
        call gspidx(ier,n3,nat,natmax,uaqsp,unam3)
c
        if (ier .gt. 0) then
          if (unam1(1:7).ne.'<blank>' .and.
     $      unam2(1:7).ne.'<blank>' .and. unam3(1:7).ne.'<blank>') then
            write (noutpt,1100) unam3(1:j4),unam1(1:j2),unam2(1:j3),
     $      unam3(1:j4)
            write (nttyo,1100) unam3(1:j4),unam1(1:j2),unam2(1:j3),
     $      unam3(1:j4)
            nerr = nerr + 1
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (n1.ne.0 .and. n2.ne.0 .and. n3.ne.0) then
c
c       Get the charges.
c
        z1 = zaqsp(n1)
        z2 = zaqsp(n2)
        z3 = zaqsp(n3)
c
        iz1 = nint(z1)
        iz2 = nint(z2)
        iz3 = nint(z3)
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
c       Check for illegal species combinations. These checks depend
c       on the charge combinations; hence, they should not be made
c       unless the all the charges have been determined.
c
        call tripck(na,nc,nerr,nn,noutpt,nttyo,n1,n2,n3,
     $  qdup12,qdup13,qdup23,unam1,unam2,unam3)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Copy the names before rearranging for storage.
c
      u1 = unam1
      u2 = unam2
      u3 = unam3
c
c     Rearrange for storage. The rules are:
c
c       If exactly one neutral is present:
c         neutral, cation, anion
c
c       If exactly two neutrals are present:
c         neutral 1, neutral 1, neutral2
c
c       If two cations are present:
c         cation 1, cation 2, anion
c
c       If two anions are present:
c         anion 1, anion 2, cation
c
c       If two species are of the same charge type,
c       store alphabetically.
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
c     Read the S-theta and psi parameters.
c
      read (ndat0s,1160,end=990,err=995) athetx(1),apsi(1,npx3)
 1160 format(13x,f8.5,13x,f8.5)
c
c     Read the first source.
c
      read (ndat0s,1170,end=990,err=995) uline
 1170 format(13x,a24)
c
c     Read the first and second temperature derivatives of the
c     S-theta and psi parameters.
c
      read (ndat0s,1180,end=990,err=995) athetx(2),athetx(3)
      read (ndat0s,1180,end=990,err=995) apsi(2,npx3),apsi(3,npx3)
 1180 format(13x,e10.3,13x,e10.3)
c
c     Read the second source.
c
      read (ndat0s,1170,end=990,err=995) uline
c
c     Check for illegal inputs in the S-theta data fields.
c
      call thetck(athetx,jpfcmx,na,nc,nerr,nn,noutpt,nttyo,
     $ n1,n2,n3,unam1,unam2,unam3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if ((nc.eq.2 .and. na.eq.1) .or.
     $  (nc.eq.1 .and. na.eq.2)) then
c
c       Find the theta index for the two cations or two anions.
c
        ifound = 0
        do nthd = 1,nthdt
          uthd1 = uthdtr(1,nthd)
          uthd2 = uthdtr(2,nthd)
          uthd3 = uthdtr(3,nthd)
          if (u1.eq.uthd1 .and. u2.eq.uthd2) then
            ifound = 1
            j2 = ilnobl(u1)
            j3 = ilnobl(u2)
            j4 = ilnobl(u3)
            j5 = ilnobl(uthd1)
            j6 = ilnobl(uthd2)
            j7 = ilnobl(uthd3)
            do j = 1,jpfcmx
              write (ux8,'(i5)') j
              call lejust(ux8)
              j8 = ilnobl(ux8)
              if (abs(atheta(j,nthd) - athetx(j)) .gt. 1.e-6) then
                write (noutpt,1200) ux8(1:j8),u1(1:j2),u2(1:j3),
     $          u3(1:j4),athetx(j),uthd1(1:j5),uthd2(1:j6),
     $          uthd3(1:j7),atheta(j,nthd)
                write (nttyo,1200) ux8(1:j8),u1(1:j2),u2(1:j3),
     $          u3(1:j4),athetx(j),uthd1(1:j5),uthd2(1:j6),
     $          uthd3(1:j7),atheta(j,nthd)
 1200           format(/' * Warning - (EQPT/rdpz3) The S-theta',
     $          ' a(',a,') coefficient for the',/7x,'triplet ',
     $          a,', ',a,', ',a,' is ',e12.5,',',/7x,'whereas the',
     $          ' value previously read for the triplet',
     $          /7x,a,', ',a,', ',a,' is ',e12.5,'. The former',
     $          ' value will be ignored.')
                nwarn = nwarn + 1
              endif
            enddo
            go to 110
          endif
        enddo
  110   continue
c
        if (ifound .eq. 0) then
c
c         Found a new pair defining a theta parameter.
c
          nthdt = nthdt + 1
          nthd = nthdt
          uthdtr(1,nthd) = u1
          uthdtr(2,nthd) = u2
          uthdtr(3,nthd) = u3
          do j = 1,jpfcmx
            atheta(j,nthd) = athetx(j)
          enddo
        endif
      endif
c
c     Go back to read the next species block.
c
      go to 100
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 write (noutpt,2000)
      write (nttyo,2000)
 2000 format(/' * Error - (EQPT/rdpz3) Unexpectedly encountered',
     $ /7x,'end-of-file while reading the DATA0 file.')
      stop
c
  995 write (noutpt,2010)
      write (nttyo,2010)
 2010 format(/' * Error - (EQPT/rdpz3) Encountered a read format',
     $ /7x,'error while reading the DATA0 file.')
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
c
c     Deallocate work array.
c
      DEALLOCATE(athetx)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
