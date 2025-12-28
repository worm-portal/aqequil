      subroutine gnenb(ipbt_asv,ikt_asv,jpdblo,jpfc_asv,nap_asv,
     $ nat_asv,nazt_asv,nbt_asv,nct_asv,ndat0s,ngt_asv,nlt_asv,
     4 nmt_asv,noutpt,npt_asv,npx2_asv,npx3_asv,nsb,nst_asv,nttyo,
     $ nxt_asv,uakey)
c
c     This subroutine makes a first pass through the DATA0 file to
c     determine the necessary dimensioning of arrays. The arrays are
c     then allocated back in the main program, and the data file is
c     then processed.
c
c     The returned array allocation size variables generally end
c     in "_asv." For example, the number of chemical elements is
c     normally nct. The associated allocation size variable is
c     nct_asv.
c
c     Note: although for example nxt_asv is described as the number
c     of solid solution phases on the data file, it is actually an
c     array allocation size variable which might have to be utilized
c     in allocating an array, even if the number of items in the set
c     is zero. Thus, the minimum value of one is generally imposed
c     on any allocation size variable.
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ipbt_asv = the number of parameters in a Pitzer alpha set
c       ndat0s   = unit number of the stripped DATA0 file
c
c     Principal output:
c
c       nap_asv  = maximum number of distinct sets of Pitzer alpha
c                    parameters
c       nazt_asv = the number of aqueous species on the data file
c                    for which hard core diameters are specified
c       nbt_asv  = the number of basis species on the data file
c       nct_asv  = the number of chemical elements on the data file
c       ikt_asv  = the maximum number of end-member component species
c                    in any solid solution on the data file
c       jpdblo   = integer flag denoting the Pitzer data block
c                    organization; -1 = classical, 0 = newlockile
c       jpfc_asv = the number of terms in the temperature function
c                    used to represent Pitzer interaction parameters
c       nat_asv  = the number of aqueous species on the data file
c       ngt_asv  = the number of gas species on the data file
c       nlt_asv  = the number of pure liquid species on the data file
c       nmt_asv  = the number of gas species on the data file
c       npt_asv  = the number of phases of all types on the data file
c       npx2_asv = the number of pairs of species not of the same
c                    charge sign for which Pitzer parameters are
c                    defined; typically, one of the pair is a cation
c                    and the other is an anion, but one or both
c                    species may also be electrically neutral
c       npx3_asv = the number of triplets of species corresponding to
c                    aqueous electrolyte mixtures for which Pitzer
c                    parameters are defined; generally, no more than
c                    two of these may have an electrical charge number
c                    that is postive, negative, or zero
c       nst_asv  = the number of species of all types on the data file
c       nxt_asv  = the number of solid-solution phases on the data file
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipbt_asv,ikt_asv,jpfc_asv,nap_asv,nat_asv,nazt_asv,
     $ nbt_asv,nct_asv,ngt_asv,nlt_asv,nmt_asv,npt_asv,npx2_asv,
     $ npx3_asv,nst_asv,nxt_asv
c
      integer ndat0s,noutpt,nttyo
c
      integer jpdblo,nsb
c
      character(len=8) uakey
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer, parameter :: nap_par = 500
c
      integer i,ier,j,ja,j2,j3,j4,j5,ikt_sum,ikt_x,n,nmax
c
      integer ilnobl
c
      logical qalpha
c
      character(len=80) udastr,uline,ux80
      character(len=16) uaqu,uaux,uele,ubas,ubdot,ugas,uliq,umin,
     $ uref,usso,uterm
      character(len=16) ustr16
      character(len=8) ux8
c
      real(8) palpha(ipbt_asv,nap_par),palphi(ipbt_asv)
c
      real(8) aax,var
c
c-----------------------------------------------------------------------
c
      data uaqu   /'aqueous species '/
      data uaux   /'auxiliary basis '/
      data ubas   /'basis species   '/
      data ubdot  /'bdot parameters '/
      data uele   /'elements        '/
      data ugas   /'gases           '/
      data uliq   /'liquids         '/
      data umin   /'solids          '/
      data uref   /'references      '/
      data usso   /'solid solutions '/
      data uterm  /'+---------------'/
c
c-----------------------------------------------------------------------
c
c     Initialize some variables to zero.
c
      ikt_asv = 0
      nazt_asv = 0
      nap_asv = 0
      npx2_asv = 0
      npx3_asv = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (uakey(1:8) .eq. 'SEDH    ') then
c
c       Have a Simple Extended Debye-Huckel model for the activity
c       coefficients of aqueous species. Determine the number of aqueous
c       species for which hard core diameters are specified.
c
  120   read (ndat0s,1000,end=990,err=995) uline
 1000   format(a)
        if (uline(1:16) .ne. ubdot(1:16)) go to 120
c
c       Skip the terminator line.
c
        read (ndat0s,1000,end=990,err=995) uline
c
c       Tally lines until next terminator.
c
        n = 0
  130   read (ndat0s,1000,end=990,err=995) uline
        if (uline(1:16) .eq. uterm(1:16)) go to 140
        n = n + 1
        go to 130
c
  140   nazt_asv = n
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (uakey(1:8) .eq. 'Pitzer  ') then
c
c       Have a Pitzer model for the activity coefficients of aqueous
c       species.
c
        if (jpdblo .eq. -1) then
c
c         Have the "classical" Pitzer data block organization.
c         Determine the number of pairs of ions corresponding to pure
c         aqueous electrolytes for which Pitzer parameters are defined
c         ('single-salt parameters').
c
  150     read (ndat0s,1000,end=990,err=995) uline
          if (uline(1:16) .ne. 'single-salt para') go to 150
c
c         Skip the terminator line.
c
          read (ndat0s,1000,end=990,err=995) uline
c
c         Determine the number of such pairs of ions (npx2_asv) by
c         counting the number of terminator lines prior to encountering
c         the 'mixture term parameters' line.
c
c         Also determine the number of distinct sets of Pitzer alpha
c         parameters (nap_asv).
c
          n = 0
          ja = 0
  160     read (ndat0s,1000,end=990,err=995) uline
c
          if (uline(1:16) .eq. 'mixture term par') then
            backspace(ndat0s)
            go to 180
          endif
c
          i = index(uline,'alpha1')
          if (i .gt. 0) then
c
c           Found the line with the alpha parameters.
c
            read (uline,1020,end=990,err=995) palphi(1),palphi(2)
 1020       format(18x,2(16x,f5.1))
c
            if (ja .gt. 1) then
              do j = 1,ja
                if (abs(palpha(1,j) - palphi(1)).le.1.e-12
     $            .and. abs(palpha(2,j) - palphi(2)).le.1.e-12) then
                  go to 170
                endif
              enddo
            endif
c
c           Not in the set; add it.
c
            ja = ja + 1
            if (ja .gt. nap_par) then
              write (ux8,'(i5)') nap_par
              call lejust(ux8)
              j2 = ilnobl(ux8)
              write (noutpt,1030) ux8(1:j2)
              write (nttyo,1030) ux8(1:j2)
 1030         format(/' * Error - (EQPT/gnenb) Have overflowed the',
     $        ' palpha array while',/7x,'attempting to find the number',
     $        ' of distinct sets of Pitzer alpha',/7x,'parameters.',
     $        ' Increase the dimensioning parameter nap_par in this',
     $        /7x,'subroutine from its present value of ',a,'.')
              stop
            endif
c
            palpha(1,ja) = palphi(1)
            palpha(2,ja) = palphi(2)
          endif
  170     continue
c
          if (uline(1:16) .eq. uterm(1:16)) n = n + 1
          go to 160
c
  180     npx2_asv = n
          nap_asv = ja
c
c         Determine the number of triplets of ions corresponding to
c         aqueous electrolyte mixtures for which Pitzer parameters are
c         defined ('mixture term parameters').
c
  190     read (ndat0s,1000,end=990,err=995) uline
          if (uline(1:16) .ne. 'mixture term par') go to 190
c
c         Skip the terminator line.
c
          read (ndat0s,1000,end=990,err=995) uline
c
c         Skip the E-theta flag line.
c
          read (ndat0s,1000,end=990,err=995) uline
c
c         Determine the number such triplets of by counting the number
c         of terminator lines prior to encountering the 'elements' line.
c
          n = 0
  200     read (ndat0s,1000,end=990,err=995) uline
          if (uline(1:16) .eq. uele(1:16)) then
            backspace(ndat0s)
            go to 210
          endif
          if (uline(1:16) .eq. uterm(1:16)) n = n + 1
          go to 200
c
  210     npx3_asv = n
c
        elseif (jpdblo .eq. 0) then
c
c         Have the "new" Pitzer data block organization.
c
c         First determine the number of distinct sets of Pitzer alpha
c         parameters (nap_asv). It is possible that no values are
c         input here, in which case only the default set holds.
c
c         Load the default set into the palpha array.
c
          nmax = ipbt_asv*nap_par
          call initaz(palpha,nmax)
          palpha(1,1) = 2.0
          palpha(2,1) = 12.0
          palpha(1,2) = 1.4
          palpha(2,2) = 12.0
          ja = 2
          nap_asv = ja
c
c         Skip down to the "ca combinations" superblock. This is
c
  220     read (ndat0s,1000,end=990,err=995) uline
          if (uline(1:15) .ne. 'ca combinations') go to 220
c
c         Skip the terminator line.
c
          read (ndat0s,1000,end=990,err=995) uline
c
  230     read (ndat0s,1000,end=990,err=995) uline
          j3 = ilnobl(uline)
c
          if (uline(1:26) .eq. "cc'a and aa'c combinations") then
            backspace(ndat0s)
            go to 250
          endif
c
          qalpha = .false.
          ux80 = uline
          call lejust(ux80)
c
c         Find the next set of alpha lines. Look for "alpha(1)".
c         That marks the start of a set.
c
          ustr16 = 'alpha(1) ='
          j4 = 10
          if (ux80(1:j4) .eq. ustr16(1:j4)) then
c
c           Found a line with the alpha(1) parameter.
c
            j5 = index(ux80,'=')
            udastr = ux80(j5 + 1:80)
            call g1dat(ier,noutpt,nttyo,udastr,var)
            if (ier .gt. 0) then
              write (noutpt,1110) uline(1:j3)
              write (nttyo,1110) uline(1:j3)
 1110         format(/' * Error - (EQPT/gnenb) Have found a line',
     $        ' starting with:',/7x,'"',a,'"',/7x,'containing an',
     $        ' expected numerical field that could not be read.',
     $        /7x,'This occurred while scanning the data file to',
     $        ' determine the number',/7x,'of distinct sets of',
     $        ' Pitzer alpha parameters.')
              stop
            endif
            palphi(1) = var
            qalpha = .true.
c
            do i = 2,ipbt_asv
              read (ndat0s,1000,end=990,err=995) uline
              j3 = ilnobl(uline)
              ux80 = uline
              call lejust(ux80)
              ustr16 = 'alpha( ) ='
              j4 = 10
              write (ustr16(7:7),'(i1)') i
              if (ux80(1:j4) .eq. ustr16(1:j4)) then
c
c               Found a line with the expected alpha(i) parameter.
c
                j5 = index(ux80,'=')
                udastr = ux80(j5 + 1:80)
                call g1dat(ier,noutpt,nttyo,udastr,var)
                if (ier .le. 0) then
                  palphi(i) = var
                else
                  write (noutpt,1110) uline(1:j3)
                  write (nttyo,1110) uline(1:j3)
                  stop
                endif
              else
c
c               Did not find a line with the expected alpha(i)
c               parameter.
c
                ustr16 = 'beta(0) ='
                j4 = 9
                if (ux80(1:j4) .ne. ustr16(1:j4)) then
c
c                 Did not find a beta(0) line either. This
c                 condition is invalid.
c
                  write (ux8,'(i1)') i
                  write (noutpt,1120) uline(1:j3),ux8(1:1)
                  write (nttyo,1120) uline(1:j3),ux8(1:1)
 1120             format(/' * Error - (EQPT/gnenb) Have found a line',
     $            ' starting with:',/7x,'"',a,'"',/7x,'that contains',
     $            ' neither the expected alpha(',a,') data nor the',
     $            /7x,'possible beta(0) data that could mark the end',
     $            ' of an abbreviated',/7x,'alpha set. This occurred',
     $            ' while scanning the data file',/7x,'to determine',
     $            ' the number of distinct sets of Pitzer',
     $            /7x,'alpha parameters.')
                  stop
                endif
              endif
            enddo
          endif
c
          if (qalpha) then
            do j = 1,ja
              do i = 1,ipbt_asv
                aax = abs(palpha(i,j) - palphi(i))
                if (aax .gt. 1.e-12) go to 240
              enddo
c
c             The current alpha combination is already in
c             the known set.
c
              go to 230
  240         continue
            enddo
c
c           The current alpha set is not in the set; add it.
c
            ja = ja + 1
            if (ja .gt. nap_par) then
c
c             Have overflowed the palpha array, which was statically
c             dimensioned.
c
              write (ux8,'(i5)') nap_par
              call lejust(ux8)
              j2 = ilnobl(ux8)
              write (noutpt,1030) ux8(1:j2)
              write (nttyo,1030) ux8(1:j2)
              stop
            endif
            do i = 1,ipbt_asv
              palpha(i,ja) = palphi(i)
            enddo
          endif
          go to 230
c
  250     nap_asv = ja
c
          rewind(ndat0s)
c
c         Skip back down to the "ca combinations" superblock.
c
  260     read (ndat0s,1000,end=990,err=995) uline
          if (uline(1:15) .ne. 'ca combinations') go to 260
c
c         Skip the terminator line.
c
          read (ndat0s,1000,end=990,err=995) uline
c
c         Determine the number of pair combinations (npx2_asv) by
c         counting the number of terminator lines prior to
c         encountering the "cc'a and aa'c combinations" line.
c         Correct this estimate for other pair superblock
c         headers.
c
          n = 0
  270     read (ndat0s,1000,end=990,err=995) uline
c
          if (uline(1:26) .eq. "cc'a and aa'c combinations") then
            backspace(ndat0s)
            go to 280
          endif
c
          if (uline(1:24) .eq. "cc' and aa' combinations") n = n - 1
          if (uline(1:22) .eq. "nc and na combinations") n = n - 1
          if (uline(1:15) .eq. "nn combinations") n = n - 1
          if (uline(1:16) .eq. "nn' combinations") n = n - 1
          if (uline(1:16) .eq. uterm(1:16)) n = n + 1
          go to 270
c
  280     npx2_asv = n
c
c         Skip down to the "cc'a and aa'c combinations" superblock.
c
  290     read (ndat0s,1000,end=990,err=995) uline
          if (uline(1:26) .ne. "cc'a and aa'c combinations") go to 290
c
c         Skip the terminator line.
c
          read (ndat0s,1000,end=990,err=995) uline
c
c         Determine the number of triplet combinations (npx3_asv) by
c         counting the number of terminator lines prior to
c         encountering the "cc'a and aa'c combinations" line.
c         Correct this estimate for other triplet superblock
c         headers.
c
          n = 0
  300     read (ndat0s,1000,end=990,err=995) uline
c
          if (uline(1:16) .eq. uele(1:16)) then
            backspace(ndat0s)
            go to 310
          endif
c
          if (uline(1:16) .eq. "nca combinations") n = n - 1
          if (uline(1:16) .eq. "nnn' combination") n = n - 1
          if (uline(1:16) .eq. uterm(1:16)) n = n + 1
          go to 300
c
  310     npx3_asv = n
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find the 'elements' line.
c
  350 read (ndat0s,1000,end=990,err=995) uline
      if (uline(1:16) .ne. uele(1:16)) go to 350
c
c     Skip the terminator line.
c
      read (ndat0s,1000,end=990,err=995) uline
c
c     Tally lines until next terminator.
c
      n = 0
  360 read (ndat0s,1000,end=990,err=995) uline
      if (uline(1:16) .eq. uterm(1:16)) go to 370
      n = n + 1
      go to 360
c
  370 nct_asv = n
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Skip to the 'basis species' line.
c
  400 read (ndat0s,1000,end=990,err=995) uline
      if (uline(1:16) .ne. ubas(1:16)) go to 400
c
c     Skip the terminator line.
c
      read (ndat0s,1000,end=990,err=995) uline
c
c     Determine the number of basis species by counting the number
c     of terminator lines prior to encountering the 'aqueous species'
c     line. In this count, note that an "extra" terminator line
c     is encountered after the "auxiliary basis species" line.
c     The number of strict basis species is also determined, by
c     noting the number of basis species prior to encountering the
c     "auxiliary basis species" line.
c
      n = 0
  410 read (ndat0s,1000,end=990,err=995) uline
      if (uline(1:16) .eq. uaux(1:16)) nsb = n
      if (uline(1:16) .eq. uaqu(1:16)) then
        backspace(ndat0s)
        go to 420
      endif
      if (uline(1:16) .eq. uterm(1:16)) n = n + 1
      go to 410
c
  420 nbt_asv = n - 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Skip to the 'aqueous species' line.
c
  450 read (ndat0s,1000,end=990,err=995) uline
      if (uline(1:16) .ne. uaqu(1:16)) go to 450
c
c     Skip the terminator line.
c
      read (ndat0s,1000,end=990,err=995) uline
c
c     Determine the number of non-basis aqueous species by counting the
c     number of terminator lines prior to encountering the 'solids'
c     line. The total number of aqueous species is the sum of the
c     basis and non-basis aqueous species.
c
      n = 0
  460 read (ndat0s,1000,end=990,err=995) uline
      if (uline(1:16) .eq. umin(1:16)) then
        backspace(ndat0s)
        go to 470
      endif
      if (uline(1:16) .eq. uterm(1:16)) n = n + 1
      go to 460
c
  470 nat_asv = nbt_asv + n
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Skip to the 'solids' line.
c
  500 read (ndat0s,1000,end=990,err=995) uline
      if (uline(1:16) .ne. umin(1:16)) go to 500
c
c     Skip the terminator line.
c
      read (ndat0s,1000,end=990,err=995) uline
c
c     Determine the number of pure mineral species by counting the
c     number of terminator lines prior to encountering the 'liquids'
c     line.
c
      n = 0
  510 read (ndat0s,1000,end=990,err=995) uline
      if (uline(1:16) .eq. uliq(1:16)) then
        backspace(ndat0s)
        go to 520
      endif
      if (uline(1:16) .eq. uterm(1:16)) n = n + 1
      go to 510
c
  520 nmt_asv = n
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Skip to the 'liquids' line.
c
  550 read (ndat0s,1000,end=990,err=995) uline
      if (uline(1:16) .ne. uliq(1:16)) go to 550
c
c     Skip the terminator line.
c
      read (ndat0s,1000,end=990,err=995) uline
c
c     Determine the number of pure liquid species by counting the
c     number of terminator lines prior to encountering the 'gases'
c     line.
c
      n = 0
  560 read (ndat0s,1000,end=990,err=995) uline
      if (uline(1:16) .eq. ugas(1:16)) then
        backspace(ndat0s)
        go to 570
      endif
      if (uline(1:16) .eq. uterm(1:16)) n = n + 1
      go to 560
c
c     Note: water as a pure liquid phase is not explicitly on the
c     DATA0 file. It is added by EQ3NR or EQ6 by cloning the liquid
c     water species belonging to the aqueous solution phase. The pure
c     liquid water phase must be accounted for in nlt_asv, so the
c     actual count of pure liquid phases explicitly on the data file
c     is incremented by one.
c
  570 nlt_asv = n + 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Skip to the 'gases' line.
c
  600 read (ndat0s,1000,end=990,err=995) uline
      if (uline(1:16) .ne. ugas(1:16)) go to 600
c
c     Skip the terminator line.
c
      read (ndat0s,1000,end=990,err=995) uline
c
c     Determine the number of gas species by counting the number of
c     terminator lines prior to encountering the 'solid solutions'
c     line.
c
      n = 0
  610 read (ndat0s,1000,end=990,err=995) uline
      if (uline(1:16) .eq. usso(1:16)) then
        backspace(ndat0s)
        go to 620
      endif
      if (uline(1:16) .eq. uterm(1:16)) n = n + 1
      go to 610
c
  620 ngt_asv = n
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Skip to the 'solid solutions' line.
c
  650 read (ndat0s,1000,end=990,err=995) uline
      if (uline(1:16) .ne. usso(1:16)) go to 650
c
c     Skip the terminator line.
c
      read (ndat0s,1000,end=990,err=995) uline
c
c     Determine the number of solid solution phases by counting the
c     number of terminator lines prior to encountering the 'references'
c     line.
c
c       ikt_x   = the number of components in the current solid
c                   solution
c       ikt_sum = the total number of components in all solid solutions
c
      n = 0
      ikt_x = 0
      ikt_sum = 0
c
  660 read (ndat0s,1000,end=990,err=995) uline
      if (uline(1:16) .eq. uref(1:16)) then
        backspace(ndat0s)
        go to 670
      endif
c
      if (ikt_x .le. 0) then
c
c       Have not yet found a 'component' line.
c
        i = index(uline,'component')
        if (i .eq. 5) then
c
c         Have found the line with the number of components.
c         Note that there is a two-component minimum here.
c
          read (uline,'(1x,i2)',err=995) ikt_x
          ikt_x = max(2,ikt_x)
          ikt_asv = max(ikt_asv,ikt_x)
          ikt_sum = ikt_sum + ikt_x
        endif
      endif
c
      if (uline(1:16) .eq. uterm(1:16)) then
c
c       Have found the end of the current solid solution block.
c
        n = n + 1
        if (ikt_x .le. 0) then
c
c         No line defining the number of components was found.
c
          ikt_x = 2
          ikt_asv = max(ikt_asv,ikt_x)
          ikt_sum = ikt_sum + ikt_x
        endif
        ikt_x = 0
      endif
      go to 660
c
  670 nxt_asv = n
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the total number of phases. Be sure to account for the
c     aqueous solution and the gas phase. Water as a pure liquid
c     phase has been implicitly added to nlt_asv.
c
      npt_asv = 2 + nmt_asv + nlt_asv + nxt_asv
c
c     Get the total number of species.
c
      nst_asv = nat_asv + nmt_asv + nlt_asv + ngt_asv + ikt_sum
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Each dimensioning parameter should have a minimum value of 1,
c     even if the corresponding number of items is 0.
c
      ikt_asv = max(1,ikt_asv)
      nap_asv = max(1,nap_asv)
      nat_asv = max(1,nat_asv)
      nazt_asv = max(1,nazt_asv)
      nbt_asv = max(1,nbt_asv)
      nct_asv = max(1,nct_asv)
      ngt_asv = max(1,ngt_asv)
      nlt_asv = max(1,nlt_asv)
      nmt_asv = max(1,nmt_asv)
      npt_asv = max(1,npt_asv)
      npx2_asv = max(1,npx2_asv)
      npx3_asv = max(1,npx3_asv)
      nst_asv = max(1,nst_asv)
      nxt_asv = max(1,nxt_asv)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Rewind the DATA0 file and exit.
c
      rewind(ndat0s)
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write a message for any read error.
c
  990 write (noutpt,2000)
      write (nttyo,2000)
 2000 format(/' * Error - (EQPT/gnenb) Unexpectedly encountered',
     $ ' end-of-file',/7x,'while scanning the DATA0 file to determine',
     $ ' the necessary array',/7x,'dimensions.')
      stop
c
  995 write (noutpt,2010)
      write (nttyo,2010)
 2010 format(/' * Error - (EQPT/gnenb) Encountered a read format',
     $ ' error while',/7x,'scanning the DATA0 file to determine',
     $ ' the necessary array dimensions.')
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
