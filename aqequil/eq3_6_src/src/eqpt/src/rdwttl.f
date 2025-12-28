      subroutine rdwttl(ipch,ipcv,jpdblo,jpfcmx,jptffl,narxt,ndata1,
     $ ndat0s,ndat1f,noutpt,nslist,ntitld,ntidmx,ntprmx,ntprt,nttyo,
     $ uakey,utitld)
c
c     This suboutine reads the title on the DATA0 file and writes it
c     on the DATA1 and DATA1F files. It checks the DATA0 file
c     header for validity. It also searches the title for embedded
c     flags and data indicating special treatment.
c
c     Possible embedded flags and data deal with the following:
c
c       1. Defining the temperature grid used on the data file.
c          By default, the classic eight-temperature grid is used.
c       2. Indicating the presence of data grids for enthalpy
c          and volume functions. By default, these data grids
c          are assumed to be absent.
c       3. Defining the temperature function used to describe
c          Pitzer interaction parameters. By default, this function
c          is the classic second-order Taylor's series centered
c          at 25C.
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
c       ipch   = enthalpy functions data grid flag:
c                  -1 = no enthalpy grids
c                   0 = enthalpy grids present
c                   1 = grids present for the enthalpy and its first
c                         pressure derivative
c                   2 = grids present for the enthalpy and its first
c                         and second pressure derivatives
c       ipcv   = volume functions data flag:
c                  -1 = no volume grids
c                   0 = volume grids present
c                   1 = grids present for the volume and its first
c                         pressure derivative
c                   2 = grids present for the volume and its first
c                         and second pressure derivatives
c       ntitld = the number of lines in the title
c       utitld = array of title lines
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer jpfcmx,ntidmx,ntprmx
c
      integer narxt(ntprmx)
c
      integer ndata1,ndat0s,ndat1f,noutpt,nslist,ntitld,nttyo
c
      integer ipch,ipcv,jpdblo,jptffl,ntprt
c
      character(len=80) utitld(ntidmx)
      character(len=8) uakey
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j,j2,j3,j4,k,n,narxti,nerr,ntpr,ntpri
c
      integer ilnobl
c
      logical qrderr
c
      character(len=56) ux56
c
      character(len=80) ulbufa,ulbufb
      character(len=72) uterm,utermc
      character(len=32) upchst,upcvst
      character(len=8) ustr,ux8a,ux8b
c
c-----------------------------------------------------------------------
c
      uterm(1:48) = '+-----------------------------------------------'
      uterm(49:72) = '------------------------'
      utermc = uterm
      utermc(1:1) = '*'
c
c     Read the data file title.
c
      ntitld = ntidmx
      do n = 1,ntitld
        read (ndat0s,1000,end=990,err=995) utitld(n)
 1000   format(a)
        if (utitld(n)(1:8) .eq. uterm(1:8)) go to 110
      enddo
  110 continue
c
c     Write title information.
c
      do n = 1,ntitld
         j2 = ilnobl(utitld(n))
         j2 = min(j2,79)
         write (ndata1) utitld(n)
         write (ndat1f,1000) utitld(n)(1:j2)
         write (nttyo,1140) utitld(n)(1:j2)
         write (noutpt,1140) utitld(n)(1:j2)
         write (nslist,1140) utitld(n)(1:j2)
 1140    format(' ',a)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the parameters for the temperature grid used for log K values
c     and standard state enthalpy and volume functions.
c
c       ntprt       = the number of temperature ranges
c       narxt(ntpr) = the number of points or coefficients in the
c                       ntpr-th temperature range
c
c     The last point of the ntpr-th temperature range is also the first
c     point of the following range, if any. This insures continuity.
c     The default values correspond to the original EQ3/6 data grid:
c
c       ntprt = 2
c       narxt(1) = 4  (0, 25, 60 100 C)
c       narxt(2) = 5  (100, 150, 200, 250, 300 C)
c
c     Note: the number of temperature ranges (ntprt) is equal to
c     the dimensioned limit (ntpr_asv), which has been determined
c     by a prior scan. The prior scan has also provided narx_asv,
c     the dimenionsed limit on the number of values in a temperature
c     range. This was found as the greatest value of any range.
c     Now, it is necessary to get the actual value for each range.
c
      ntprt = ntprmx
      call initiz(narxt,ntprmx)
      if (ntprt .eq. 2) then
        narxt(1) = 4
        narxt(2) = 5
      endif
c
      do ntpr = 1,ntprt
        do n = 1,ntitld
          ulbufa = ' '
          ulbufb = ' '
          j = 0
c
c         Check for a keystring.
c
          i = index(utitld(n),'NO. OF POINTS IN RANGE')
          if (i .gt. 0) then
            j = i + 22
          else
             i = index(utitld(n),'NUMBER OF POINTS IN RANGE')
            if (i .gt. 0) then
              j = i + 25
            else
              i = index(utitld(n),'NARXT')
              if (i .gt. 0) j = i + 5
            endif
          endif
c
          if (j .gt. 0) then
c
c           Extract a number acting as a subscript to the keystring.
c
            ulbufb = utitld(n)(j:80)
            call lejust(ulbufb)
            j = index(ulbufb,' ')
            j = min(j,9)
            k = j - 1
            if (k .gt. 0) then
              ustr = ulbufb(1:k)
              call chrint(ntpri,nttyo,qrderr,ustr)
              if (qrderr) go to 997
              if (ntpri .ne. ntpr) j = 0
            else
              j = 0
            endif
          endif
c
          if (j .gt. 0) then
c
c           Check for an equal sign following the keystring.
c
            ulbufa = ulbufb(j:80)
            call lejust(ulbufa)
            i = index(ulbufa,'=')
            j = 2
            if (i .ne. 1) j = 0
          endif
c
          if (j .gt. 0) then
c
c           Extract a number matching the keystring.
c
            ulbufb = ulbufa(j:80)
            call lejust(ulbufb)
            ustr = ulbufb(1:8)
            call chrint(narxti,nttyo,qrderr,ustr)
            if (qrderr) go to 997
            narxt(ntpr) = narxti
            go to 220
          endif
        enddo
  220   continue
c
      enddo
c
      nerr = 0
      ntpr = 1
      if (narxt(ntpr) .lt. 2) then
        write (ux8a,'(i5)') ntpr
        call lejust(ux8a)
        j3 = ilnobl(ux8a)
        write (ux8b,'(i5)') narxt(ntpr)
        call lejust(ux8b)
        j4 = ilnobl(ux8b)
        write (noutpt,1320) ux8a(1:j3),ux8b(1:j4)
        write (nttyo,1320) ux8a(1:j3),ux8b(1:j4)
 1320   format(/' * Error - (EQPT/rdwttl) The number of points in',
     $  /7x,'temperature range ',a,' is ',a,'. The number of',
     $  /7x,'points in the first or last temperature range must be',
     $  /7x,'at least 2.')
        nerr = nerr + 1
      endif
c
      do ntpr = 2,ntprt - 1
        if (narxt(ntpr) .lt. 3) then
          write (ux8a,'(i5)') ntpr
          call lejust(ux8a)
          j3 = ilnobl(ux8a)
          write (ux8b,'(i5)') narxt(ntpr)
          call lejust(ux8b)
          j4 = ilnobl(ux8b)
          write (noutpt,1330) ux8a(1:j3),ux8b(1:j4)
          write (nttyo,1330) ux8a(1:j3),ux8b(1:j4)
 1330     format(/' * Error - (EQPT/rdwttl) The number of points in',
     $    /7x,'temperature range ',a,' is ',a,'. The number of',
     $    /7x,'points in an interior temperature range must be at',
     $    /7x,'least 3.')
          nerr = nerr + 1
        endif
      enddo
c
      ntpr = ntprt
      if (ntpr.gt.1 .and. narxt(ntpr).lt.2) then
        write (ux8a,'(i5)') ntpr
        call lejust(ux8a)
        j3 = ilnobl(ux8a)
        write (ux8b,'(i5)') narxt(ntpr)
        call lejust(ux8b)
        j4 = ilnobl(ux8b)
        write (noutpt,1320) ux8a(1:j3),ux8b(1:j4)
        write (nttyo,1320) ux8a(1:j3),ux8b(1:j4)
        nerr = nerr + 1
      endif
c
      if (nerr .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the grid parameters.
c
      ux56 = 'Number of ranges in the logK temperature grid'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1350) ux56(1:j2)
 1350 format(a)
c
      write (ndata1) ntprt
      write (ndat1f,1360) ntprt
 1360 format(i5)
c
      do ntpr = 1,ntprt
        ux56 = 'Number of points in range '
        j2 = ilnobl(ux56)
        write (ux8a,'(i5)') ntpr
        call lejust(ux8a)
        j3 = ilnobl(ux8a)
        ux56(j2 + 2: j2 + 1 + j3) = ux8a(1:j3)
c
        write (ndata1) ux56
        j2 = ilnobl(ux56)
        write (ndat1f,1350) ux56(1:j2)
c
        n = narxt(ntpr)
        write (ndata1) n
        write (ndat1f,1360) n
      enddo
c
      write (ux8a,'(i5)') ntprt
      call lejust(ux8a)
      j3 = ilnobl(ux8a)
      write (noutpt,1380) ux8a(1:j3)
      write (nttyo,1380) ux8a(1:j3)
 1380 format(/' Number of logK temperature grid ranges= ',a)
      do ntpr = 1,ntprt
        write (ux8a,'(i5)') ntpr
        call lejust(ux8a)
        j2 = ilnobl(ux8a)
        n = narxt(ntpr)
        write (ux8b,'(i5)') n
        call lejust(ux8a)
        j3 = ilnobl(ux8a)
 1390   format('   Number of points in range ',a,'= ',a)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set flags for enthalpy and volume functions data grids.
c
      ipch = -1
      do n = 1,ntitld
        ulbufa = ' '
        j = 0
c
c       Check for a keystring.
c
        i = index(utitld(n),'ENTHALPY')
        if (i .gt. 0) j = i + 8
c
        if (j .gt. 0) then
c
c         Check for an equal sign following the keystring.
c
          ulbufa = utitld(n)(j:80)
          call lejust(ulbufa)
          i = index(ulbufa,'=')
          j = 2
          if (i .ne. 1) j = 0
        endif
c
        if (j .gt. 0) then
c
c         Extract a string input matching the keystring.
c
          ustr = ulbufa(2:9)
          call lejust(ustr)
          k = index(ustr,'ON')
          if (k .le. 0) k = index(ustr,'PRESENT')
          if (k .le. 0) k = index(ustr,'ACTIVE')
          if (k .gt. 0) ipch = 0
          go to 300
        endif
      enddo
  300 continue
c
      do n = 1,ntitld
        ulbufa = ' '
        j = 0
c
c       Check for a keystring.
c
        i = index(utitld(n),'dH/dP')
        if (i .gt. 0) j = i + 5
c
        if (j .gt. 0) then
c
c         Check for an equal sign following the keystring.
c
          ulbufa = utitld(n)(j:80)
          call lejust(ulbufa)
          i = index(ulbufa,'=')
          j = 2
          if (i .ne. 1) j = 0
        endif
c
        if (j .gt. 0) then
c
c         Extract a string input matching the keystring.
c
          ustr = ulbufa(2:9)
          call lejust(ustr)
          k = index(ustr,'ON')
          if (k .le. 0) k = index(ustr,'PRESENT')
          if (k .le. 0) k = index(ustr,'ACTIVE')
          if (k .gt. 0) then
            if (ipch .eq. 0) then
              ipch = 1
            else
              write (noutpt,1400)
              write (nttyo,1400)
 1400         format(/" * Error - (EQPT/rdwttl) Can't use dH/dP data",
     $        ' without',/7x,'enthalpy (H) data.')
              stop
            endif
          endif
          go to 310
        endif
      enddo
  310 continue
c
      do n = 1,ntitld
        ulbufa = ' '
        j = 0
c
c       Check for a keystring.
c
        i = index(utitld(n),'d2H/dP2')
        if (i .gt. 0) j = i + 7
c
        if (j .gt. 0) then
c
c         Check for an equal sign following the keystring.
c
          ulbufa = utitld(n)(j:80)
          call lejust(ulbufa)
          i = index(ulbufa,'=')
          j = 2
          if (i .ne. 1) j = 0
        endif
c
        if (j .gt. 0) then
c
c         Extract a string input matching the keystring.
c
          ustr = ulbufa(2:9)
          call lejust(ustr)
          k = index(ustr,'ON')
          if (k .le. 0) k = index(ustr,'PRESENT')
          if (k .le. 0) k = index(ustr,'ACTIVE')
          if (k .gt. 0) then
            if (ipch .eq. 1) then
              ipch = 2
            else
              write (noutpt,1410)
              write (nttyo,1410)
 1410         format(/" * Error - (EQPT/rdwttl) Can't use d2H/dP2",
     $        ' data without',/7x,'enthalpy (H) and dH/dP data.')
              stop
            endif
          endif
          go to 320
        endif
      enddo
  320 continue
c
      ipcv = -1
      do n = 1,ntitld
        ulbufa = ' '
        j = 0
c
c       Check for a keystring.
c
        i = index(utitld(n),'VOLUME')
        if (i .gt. 0) j = i + 6
c
        if (j .gt. 0) then
c
c         Check for an equal sign following the keystring.
c
          ulbufa = utitld(n)(j:80)
          call lejust(ulbufa)
          i = index(ulbufa,'=')
          j = 2
          if (i .ne. 1) j = 0
        endif
c
        if (j .gt. 0) then
c
c         Extract a string input matching the keystring.
c
          ustr = ulbufa(2:9)
          call lejust(ustr)
          k = index(ustr,'ON')
          if (k .le. 0) k = index(ustr,'PRESENT')
          if (k .le. 0) k = index(ustr,'ACTIVE')
          if (k .gt. 0) ipcv = 0
          go to 330
        endif
      enddo
  330 continue
c
      do n = 1,ntitld
        ulbufa = ' '
        j = 0
c
c       Check for a keystring.
c
        i = index(utitld(n),'dV/dP')
        if (i .gt. 0) j = i + 5
c
        if (j .gt. 0) then
c
c         Check for an equal sign following the keystring.
c
          ulbufa = utitld(n)(j:80)
          call lejust(ulbufa)
          i = index(ulbufa,'=')
          j = 2
          if (i .ne. 1) j = 0
        endif
c
        if (j .gt. 0) then
c
c         Extract a string input matching the keystring.
c
          ustr = ulbufa(2:9)
          call lejust(ustr)
          k = index(ustr,'ON')
          if (k .le. 0) k = index(ustr,'PRESENT')
          if (k .le. 0) k = index(ustr,'ACTIVE')
          if (k .gt. 0) then
            if (ipcv .eq. 0) then
              ipcv = 1
            else
              write (noutpt,1420)
              write (nttyo,1420)
 1420         format(/" * Error - (EQPT/rdwttl) Can't use dV/dP data",
     $        ' without',/7x,'volume (V) data.')
              stop
            endif
          endif
          go to 340
        endif
      enddo
  340 continue
c
      do n = 1,ntitld
        ulbufa = ' '
        j = 0
c
c       Check for a keystring.
c
        i = index(utitld(n),'d2V/dP2')
        if (i .gt. 0) j = i + 7
c
        if (j .gt. 0) then
c
c         Check for an equal sign following the keystring.
c
          ulbufa = utitld(n)(j:80)
          call lejust(ulbufa)
          i = index(ulbufa,'=')
          j = 2
          if (i .ne. 1) j = 0
        endif
c
        if (j .gt. 0) then
c
c         Extract a string input matching the keystring.
c
          ustr = ulbufa(2:9)
          call lejust(ustr)
          k = index(ustr,'ON')
          if (k .le. 0) k = index(ustr,'PRESENT')
          if (k .le. 0) k = index(ustr,'ACTIVE')
          if (k .gt. 0) then
            if (ipcv .eq. 1) then
              ipcv = 2
            else
              write (noutpt,1430)
              write (nttyo,1430)
 1430         format(/" * Error - (EQPT/rdwttl) Can't use d2V/dP2",
     $        ' data without',/7x,'volume (V) and dV/dP data.')
              nerr = nerr + 1
            endif
          endif
          go to 350
        endif
      enddo
  350 continue
c
      if (ipch.gt.0 .and. ipcv.lt.0) then
        write (noutpt,1440)
        write (nttyo,1440)
 1440   format(/" * Error - (EQPT/rdwttl) Can't have pressure",
     $  ' corrections',
     $  /7x,'for enthalpy functions without having such corrections',
     $  /7x,'for log K functions. At a minimum, volume functions',
     $  /7x,'must be present.')
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write flags for enthalpy and volume functions grids.
c
      ux56 = 'Enthalpy functions flag'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1450) ux56(1:j2)
 1450 format(a)
c
      write (ndata1) ipch
      write (ndat1f,1460) ipch
 1460 format(i5)
c
      upchst = 'Not present'
      if (ipch .eq. 0) then
        upchst = 'Present, no dH/dP data'
      elseif (ipch .eq. 1) then
        upchst = 'Present, with dH/dP data'
      elseif (ipch .ge. 2) then
        upchst = 'Present, with dH/dP-dnH/dPn data'
        write (upchst,'(21x,i1,4x,i1)') ipch
      endif
      j2 = ilnobl(upchst)
c
      write (noutpt,1470) ipch,upchst(1:j2)
      write (nttyo,1470) ipch,upchst(1:j2)
 1470 format(/' Enthalpy functions flag= ',i2,' (',a,')')
c
      ux56 = 'Volume functions flag'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1450) ux56(1:j2)
c
      write (ndata1) ipcv
      write (ndat1f,1460) ipcv
c
      upcvst = 'Not present'
      if (ipcv .eq. 0) then
        upcvst = 'Present, no dV/dP data'
      elseif (ipcv .eq. 1) then
        upcvst = 'Present, with dV/dP data'
      elseif (ipcv .ge. 2) then
        upcvst = 'Present, with dV/dP-dnV/dPn data'
        write (upcvst,'(21x,i1,4x,i1)') ipcv
      endif
      j2 = ilnobl(upcvst)
c
      write (noutpt,1480) ipcv,upcvst(1:j2)
      write (nttyo,1480) ipcv,upcvst(1:j2)
 1480 format(' Volume functions flag  = ',i2,' (',a,')')
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (uakey(1:8) .eq. 'Pitzer  ') then
c
c       Write the flag for Pitzer data block organization.
c
        ux56 = 'Pitzer data block organization flag'
        write (ndata1) ux56
        j2 = ilnobl(ux56)
        write (ndat1f,1450) ux56(1:j2)
c
        write (ndata1) jpdblo
        write (ndat1f,1460) jpdblo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the temperature function used to represent Pitzer
c     interaction parameters. This is determined by the value
c     of jptffl, which was read previously by ggridp.f.
c
      if (uakey(1:8) .eq. 'Pitzer  ') then
        if (jpdblo .eq. -1) then
          if (jptffl .ne. -1) then
            write (noutpt,1510)
            write (nttyo,1510)
 1510       format(/' * Note - (EQPT/rdwttl) Resetting the jptffl',
     $      ' Pitzer data',/7x,'temperature function flag to -1',
     $      " (classic 25C-centric Taylor's series",/7x,'truncated',
     $      ' at second order), because the Pitzer data block',
     $      /7x,'organization flag (jpdblo) is set to -1',
     $      ' ("Classical").',/7x,'Other Pitzer data temperature',
     $      ' functions are not permitted.')
            jptffl = -1
          endif
        endif
c
        if (jptffl .eq. -1) then
          if (jpfcmx .ne. 3) then
            write (noutpt,1520)
            write (nttyo,1520)
 1520       format(/' * Error - (EQPT/rdwttl) The number of terms',
     $      ' (jpfcmx) used in',/7x,'the Pitzer data temperature',
     $      ' function is not 3, as is required',/7x,'in the case',
     $      " of the classical 25C-centric Taylor's series",
     $      /7x,'truncated at second order (jptffl = -1).')
            nerr = nerr + 1
          endif
        endif
c
        if (jptffl .eq. 1) then
          if (jpfcmx .ne. 8) then
            write (noutpt,1530)
            write (nttyo,1530)
 1530       format(/' * Error - (EQPT/rdwttl) The number of terms',
     $      ' (jpfcmx) used in',/7x,'the Pitzer data temperature',
     $      ' function is not 8, as is required',/7x,'in the case',
     $      ' of the case of the Greenberg-Moller combination',
     $      /7x,'temperature function (jptffl = 1).')
            nerr = nerr + 1
          endif
        endif
c
        if (jptffl .eq. 0) then
          if (jpfcmx.lt.1 .or. jpfcmx.gt.5) then
            write (noutpt,1540)
            write (nttyo,1540)
 1540       format(/' * Error - (EQPT/rdwttl) The number of terms',
     $      ' terms (jpfcmx)',/7x,'used in the LLNL maximal five-term',
     $      ' temperature equation (jptffl = 0)',/7x,'is not in the',
     $      ' allowed range of of 1-5. ')
            nerr = nerr + 1
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (uakey(1:8) .eq. 'Pitzer  ') then
c
c       Write the flag for the temperature function used to represent
c       Pitzer interaction coefficients.
c
        ux56 = 'Pitzer parameter temperature function'
        write (ndata1) ux56
        j2 = ilnobl(ux56)
        write (ndat1f,1450) ux56(1:j2)
c
        write (ndata1) jptffl
        write (ndat1f,1460) jptffl
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (ndat1f,1220) utermc(1:72)
 1220 format(a)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nerr .gt. 0) stop
c
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 write (noutpt,2000)
      write (nttyo,2000)
 2000 format(/' * Error - (EQPT/rdwttl) Unexpectedly encountered',
     $ /7x,'end-of-file while reading the title of the DATA0 file.')
      if (n .le. 1) then
        write (noutpt,2010)
        write (nttyo,2010)
 2010   format(/7x,'This occurred while trying to read the first line',
     $  /7x,'of the title.')
      else
        ulbufa = utitld(n - 1)
        call lejust(ulbufa)
        j2 = ilnobl(ulbufa)
        j2 = min(j2,70)
        write (noutpt,2020) ulbufa(1:j2)
        write (nttyo,2020) ulbufa(1:j2)
 2020   format(/7x,'This occurred while trying to read the line',
     $  ' following',/7x,'the line:',/7x,'"',a,'"')
      endif
      stop
c
  995 write (noutpt,2030)
      write (nttyo,2030)
 2030 format(/' * Error - (EQPT/rdwttl) Encountered a read format',
     $ /7x,'error while reading the DATA0 file.')
      if (n .le. 1) then
        write (noutpt,2010)
        write (nttyo,2010)
      else
        ulbufa = utitld(n - 1)
        call lejust(ulbufa)
        j2 = ilnobl(ulbufa)
        j2 = min(j2,70)
        write (noutpt,2020) ulbufa(1:j2)
      endif
      stop
c
  997 ulbufa = utitld(n)
      call lejust(ulbufa)
      j2 = ilnobl(ulbufa)
      j2 = min(j2,70)
      write (noutpt,2050) ulbufa(1:j2)
      write (nttyo,2050) ulbufa(1:j2)
 2050 format(/' * Error - (EQPT/rdwttl) Encountered a read format',
     $ /7x,'error while reading data embedded in the title of the',
     $ /7x,'DATA0 file. This occurred while attempting to process',
     $ ' the line:',/7x,'"',a,'"')
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
