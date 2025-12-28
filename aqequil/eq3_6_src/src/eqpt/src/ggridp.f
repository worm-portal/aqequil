      subroutine ggridp(ipch_asv,ipcv_asv,itgenf,jpdblo,jpfc_asv,
     $ jptffl,narx_asv,ndat0s,ndb_asv,noutpt,ntid_asv,ntpr_asv,
     $ nttyo,q500fl,uakey)
c
c     This suboutine gets the necessary dimensioning parameters for the
c     temperature grid on which thermodynamic data are represented.
c     These parameters are presently embedded in the data file title,
c     in order to maintain maximum forward and backward compatibility
c     of data files between Versions 7 and 8 of EQ3/6. There is
c     compatibility if the temperature grid matches the old "standard"
c     form and if there are no data grids for volume and enthalpy data.
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
c       uakey  = string identifying the type of aqueous species
c                activity coefficient model ("SEDH" or "Pitzer")
c
c
c     Principal output:
c
c       ipch_asv = the maximum order of pressure corrections to the
c                    enthalpy functions; ipch_asv = min(1,ipch)
c       ipcv_asv = the maximum order of pressure corrections to the
c       itgenf   = temperature grid enforcement parameter:
c                    -1 = do nothing
c                     0 = warn if any sparsely filled ranges in
c                         "log K" temperature grid data
c                     1 = declare error condition if any sparsely
c                           filled ranges in "log K" temperature
c                           grid data
c       jpdblo   = Pitzer data block organization flag; -1 =
c                    "classical", 0 = new
c                    volume functions; ipcv_asv = min(1,ipcv)
c       jpfc_asv = the number of terms in the temperature function
c                    used to represent Pitzer interaction parameters;
c                    must be at least 1
c       jptffl   = integer flag denoting the temperature function
c                    used to represent Pitzer interaction parameters;
c                    -1 = classical truncated Taylor's series,
c                    0 = LLNL5TERM, 1 = Greenberg and Moller (1989)
c                    eight-term equation
c       narx_asv = the maximum number of points in any range of the
c                    temperature grid
c       ntid_asv = the number of lines in the data file title
c       ntpr_asv = the number of temperature ranges on the temperature
c                    grid (ntpr_asv = ntprt)
c       ndb_asv  = the maximum number of distinct points on the
c                    temperature grid;
c                    ndb_asv = ntpr_asv*(narx_asv - 1) + 1
c       q500fl   = logical flag; if .true., treat instances of "500."
c                    in "log K" temperature grids as a "no data"
c                    condition"
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ndat0s,noutpt,nttyo
c
      integer ipch_asv,ipcv_asv,jpfc_asv,narx_asv,ndb_asv,ntid_asv,
     $ ntpr_asv
c
      integer itgenf,jpdblo,jptffl
c
      character(len=8) uakey
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j,j2,k,n,narxti,narx_x,ntpr,ntpri
c
      integer ilnobl
c
      logical q500fl,qnarxr,qrderr
c
      character(len=80), dimension(:), allocatable :: ux80ar
c
      character(len=80) ulbufa,ulbufb,uline
      character(len=16) ustr16,uterm
      character(len=8) ustr
c
c-----------------------------------------------------------------------
c
      data uterm  /'+---------------'/
c
c-----------------------------------------------------------------------
c
c     Determine the number of lines in the data file title.
c     Count that first terminator line as the last line of the title.
c
      n = 0
  100 read (ndat0s,1000,end=990,err=995) uline
 1000 format(a)
      n = n + 1
      if (uline(1:16) .ne. uterm(1:16)) go to 100
c
      ntid_asv = n
      rewind(ndat0s)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Allocate an array to hold the lines in the title.
c
      ALLOCATE(ux80ar(ntid_asv))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the title into the holding array.
c
      do n = 1,ntid_asv
        read (ndat0s,1000,end=990,err=995) ux80ar(n)
      enddo
c
      rewind(ndat0s)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set parameters for the temperature grid used to represent log K
c     values and standard state enthalpy and volume functions.
c
c     Values corresponding to the "classical" EQ3/6 temperature grid
c     will be used as defaults. These are intended to support the use
c     of version 7 format data files by the version 8 and higher
c     software.
c
c       ntprt       = the number of temperature ranges
c       narxt(ntpr) = the number of points or coefficients in the
c                       ntpr-th temperature range
c
c     The old standard grid corresponds to:
c
c       ntprt = 2
c       narxt(1) = 4  (0, 25, 60 100 C)
c       narxt(2) = 5  (100, 150, 200, 250, 300 C)
c       ipch = -1   (no enthalpy function grids are present)
c       ipcv = -1   (no volume function grids are present)
c
c     Note that the last point of the ntpr-th temperature range
c     is also the first point of the following range, if any.
c     This insures continuity.
c
c     Hence the matching default dimensioning is as follows. Note
c     that if ipch < 1, ipch_asv must be set to one, as an array
c     can't have a null or negative dimension. The same applies to
c     ipcv and ipcv_asv.
c
      ntpr_asv = 2
      narx_asv = 5
      ipch_asv = 1
      ipcv_asv = 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the number of temperature ranges (ntpr_asv).
c
      do n = 1,ntid_asv
        ulbufa = ' '
        ulbufb = ' '
        j = 0
c
c       Check for a keystring.
c
        i = index(ux80ar(n),'NO. OF TEMPERATURE RANGES')
        if (i .gt. 0) then
          j = i + 25
        else
          i = index(ux80ar(n),'NUMBER OF TEMPERATURE RANGES')
          if (i .gt. 0) then
            j = i + 28
          else
            i = index(ux80ar(n),'NTPRT')
            if (i .gt. 0) j = i + 5
          endif
        endif
c
        if (j .gt. 0) then
c
c         Check for an equal sign following the keystring.
c
          ulbufa = ux80ar(n)(j:80)
          call lejust(ulbufa)
          i = index(ulbufa,'=')
          j = 2
          if (i .ne. 1) j = 0
        endif
c
        if (j .gt. 0) then
c
c         Extract a number matching the keystring.
c
          ulbufb = ulbufa(2:80)
          call lejust(ulbufb)
          ustr = ulbufb(1:8)
          call chrint(ntpr_asv,nttyo,qrderr,ustr)
          if (qrderr) go to 997
          go to 200
        endif
      enddo
c
  200 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the maximum number of points in any temperature range
c     (narx_asv). Note that this is done by finding the number of
c     points in each range and determining the maximum. However, an
c     array to hold the actual values has not yet been allocated.
c     The actual values will be redetermined and stored by rdwttl.f.
c
      qnarxr = .false.
      narx_x = 0
      narxti = 0
      do ntpr = 1,ntpr_asv
        do n = 1,ntid_asv
          ulbufa = ' '
          ulbufb = ' '
          j = 0
c
c         Check for a keystring.
c
          i = index(ux80ar(n),'NO. OF POINTS IN RANGE')
          if (i .gt. 0) then
            j = i + 22
          else
             i = index(ux80ar(n),'NUMBER OF POINTS IN RANGE')
            if (i .gt. 0) then
              j = i + 25
            else
              i = index(ux80ar(n),'NARXT')
              if (i .gt. 0) j = i + 5
            endif
          endif
c
          if (j .gt. 0) then
c
c           Extract a number acting as a subscript to the keystring.
c
            ulbufb = ux80ar(n)(j:80)
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
            qnarxr = .true.
            narx_x = max(narx_x,narxti)
            go to 220
          endif
        enddo
  220   continue
c
      enddo
c
c     If a value was read for any temperature range, use it. But
c     do not use a number less than 1, as an array can not have a
c     null or negative dimension.
c
      if (qnarxr) narx_asv = narx_x
      narx_asv = max(1,narx_asv)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the maximum possible number of points on the temperature
c     grid for the current grid parameters.
c
      ndb_asv = ntpr_asv*(narx_asv - 1) + 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set flag for dealing with sparsely filled ranges in the
c     "log K" temperature grids.
c       itgenf:
c         -1 = do nothing
c          0 = generate warnings for sparsely filled ranges
c          1 = generate errors for sparsely filled ranges
c
      itgenf = -1
      do n = 1,ntid_asv
        ulbufa = ' '
        j = 0
c
c       Check for a keystring.
c
        i = index(ux80ar(n),'SPARSE GRID RANGE CONDITION')
        if (i .gt. 0) j = i + 27
c
        if (j .gt. 0) then
c
c         Check for an equal sign following the keystring.
c
          ulbufa = ux80ar(n)(j:80)
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
          ustr16 = ulbufa(2:17)
          call lejust(ustr16)
          k = index(ustr16,'IGNORE')
          if (k .gt. 0) then
            itgenf = -1
          else
            k = index(ustr16,'WARN')
            if (k .gt. 0) then
              itgenf = 0
            else
              k = index(ustr16,'ERROR')
              if (k .gt. 0) then
                itgenf = 1
              endif
            endif
          endif
          go to 300
        endif
      enddo
  300 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set flag for dealing with values of "500." on the "log K"
c     temperature grids.
c       q500fl:
c         .false. = do nothing
c          .true. = treat values of 500. as "no data"
c
      q500fl = .true.
      do n = 1,ntid_asv
        ulbufa = ' '
        j = 0
c
c       Check for a keystring.
c
        i = index(ux80ar(n),'INTERPRET 500 AS NO DATA')
        if (i .gt. 0) j = i + 24
c
        if (j .gt. 0) then
c
c         Check for an equal sign following the keystring.
c
          ulbufa = ux80ar(n)(j:80)
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
          ustr16 = ulbufa(2:17)
          call lejust(ustr16)
          k = index(ustr16,'NO')
          if (k .le. 0) k = index(ustr16,'FALSE')
          if (k .gt. 0) then
            q500fl = .false.
          else
            k = index(ustr16,'YES')
            if (k .le. 0) k = index(ustr16,'TRUE')
            if (k .gt. 0) then
              q500fl = .true.
            endif
          endif
          go to 310
        endif
      enddo
  310 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the dimension corresponding to the order of the pressure
c     dependence of the enthalpy functions. An order of -1 means that
c     No enthalpy functions data are present on the data file. An order
c     of 0 means that enthalpy functions data are present, but no data
c     for the pressure dependence of the same are present. However,
c     the dimensioning variable ipch_asv must be set to a value of
c     at least one, as an array can not have a null or negative
c     dimension. Thus, it is necessary here to check only for data
c     corresponding to order 2 or higher.
c
      do n = 1,ntid_asv
        ulbufa = ' '
        j = 0
c
c       Check for a keystring.
c
        i = index(ux80ar(n),'d2H/dP2')
        if (i .gt. 0) j = i + 7
c
        if (j .gt. 0) then
c
c         Check for an equal sign following the keystring.
c
          ulbufa = ux80ar(n)(j:80)
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
          if (k .gt. 0) ipch_asv = 2
          go to 320
        endif
      enddo
  320 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the dimension corresponding to the order of the pressure
c     dependence of the volume functions. This process and the dimension
c     itself are analogous to those described above for the pressure
c     dependence of the enthalpy functions.
c
      do n = 1,ntid_asv
        ulbufa = ' '
        j = 0
c
c       Check for a keystring.
c
        i = index(ux80ar(n),'d2V/dP2')
        if (i .gt. 0) j = i + 7
c
        if (j .gt. 0) then
c
c         Check for an equal sign following the keystring.
c
          ulbufa = ux80ar(n)(j:80)
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
          if (k .gt. 0) ipcv_asv = 2
          go to 350
        endif
      enddo
  350 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set flag for Pitzer data block organization.
c       jpdblo:
c         = -1 Classical: two superblocks, one for pure aqueous
c              electrolytes, one for mixtures of two such electrolytes
c              having a common ion. A given theta parameter can appear
c              multiple times. Data for parameters involving electrically
c              neutral species are folded into this structure. Checks are
c              in place to sort things out.
c         =  0 New: multiple superblocks, each corresponding to a
c              charge-based type of species pair or triplet. This
c              structure is more rational. Theta parameters for example
c              have their own superblock.
c
      jpdblo = -1
      if (uakey(1:8) .eq. 'Pitzer  ') then
        do n = 1,ntid_asv
          ulbufa = ' '
          j = 0
c
c         Check for a keystring.
c
          i = index(ux80ar(n),'PITZER DATA BLOCK ORG.')
          if (i .gt. 0) j = i + 22
c
          if (j .gt. 0) then
c
c           Check for an equal sign following the keystring.
c
            ulbufa = ux80ar(n)(j:80)
            call lejust(ulbufa)
            i = index(ulbufa,'=')
            j = 2
            if (i .ne. 1) j = 0
          endif
c
          if (j .gt. 0) then
c
c           Extract a string input matching the keystring.
c
            ustr16 = ulbufa(2:17)
            call lejust(ustr16)
            k = index(ustr16,'NEW')
            if (k .le. 0) k = index(ustr16,'BYSPECIES')
            if (k .gt. 0) jpdblo = 0
            go to 400
          endif
        enddo
  400   continue
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the identification of the temperature function used to
c     represent Pitzer interaction parameters. Note that the number
c     of terms or coefficients (jpfcmx) is a dimensioning parameter
c     that may also be set by an option embedded in the data file
c     title.
c
c     Temperature function flag.
c       jptffl:
c         = -1 Classic: a Taylor's series centered at 25C and
c              truncated at second order.
c
c                a1 = value at 25C
c                a2 = first derivative at 25C
c                a3 = second derivative at 25C
c
c              This is currently required if the Pitzer data block
c              organization flag is set to "Classical" (jpdblo = -1).
c              The number of terms used (jpfcmx) is fixed at 3.
c
c         =  0 The function:
c
c                x(T) = a1 + a2/(T - Tr) + a3*ln(T/Tr)
c                            a4(T - Tr) + a5(T**2 - Tr**2)
c
c              This is also 25C centric (Tr = 25C).
c              The number of terms used (jpfcmx) may vary from
c              1 to 5. This always refers to the first jpfcmx terms.
c
c         =  1 The Greenberg-Moller (1988) combination temperature
c              function. This is not 25C centric. The number of terms
c              used (mx) is fixed at 8.
c
      jptffl = -1
      if (uakey(1:8) .eq. 'Pitzer  ') then
        do n = 1,ntid_asv
          ulbufa = ' '
          j = 0
c
c         Check for a keystring.
c
          i = index(ux80ar(n),'PITZER TEMP FUNCTION')
          if (i .gt. 0) j = i + 20
c
          if (j .gt. 0) then
c
c           Check for an equal sign following the keystring.
c
            ulbufa = ux80ar(n)(j:80)
            call lejust(ulbufa)
            i = index(ulbufa,'=')
            j = 2
            if (i .ne. 1) j = 0
          endif
c
          if (j .gt. 0) then
c
c           Extract a string input matching the keystring.
c
            ustr16 = ulbufa(2:17)
            call lejust(ustr16)
            k = index(ustr16,'NEW')
            if (k .le. 0) k = index(ustr16,'Livermore')
            if (k .le. 0) k = index(ustr16,'LIVERMORE')
            if (k .le. 0) k = index(ustr16,'LLNL')
            if (k .le. 0) k = index(ustr16,'L5TERM')
            if (k .le. 0) k = index(ustr16,'5TERM')
            if (k .gt. 0) then
              jptffl = 0
              go to 410
            endif
            k = index(ustr16,'TEQUIL')
            if (k .le. 0) k = index(ustr16,'T8TERM')
            if (k .le. 0) k = index(ustr16,'GM8TERM')
            if (k .le. 0) k = index(ustr16,'8TERM')
            if (k .gt. 0) then
              jptffl = 1
              go to 410
            endif
          endif
        enddo
  410   continue
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the number of terms or coefficients (jpfc_asv) in the
c     temperature function used to represent Pitzer parameters.
c
      jpfc_asv = 1
      if (uakey(1:8) .eq. 'Pitzer  ') then
c
c       Three is the default value for a Pitzer model.
c
        jpfc_asv = 3
        do n = 1,ntid_asv
          ulbufa = ' '
          j = 0
c
c         Check for a keystring.
c
          i = index(ux80ar(n),'NO. OF PITZER TEMP FUNC TERMS')
          if (i .gt. 0) j = i + 29
c
          if (j .gt. 0) then
c
c           Check for an equal sign following the keystring.
c
            ulbufa = ux80ar(n)(j:80)
            call lejust(ulbufa)
            i = index(ulbufa,'=')
            j = 2
            if (i .ne. 1) j = 0
          endif
c
          if (j .gt. 0) then
c
c           Extract a string input matching the keystring.
c
            ustr = ulbufa(2:9)
            call lejust(ustr)
            read (ustr,'(i8)') jpfc_asv
            go to 420
          endif
        enddo
  420   continue
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Deallocate the array holding the data file title.
c
      DEALLOCATE(ux80ar)
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write messages for errors associated with reading the DATA0
c     file.
c
  990 write (noutpt,2000)
      write (nttyo,2000)
 2000 format(/' * Error - (EQPT/ggridp) Unexpectedly encountered',
     $ /7x,'end-of-file while reading the title of the DATA0 file.')
      if (n .le. 1) then
        write (noutpt,2010)
        write (nttyo,2010)
 2010   format(/7x,'This occurred while trying to read the first line',
     $  /7x,'of the title.')
      else
        ulbufa = ux80ar(n - 1)
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
 2030 format(/' * Error - (EQPT/ggridp) Encountered a read format',
     $ /7x,'error while reading the DATA0 file.')
      if (n .le. 1) then
        write (noutpt,2010)
        write (nttyo,2010)
      else
        ulbufa = ux80ar(n - 1)
        call lejust(ulbufa)
        j2 = ilnobl(ulbufa)
        j2 = min(j2,70)
        write (noutpt,2020) ulbufa(1:j2)
      endif
      stop
c
  997 ulbufa = ux80ar(n)
      call lejust(ulbufa)
      j2 = ilnobl(ulbufa)
      j2 = min(j2,70)
      write (noutpt,2050) ulbufa(1:j2)
      write (nttyo,2050) ulbufa(1:j2)
 2050 format(/' * Error - (EQPT/ggridp) Encountered a read format',
     $ /7x,'error while reading data embedded in the title of the',
     $ /7x,'DATA0 file. This occurred while attempting to process',
     $ ' the line:',/7x,'"',a,'"')
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
