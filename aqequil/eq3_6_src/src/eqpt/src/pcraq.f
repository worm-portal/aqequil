      subroutine pcraq(aamatr,apr,atwt,avgrid,cdrs,cdrsi,cess,cessi,
     $ cof,dhfe,dhfs,dvfe,dvfs,eps100,gmmatr,ipch,ipchmx,ipcv,ipcvmx,
     $ ipivot,itgenf,mtotr,nacdpr,narxmx,narxt,nat,natmax,nbt,nbtmx1,
     $ nbtmx2,nch,nco,nct,nctmax,ndata1,ndat0s,ndat1f,ndbmax,ndbptg,
     $ ndbptl,nentei,nentri,nerr,nmodwr,noutpt,nsb,nslist,ntprmx,ntprt,
     $ nttyo,nwarn,qelect,q500fl,tempc,tempcs,tmpcmx,uaqsp,udbfmt,
     $ udbval,udrsi,uelem,uessi,uspec,xdbval,xhfe,xhfs,xlke,xlks,xvfe,
     $ xvfs,xvec,yvec,zaqsp,zchar)
c
c     This subroutine reads data on aqueous species from the stripped
c     DATA0 file, processes this data, and writes the results on the
c     DATA1 and DATA1F files. The counter nerr is incremented by
c     one for each error encountered. The counter nwarn is incremented
c     similarly for each warning.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       narxt  = array of numbers of coefficients in the temperature
c                  ranges
c       nbt    = the number of basis species
c       nct    = the number of chemical elements
c       ndata1 = unit number of the DATA1 file
c       ndat0s = unit number of the stripped DATA0 file
c       ndat1f = unit number of the DATA1F file
c       ndbptg = the number of distinct points on the temperature grid
c       ndbptl = the maximum number of points on the temperature grid
c                  per line
c       ntprt  = the number of temperature ranges on the standard
c                  temperature grid
c       qelect = flag denoting the use of "e-" instead of "O2(g)"
c                  in writing chemical reactions
c       tempc  = array of temperatures (on the temperature grid)
c       tempcs = array of scaled temperatures (on the temperature grid)
c       tmpcmx = the max norm of the tempeatures on the temperature grid
c       udbfmt = the format for reading a line of data on the
c                  temperature grid
c
c
c     Principal output:
c
c       atwt   = array of atomic weights
c       cdrs   = array of reaction coefficients
c       cess   = array of elemental composition coefficients
c       nat    = the number of aqueous species
c       nerr   = cumulative error counter
c       nwarn  = cumulative warning counter
c       uaqsp  = array of aqueous species names
c       uelem  = array of names of chemical elements
c       uspec  = array of names of species
c       xlke   = array of log K values for the "Eh" reaction (on the
c                  temperature grid)
c       xlks   = array of log K values for reactions
c       zaqsp  = array of charge numbers of aqueous species
c       zchar  = array of electrical charge numbers of species
c
c     Workspace:
c
c       aamatr = matrix used to calculate the polynomial coefficients
c       apr    = array of polynomial coefficients (for all ranges)
c       avgrid = array containing the data on the temperature grid
c       cof    = array of fitted polynomial coefficients (for a
c                  single temperature range)
c       gmmatr = a copy of the amatr matrix
c       ipivot = the pivot vector, used in solving matrix equations
c       nacdpr = array containging the number of actual data points
c                  by range on the "log K" temperature grid; excludes
c       udbval = string array for reading data on the temperature grid
c       xdbval = holding space array for data on the temperature grid
c       xvec   = array of scaled temperatures corresponding to the
c                  data in the yvec array
c       yvec   = array of data to be fitted (for a single temperature
c                  range)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipchmx,ipcvmx,narxmx,natmax,nbtmx1,nbtmx2,nctmax,
     $ ndbmax,ntprmx
c
      integer ndata1,ndat0s,ndat1f,noutpt,nslist,nttyo
c
      integer ipivot(narxmx),nacdpr(ntprmx),narxt(ntprmx),
     $ nentei(nctmax),nentri(nbtmx1)
c
      integer ipch,ipcv,itgenf,nat,nbt,nch,nco,nct,ndbptg,ndbptl,nerr,
     $ nmodwr,nsb,ntprt,nwarn
c
      logical qelect,q500fl
c
      character(len=24) uaqsp(natmax),udrsi(nbtmx1),uspec(nbtmx1)
      character(len=16) udbval(ndbmax)
      character(len=8) uelem(nctmax),uessi(nctmax)
c
      character(len=16) udbfmt
c
      real(8) atwt(nctmax),cdrs(nbtmx2,nbtmx1),cdrsi(nbtmx1),
     $ cess(nctmax,nbtmx1),cessi(nctmax),dhfe(narxmx,ntprmx,ipchmx),
     $ dhfs(narxmx,ntprmx,ipchmx,nbtmx1),dvfe(narxmx,ntprmx,ipcvmx),
     $ dvfs(narxmx,ntprmx,ipcvmx,nbtmx1),xhfe(narxmx,ntprmx),
     $ xhfs(narxmx,ntprmx,nbtmx1),xlke(narxmx,ntprmx),
     $ xlks(narxmx,ntprmx,nbtmx1),xvfe(narxmx,ntprmx),
     $ xvfs(narxmx,ntprmx,nbtmx1),zchar(nbtmx1),zaqsp(natmax)
c
      real(8) apr(narxmx,ntprmx),avgrid(narxmx,ntprmx)
      real(8) tempc(narxmx,ntprmx),tempcs(narxmx,ntprmx),tmpcmx(ntprmx)
      real(8) aamatr(narxmx,narxmx),gmmatr(narxmx,narxmx)
      real(8) cof(narxmx),xvec(narxmx),yvec(narxmx)
      real(8) xdbval(ndbmax),mtotr(nctmax)
c
      real(8) eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ii,ipc,j,j2,j3,j4,j5,k,n,nbt1,nbw,nc,ncount,ncts,ndrsts,
     $ nmodx,nn,nnx,ns,nse,nsm1,nspnx,nt,ntpr,nxm
c
      integer ilnobl
c
      logical qblkes,qblkrs,qend,qerr,qnofes,qnofrs,qzeres,qzerrs,
     $ q500nd
c
      character(len=24) uspn(2)
c
      character(len=80) ulbufa,ulbufb,uline,ux80
      character(len=72) uterm,utermc
      character(len=56) ustrgr
      character(len=24) uaq24,uend24,ublk24,unone,usblkf,ux24
      character(len=16) ux16
      character(len=8) uaqu,uaux,ux8,ux8a,ux8b,ux8c,uendit,usol
c
      real(8) mwtss,zx
c
c-----------------------------------------------------------------------
c
      data ublk24 /'                        '/
      data uaq24 /'aqueous                 '/
      data unone /'none                    '/
      data uendit /'endit.  '/
      data uaux /'auxiliar'/
      data uaqu /'aqueous '/
      data usol /'solids  '/
c
c-----------------------------------------------------------------------
c
      uterm(1:48) = '+-----------------------------------------------'
      uterm(49:72) = '------------------------'
      utermc = uterm
      utermc(1:1) = '*'
c
c     The variable nn must be set to zero here to provide a proper
c     diagnostic if an early end-of-file or read format error occurs.
c
      nn = 0
c
      q500nd = q500fl
c
c     Set the expected index of H2O.
c
      nbw = 1
c
c     Skip 2 lines to the 'basis species' delimiter.
c
      read (ndat0s,1000,end=990,err=995) uline
      read (ndat0s,1000,end=990,err=995) uline
 1000 format(a)
c
c     Write label.
c
      usblkf = uaq24
      j3 = ilnobl(uaq24)
      write (ndata1) uaq24
      write (ndat1f,1010) uaq24(1:j3)
 1010 format(a)
      write (ndat1f,1010) utermc(1:72)
      write (noutpt,1020) uaq24(1:j3)
      write (nttyo,1020) uaq24(1:j3)
      write (nslist,1020) uaq24(1:j3)
 1020 format(//1x,a,/)
c
c     Initialize some variables.
c
      qelect = .false.
      uend24(1:8) = uendit(1:8)
      uend24(9:24) = ublk24(9:24)
      uspn(1) = unone
      uspn(2) = ublk24
      nspnx = 1
      nnx = -1
      nn = 1
      nmodx = 0
      nbt1 = nbt + 1
      nat = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Note: the main loop returns here.
c
  100 ns = min(nn,nbt1)
      ncts = 0
      ndrsts = 0
c
      do nc = 1,nct
        cess(nc,ns) = 0.
        cessi(nc) = 0.
      enddo
c
      do nse = 1,nbt1
        cdrs(nse,ns) = 0.
        cdrsi(nse) = 0.
      enddo
      cdrs(nbt + 2,ns) = 0.
c
      call initcv(uessi,nct,' ')
      call initcv(udrsi,nbt1,' ')
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the first line of the block.
c
  110 continue
      read (ndat0s,1000,end=990,err=995) uline
      ux24 = uline(1:24)
      ux8 = ux24(1:8)
c
      if (ux8(1:8).eq.uaux(1:8) .or. ux8(1:8).eq.uaqu(1:8)) then
c
c       Skip the terminator line, read the next line.
c
        read (ndat0s,1000,end=990,err=995) uline
        go to 110
c
      elseif (ux8(1:8) .eq. usol(1:8)) then
c
c       Skip the terminator line.
c
        read (ndat0s,1000,end=990,err=995) uline
c
c       Finish and exit here.
c
        if (uspn(1)(1:24) .eq. unone(1:24)) then
          nnx = 0
          nn = 0
          write (nttyo,1140) nnx,unone(1:4)
          write (nslist,1140) nnx,unone(1:4)
 1140     format(1x,i5,2x,a)
          write (noutpt,1150) nn,unone(1:4)
 1150     format(1x,i5,1x,a)
        else
          nxm = nspnx - 1
          if (nxm .ge. 1) then
            nnx = nnx + 2
            if (nxm .eq. 1) then
              j3 = ilnobl(uspn(1))
              write (nttyo,1140) nnx,uspn(1)(1:j3)
              write (nslist,1140) nnx,uspn(1)(1:j3)
            else
              j3 = ilnobl(uspn(2))
              write (nttyo,1160) nnx,uspn(1),uspn(2)(1:j3)
              write (nslist,1160) nnx,uspn(1),uspn(2)(1:j3)
 1160         format(1x,i5,2x,a24,2x,a)
            endif
          elseif (nmodx .ne. 1) then
            j3 = ilnobl(uspn(2))
            write (nttyo,1160) nnx,uspn(1),uspn(2)(1:j3)
          endif
        endif
        write (ndata1) uend24,ublk24,ublk24
        j3 = ilnobl(uend24)
        write (ndat1f,1010) uend24(1:j3)
        write (ndat1f,1010) utermc(1:72)
        if (nerr .gt. 0) stop
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Have a species block.
c
      nat = nat + 1
c
      j2 = ilnobl(ux24)
      if (j2 .le. 0) then
        write (ux8a,'(i5)') nat
        call lejust(ux8a)
        j3 = ilnobl(ux8a)
        j4 = ilnobl(usblkf)
        write (noutpt,1240) ux8a(1:j3),usblkf(1:j4)
        write (nttyo,1240) ux8a(1:j3),usblkf(1:j4)
 1240   format(/' * Error - (EQPT/pcraq) Have encountered a blank',
     $  ' species name',/7x,'for species block ',a,' of the ',a,
     $  ' superblock.')
        if (nat .gt. 1) then
          ux24 = uaqsp(nat - 1)
          j5 = ilnobl(ux24)
          if (j5 .gt. 0) then
            write (noutpt,1250) ux24(1:j5)
            write (nttyo,1250) ux24(1:j5)
 1250       format(7x,'This block follows the one for ',a,'.')
          endif
        endif
        ux24 = '<blank>'
        nerr = nerr + 1
      endif
c
      uaqsp(nat) = ux24
      uspec(ns) = ux24
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Save the species name for screen and SLIST output.
c
      uspn(nspnx) = uspec(ns)
c
c     Write the species name on the OUTPUT file.
c
      j2 = ilnobl(uspec(ns))
      write (noutpt,1150) nn,uspec(ns)(1:j2)
c
c     Write species names two per line on the SLIST file.
c
      nspnx = nspnx + 1
      if (nspnx .gt. 2) then
        nspnx = 1
        nnx = nnx + 2
        nmodx = mod(nnx,nmodwr)
        if (nmodwr .eq. 1) nmodx = 1
        j3 = ilnobl(uspn(2))
        if (nmodx .eq. 1) write (nttyo,1160) nnx,uspn(1),uspn(2)(1:j3)
        write (nslist,1160) nnx,uspn(1),uspn(2)(1:j3)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Skip the obsolete 'sp.type =' and 'revised =' lines.
c     Read the electrical charge.
c
  130 read (ndat0s,1000,end=990,err=995) uline
      ii = index(uline,'charge')
      if (ii .eq. 0) go to 130
c
      ux80 = uline(ii + 6:80)
      ii = index(ux80,'=')
      ux16 = ux80(ii + 1:80)
      call lejust(ux16)
      read (ux16,'(f5.1)',err=995) zx
c1280 format(14x,f5.1)
      zaqsp(nat) = zx
      zchar(ns) = zx
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the number of chemical elements composing the species.
c
      read (ndat0s,1000,end=990,err=995) uline
      ux80 = uline
      call lejust(ux80)
      ii = index(ux80,'element(s)')
      if (ii .gt. 1) then
        ux16 = ux80(1:ii - 1)
        read (ux16,'(i5)',err=995) ncts
c1330   format(4x,i2)
      else
        ii = index(ux80,' ')
        if (ii .gt. 1) then
          ux16 = ux80(1:ii - 1)
          read (ux16,'(i5)',err=995) ncts
        else
          ncts = 0
        endif
      endif
c
      if (ncts .gt. nct) then
        j2 = ilnobl(uspec(ns))
        write (ux8a,'(i5)') ncts
        call lejust(ux8a)
        j3 = ilnobl(ux8a)
        write (ux8c,'(i5)') nct
        call lejust(ux8c)
        j5 = ilnobl(ux8c)
        write (noutpt,1350) uspec(ns)(1:j2),ux8a(1:j3),ux8c(1:j5)
        write (nttyo,1350) uspec(ns)(1:j2),ux8a(1:j3),ux8c(1:j5)
 1350   format(/' * Error - (EQPT/pcraq) Species "',a,'" is',
     $  ' composed of ',a,/7x,'chemical elements, but there',
     $  ' are only ',a,' elements on the data file.')
        nerr = nerr + 1
      endif
c
      if (ncts .le. 0) then
        if (uspec(ns)(1:3).ne.'e- ' .and. uspec(ns)(1:3).ne.'E- ') then
          j2 = ilnobl(uspec(ns))
          write (noutpt,1360) uspec(ns)(1:j2)
          write (nttyo,1360) uspec(ns)(1:j2)
 1360     format(/' * Error - (EQPT/pcraq) Species "',a,'" is not',
     $    ' composed of',/7x,'any chemical elements.')
          nerr = nerr + 1
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the names of the elements and the corresponding
c     compositional coefficients.
c
      n = 0
      do i = 1,ncts/3
        read (ndat0s,1000,end=990,err=995) uline
        ulbufa = uline
        do k = 1,3
          call lejust(ulbufa)
          ii = index(ulbufa,' ')
          if (ii .gt. 1) then
            ux16 = ulbufa(1:ii - 1)
            read (ux16,'(f8.4)',err=995) cessi(n + k)
          else
            cessi(n + k) = 0.
          endif
          ulbufb = ulbufa(ii + 1:80)
          call lejust(ulbufb)
          uessi(n + k) = ulbufb(1:8)
          if (k .lt. 3) ulbufa = ulbufb(9:80)
        enddo
c       read (uline,1370,err=995) (cessi(n + k),uessi(n + k), k = 1,3)
c1370   format((4x,3(f8.4,1x,a8,5x)))
        n = n + 3
      enddo
c
      j = mod(ncts,3)
      if (j .gt. 0) then
        read (ndat0s,1000,end=990,err=995) uline
        ulbufa = uline
        do k = 1,j
          call lejust(ulbufa)
          ii = index(ulbufa,' ')
          if (ii .gt. 1) then
            ux16 = ulbufa(1:ii - 1)
            read (ux16,'(f8.4)',err=995) cessi(n + k)
          else
            cessi(n + k) = 0.
          endif
          ulbufb = ulbufa(ii + 1:80)
          call lejust(ulbufb)
          uessi(n + k) = ulbufb(1:8)
          if (k .lt. j) ulbufa = ulbufb(9:80)
        enddo
c       read (uline,1370,err=995) (cessi(n + k),uessi(n + k), k = 1,j)
        n = n + j
      endif
c
c     Check for blank element names, duplicate element names, and
c     zero-valued stoichiometric coefficients.
c
      call elesck(cessi,nbtmx1,nctmax,ncts,nentei,nerr,noutpt,
     $ ns,nttyo,qblkes,qzeres,uessi,usblkf,uspec)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ns .eq. nbw) then
c
c       Make sure that H2O is the first aqueous species.
c
        if (uspec(nbw)(1:4) .ne. 'H2O ') then
          j2 = ilnobl(uspec(nbw))
          write (noutpt,1420) uspec(nbw)(1:j2)
          write (nttyo,1420) uspec(nbw)(1:j2)
 1420     format(/' * Error - (EQPT/pcraq) The first strict basis',
     $    ' species on the data file',/7x,'is ',a,'. The first',
     $    ' strict basis species must be H2O.',/7x,'This is an',
     $    ' idiosyncrasy of EQ3/6.')
          nerr = nerr + 1
        endif
      endif
c
      if (ns .le. nsb) then
c
c       A strict basis species may be composed of only one chemical
c       element other than O or H.
c
        ncount = 0
        do i = 1,ncts
          if (cessi(i) .gt. 0.) then
            if (uessi(i)(1:2).ne.'O ' .and. uessi(i)(1:2).ne.'H ' .and.
     $      uessi(i)(1:7).ne.'<blank>') ncount = ncount + 1
          endif
        enddo
c
        if (ncount .gt. 1) then
          j2 = ilnobl(uspec(ns))
          write (noutpt,1430) uspec(ns)(1:j2)
          write (nttyo,1430) uspec(ns)(1:j2)
 1430     format(/' * Error - (EQPT/pcraq) The strict basis species ',
     $    a,' is composed',/7x,'of more than one chemical element',
     $    ' other than O and H.')
          nerr = nerr + 1
        endif
      endif
c
      if (ns .le. nct) then
c
c       Make sure that there is a 1:1 mapping to the chemical elements
c       for the first nct members of the strict basis set.
c
c       H2O (solvent water) is a special case. H2O always maps to O.
c       Since the first chemical element is required to be O, the first
c       aqueous species is required to  be H2O. Here it is sufficient
c       to make sure that that indeed is the case.
c
        if (ns .eq. 1) then
          if (uspec(ns)(1:4) .ne. 'H2O ') then
            j2 = ilnobl(uspec(ns))
            write (noutpt,1440) uspec(ns)(1:j2)
            write (nttyo,1440) uspec(ns)(1:j2)
 1440       format(/' * Error - (EQPT/pcraq) The first aqueous species',
     $      ' is ',a,'.',/7x,'The first aqueous species must be H2O.')
            nerr = nerr + 1
          endif
        endif
c
c       The hydrogen ion is also a special case. H+ always maps to H.
c
        if (uspec(ns)(1:3) .eq. 'H+ ') then
          nc = ns
          if (uelem(nc)(1:2) .ne. 'H ') then
            j2 = ilnobl(uspec(ns))
            write (noutpt,1450)
            write (nttyo,1450)
 1450       format(/' * Error - (EQPT/pcraq) The species H+ is',
     $      ' required to match',/7x,'the chemical element H. If the',
     $      ' species is the n-th member of',/7x,'the strict',
     $      ' basis set, the corresponding chemical element must be',
     $      /7x,'the n-th chemical element (there is one exception',
     $      ' to this rule,',/7x,"the redox species, which doesn't",
     $      ' correspond to any chemical element).')
            nerr = nerr + 1
          endif
        endif
c
c       Any strict basis species, other than the redox species,
c       must be composed of the chemical element to which it
c       formally corresponds.
c
        nc = ns
        do i = 1,ncts
          if (uessi(i)(1:8) .eq. uelem(nc)(1:8)) then
            if (cessi(i) .gt. 0.) go to 120
          endif
        enddo
c
        j2 = ilnobl(uspec(ns))
        j3 = ilnobl(uelem(nc))
        write (noutpt,1460) uspec(ns)(1:j2),uelem(nc)(1:j3)
        write (nttyo,1460) uspec(ns)(1:j2),uelem(nc)(1:j3)
 1460   format(/' * Error - (EQPT/pcraq) The strict basis',
     $  ' species ',a,' is not',/7x,'composed of the chemical',
     $  ' element ',a,' to which it must correspond.',
     $  /7x,'If the species is the n-th member of the strict',
     $  ' basis set,',/7x,'the, corresponding chemical',
     $  ' element must be the n-th chemical element.',
     $  /7x,'There is one exception to this rule, the',
     $  ' redox species. It must be',/7x,'the last strict',
     $  " basis species, if present. It doesn't correspond",
     $  /7x,'to any chemical element.')
        nerr = nerr + 1
  120   continue
      endif
c
c     Check out the redox species.
c
      if (ns .eq. nsb) then
        if (uspec(ns)(1:6).ne.'O2(g) ' .and.
     $    uspec(ns)(1:3).ne.'e- ' .and.
     $    uspec(ns)(1:3).ne.'E- ') then
          j2 = ilnobl(uspec(ns))
          write (noutpt,1470) uspec(ns)(1:j2)
          write (nttyo,1470) uspec(ns)(1:j2)
 1470     format(/' * Error - (EQPT/pcraq) The redox species is ',a,
     $    '. This must',/7x,'be either O2(g) or e-.')
          nerr = nerr + 1
        endif
      endif
c
      if (ns .eq. nsb) then
        if (nerr .gt. 0) stop
      endif
c
      if (ns .eq. nbt) then
        if (nerr .gt. 0) stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ns .gt. nsb) then
c
c       Read number of species in the associated reaction.
c
        read (ndat0s,1000,end=990,err=995) uline
        ux80 = uline
        call lejust(ux80)
        ii = index(ux80,'species')
        if (ii .gt. 1) then
          ux16 = ux80(1:ii - 1)
          read (ux16,'(i5)',err=995) ndrsts
c1500   format(4x,i2)
        else
          ii = index(ux80,' ')
          if (ii .gt. 1) then
            ux16 = ux80(1:ii - 1)
            read (ux16,'(i5)',err=995) ndrsts
          else
            ndrsts = 0
          endif
        endif
c
        if (ndrsts .gt. nbt1) then
          j2 = ilnobl(uspec(ns))
          write (ux8a,'(i5)') ndrsts
          call lejust(ux8a)
          j3 = ilnobl(ux8a)
          write (ux8b,'(i5)') nbt
          call lejust(ux8b)
          j4 = ilnobl(ux8b)
          write (ux8c,'(i5)') nbt1
          call lejust(ux8c)
          j5 = ilnobl(ux8c)
          write (noutpt,1520) uspec(ns)(1:j2),ux8a(1:j3),ux8b(1:j4),
     $    ux8c(1:j5)
          write (nttyo,1520) uspec(ns)(1:j2),ux8a(1:j3),ux8b(1:j4),
     $    ux8c(1:j5)
 1520     format(/' * Error - (EQPT/pcraq) The reaction for the',
     $    /7x,'destruction of species "',a,'" includes ',a,
     $    /7x,'species, but there are only ',a,' basis species on the',
     $    /7x,'data file, so only ',a,' species may appear in the ',
     $    'reaction.')
          nerr = nerr + 1
        endif
c
        if (ndrsts .lt. 2) then
          j2 = ilnobl(uspec(ns))
          j4 = ilnobl(usblkf)
          write (noutpt,1530) uspec(ns)(1:j2),usblkf(1:j4)
          write (nttyo,1530) uspec(ns)(1:j2),usblkf(1:j4)
 1530     format(/' * Error - (EQPT/pcraq) The species ',a,
     $    ' appearing',/7x,'on the data file in the ',a,' superblock',
     $    ' has fewer than',/7x,'two species in its associated',
     $    ' reaction. This is not a',/7x,'valid reaction.')
          nerr = nerr + 1
        endif
c
c       Read the names of the species in the reaction and the
c       corresponding reaction coefficients.
c
        n = 0
        do i = 1,ndrsts/2
          read (ndat0s,1000,end=990,err=995) uline
          ulbufa = uline
          do k = 1,2
            call lejust(ulbufa)
            ii = index(ulbufa,' ')
            if (ii .gt. 1) then
              ux16 = ulbufa(1:ii - 1)
              read (ux16,'(f10.4)',err=995) cdrsi(n + k)
            else
              cdrsi(n + k) = 0.
            endif
            ulbufb = ulbufa(ii + 1:80)
            call lejust(ulbufb)
            udrsi(n + k) = ulbufb(1:24)
            if (k .lt. 2) ulbufa = ulbufb(25:80)
          enddo
c         read (uline,1540,err=995) (cdrsi(n + k),udrsi(n + k), k = 1,2)
c1540     format((2(1x,f10.4,2x,a24)))
          n = n + 2
        enddo
c
        j = mod(ndrsts,2)
        if (j .gt. 0) then
          read (ndat0s,1000,end=990,err=995) uline
          ulbufa = uline
          do k = 1,j
            call lejust(ulbufa)
            ii = index(ulbufa,' ')
            if (ii .gt. 1) then
              ux16 = ulbufa(1:ii - 1)
              read (ux16,'(f10.4)',err=995) cdrsi(n + k)
            else
              cdrsi(n + k) = 0.
            endif
            ulbufb = ulbufa(ii + 1:80)
            call lejust(ulbufb)
            udrsi(n + k) = ulbufb(1:24)
            if (k .lt. j) ulbufa = ulbufb(25:80)
          enddo
c         read (uline,1540,err=995) (cdrsi(n + k),udrsi(n + k), k = 1,j)
          n = n + j
        endif
c
c       Check for blank species names, duplicate species names, and
c       zero-valued reaction coefficients. Check that the reaction
c       coefficient of the species with which the reaction is
c       associated has a negative value.
c
        call rxnsck(nbtmx1,cdrsi,nct,ndrsts,nentri,nerr,noutpt,
     $  ns,nsb,nttyo,qblkrs,qzerrs,udrsi,usblkf,uspec)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c       Read the log K grid for the current species.
c       Return the data in the xdbval holding array.
c
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $  udbval,xdbval)
        if (qend) go to 990
        if (qerr) go to 995
c
c       Load the data into the xlks array.
c
c       Calling sequence substitutions:
c         ns for ipc
c         nbtmx1 for ipcmax
c         xlks for zdbval
c
        call ldbar3(ns,nbtmx1,nacdpr,narxmx,narxt,ndbmax,
     $  ntprmx,ntprt,xdbval,xlks)
c
        if (itgenf .ge. 0) then
c
c         Test the grid ranges for sparseness of actual data.
c
          ustrgr = 'the log K for ' // uspec(ns)
          call tegrid(itgenf,nacdpr,narxt,nerr,noutpt,ntprmx,
     $    ntprt,nttyo,nwarn,ustrgr)
        endif
      endif
c
c     In the case of a strict basis species (a species which has no
c     real associated reaction), xhfs is the partial molar enthalpy
c     of formation and xvfs is the partial molar volume. In the case
c     of all other species, there is a real associated reaction and
c     xhfs is then the partial molar enthalpy of reaction and xvfs
c     is the partial molar volume of reaction.
c
      if (ipch .ge. 0) then
c
c       Read the enthalpy function grid for the current species.
c       Return the data in the xdbval holding array.
c
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $  udbval,xdbval)
        if (qend) go to 990
        if (qerr) go to 995
c
c       Load the data into the xhfs array.
c
c       Calling sequence substitutions:
c         ns for ipc
c         nbtmx1 for ipcmax
c         xhfs for zdbval
c
        call ldbar3(ns,nbtmx1,nacdpr,narxmx,narxt,ndbmax,
     $  ntprmx,ntprt,xdbval,xhfs)
c
        if (itgenf .ge. 0) then
c
c         Test the grid ranges for sparseness of actual data.
c
          ustrgr = 'the enthalpy function for ' // uspec(ns)
          call tegrid(itgenf,nacdpr,narxt,nerr,noutpt,ntprmx,
     $    ntprt,nttyo,nwarn,ustrgr)
        endif
c
        do ipc = 1,ipch
c
c         Read the enthalpy function derivative grid for the current
c         species. Return the data in the xdbval holding array.
c
          call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $    udbval,xdbval)
          if (qend) go to 990
          if (qerr) go to 995
c
c         Load the data into the dhfs array.
c
c         Calling sequence substitutions:
c           ipchmx for ipcmax
c           dhfs for zdbval
c
          call ldbar4(ipc,ipchmx,nacdpr,narxmx,narxt,nbtmx1,
     $    ndbmax,ns,ntprmx,ntprt,xdbval,dhfs)
        enddo
      endif
c
      if (ipcv .ge. 0) then
c
c       Read the volume function grid for the current species.
c       Return the data in the xdbval holding array.
c
        call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $  udbval,xdbval)
        if (qend) go to 990
        if (qerr) go to 995
c
c       Load the data into the xvfs array.
c
c       Calling sequence substitutions:
c         ns for ipc
c         nbtmx1 for ipcmax
c         xvfs for zdbval
c
        call ldbar3(ns,nbtmx1,nacdpr,narxmx,narxt,ndbmax,
     $  ntprmx,ntprt,xdbval,xvfs)
c
        if (itgenf .ge. 0) then
c
c         Test the grid ranges for sparseness of actual data.
c
          ustrgr = 'the volume function for ' // uspec(ns)
          call tegrid(itgenf,nacdpr,narxt,nerr,noutpt,ntprmx,
     $    ntprt,nttyo,nwarn,ustrgr)
        endif
c
        do ipc = 1,ipcv
c
c         Read the volume function derivative grid for the current
c         species. Return the data in the xdbval holding array.
c
          call rdgrid(ndat0s,ndbmax,ndbptg,ndbptl,qend,qerr,q500nd,
     $    udbval,xdbval)
          if (qend) go to 990
          if (qerr) go to 995
c
c         Load the data into the dvfs array.
c
c         Calling sequence substitutions:
c           ipcvmx for ipcmax
c           dvfs for zdbval
c
          call ldbar4(ipc,ipcvmx,nacdpr,narxmx,narxt,nbtmx1,
     $    ndbmax,ns,ntprmx,ntprt,xdbval,dvfs)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the redox species in the strict basis set. Set the
c     appropriate flags.
c
      if (ns .eq. nsb) then
        qelect = uspec(ns)(1:3).eq.'e- '
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     If the strict basis redox species is e-, change it to O2(g).
c
      if (ns .eq. nsb) then
        if (qelect) then
          uspec(ns) = 'O2(g)'
          ncts = 1
          uessi(1) = 'O'
          cessi(1) = 2
          zchar(ns) = 0.
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     If the reaction is written in terms of e-, rewrite it in terms
c     of O2(g).
c
      if (qelect) then
        if (ns .gt. nsb) then
          call etoo2(cdrsi,dhfe,dhfs,dvfe,dvfs,ipch,ipchmx,ipcv,
     $    ipcvmx,narxmx,narxt,nbtmx1,ndrsts,ns,ntprmx,ntprt,udrsi,
     $    xhfe,xhfs,xlke,xlks,xvfe,xvfs)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Process the species data.
c
      mwtss = 0.
      qnofes = .false.
      do n = 1,ncts
c
c       Search for element name in the uelem array.
c
        ux8 = uessi(n)(1:8)
        if (ux8(1:7) .ne. '<blank>') then
c
          do nc = 1,nct
            if (ux8(1:8) .eq. uelem(nc)(1:8)) then
              cess(nc,ns) = cess(nc,ns) + cessi(n)
              mwtss = mwtss + atwt(nc)*cessi(n)
              go to 220
            endif
          enddo
c
c         Error, not found.
c
          if (ns .gt. nsb) then
            j2 = ilnobl(uspec(ns))
            j3 = ilnobl(ux8)
            write (noutpt,1760) ux8(1:j3),uspec(ns)(1:j2)
            write (nttyo,1760) ux8(1:j3),uspec(ns)(1:j2)
 1760       format(/' * Error - (EQPT/pcraq) Unrecognized chemical',
     $      ' element "',a,'" is listed',/7x,'in the composition of',
     $      ' the species ',a,'.')
            nerr = nerr + 1
            qnofes = .true.
          endif
        endif
  220   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (ndata1) uspec(ns),ncts,ndrsts
      write (ndat1f,1800) uspec(ns),ncts,ndrsts
 1800 format(a24,2x,2i5)
c
      write (ndata1) mwtss,zchar(ns)
      write (ndat1f,1810) mwtss,zchar(ns)
 1810 format(5x,f10.3,f5.0)
c
      if (ncts .gt. 0) then
        write (ndata1) (cessi(n),uessi(n), n = 1,ncts)
        n = 1
  230   continue
        if (n .eq. ncts) then
          j4 = ilnobl(uessi(n))
          write (ndat1f,1820) cessi(n),uessi(n)(1:j4)
 1820     format(1x,f10.4,2x,a)
          n = n + 1
        else
          j4 = ilnobl(uessi(n + 1))
          write (ndat1f,1830) cessi(n),uessi(n),
     $    cessi(n + 1),uessi(n + 1)(1:j4)
 1830     format(1x,f10.4,2x,a8,1x,f10.4,2x,a)
          n = n + 2
        endif
        if (n .le. ncts) go to 230
      endif
c
      if (ns .le. nsb) go to 270
c
      write (ndata1) (cdrsi(n),udrsi(n), n = 1,ndrsts)
      n = 1
c
  240 continue
      if (n .eq. ndrsts) then
        j4 = ilnobl(udrsi(n))
        write (ndat1f,1860) cdrsi(n),udrsi(n)(1:j4)
 1860   format(1x,f10.4,2x,a)
        n = n + 1
      else
        j4 = ilnobl(udrsi(n + 1))
        write (ndat1f,1870) cdrsi(n),udrsi(n),
     $  cdrsi(n + 1),udrsi(n + 1)(1:j4)
 1870   format(1x,f10.4,2x,a24,1x,f10.4,2x,a)
        n = n + 2
      endif
      if (n .le. ndrsts) go to 240
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      qnofrs = .false.
      cdrs(nbt1,ns) = cdrsi(1)
      do n = 2,ndrsts
c
c       Search for species name in the uspec array.
c
        ux24 = udrsi(n)
        if (ux24(1:7) .ne. '<blank>') then
c
          do nse = 1,nbt
            if (ux24 .eq. uspec(nse)) then
              cdrs(nse,ns) = cdrs(nse,ns) + cdrsi(n)
              go to 250
            endif
          enddo
c
c         Error, not found.
c
          j2 = ilnobl(uspec(ns))
          j3 = ilnobl(ux24)
          write (noutpt,1880) uspec(ns)(1:j2),ux24(1:j3)
          write (nttyo,1880) uspec(ns)(1:j2),ux24(1:j3)
 1880     format(/' * Error - (EQPT/pcraq) The reaction which destroys',
     $    /7x,'non-basis species ',a,' is written in terms of an',
     $     /7x,'unrecognized basis species "',a,'".')
          nerr = nerr + 1
          qnofrs = .true.
  250     continue
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Test the reaction for mass and charge balance.
c     Skip if there are already obvious problems.
c
      if (ndrsts .ge. 2) then
        if (.not.qblkes .and. .not.qnofes) then
          if (.not.qblkrs .and. .not.qnofrs) then
            call rxnchk(cdrs,cess,mtotr,nbt,nbtmx1,nbtmx2,nco,
     $      nct,nctmax,nerr,noutpt,ns,nsb,nttyo,uelem,uspec,zchar)
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Fit interpolating polynomials to the log K grid.
c
      do ntpr = 1,ntprt
        do n = 1,narxt(ntpr)
          avgrid(n,ntpr) = xlks(n,ntpr,ns)
        enddo
      enddo
c
      call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $ narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $ xvec,yvec)
c
      ux24 = 'Log K'
      j2 = ilnobl(ux24)
      write (ndata1) ux24
      write (ndat1f,1010) ux24(1:j2)
c
      do ntpr = 1,ntprt
        nt = narxt(ntpr)
        write (ndata1) (apr(i,ntpr), i = 1,nt)
        write (ndat1f,1890) (apr(i,ntpr), i = 1,nt)
 1890   format( 5(1pe16.9) )
      enddo
c
  270 continue
c
      if (ipch .ge. 0) then
c
c       Fit interpolating polynomials to the enthalpy function grid.
c
        do ntpr = 1,ntprt
          do n = 1,narxt(ntpr)
            avgrid(n,ntpr) = xhfs(n,ntpr,ns)
          enddo
        enddo
c
        call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $  narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $  xvec,yvec)
c
        ux24 = 'Enthalpy'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
c
        do ntpr = 1,ntprt
          nt = narxt(ntpr)
          write (ndata1) (apr(i,ntpr), i = 1,nt)
          write (ndat1f,1890) (apr(i,ntpr), i = 1,nt)
        enddo
c
c       Fit interpolating polynomials to the grids for its
c       pressure derivatives.
c
        do ipc = 1,ipch
c
          do ntpr = 1,ntprt
            do n = 1,narxt(ntpr)
              avgrid(n,ntpr) = dhfs(n,ntpr,ipc,ns)
            enddo
          enddo
c
          call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $    narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $    xvec,yvec)
c
          ux24 = 'dhfs( )'
          write (ux24(6:6),'(i1)') ipc
          j2 = ilnobl(ux24)
          write (ndata1) ux24
          write (ndat1f,1010) ux24(1:j2)
c
          do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1890) (apr(i,ntpr), i = 1,nt)
          enddo
        enddo
      endif
c
      if (ipcv .ge. 0) then
c
c       Fit interpolating polynomials to the volume function grid.
c
        do ntpr = 1,ntprt
          do n = 1,narxt(ntpr)
            avgrid(n,ntpr) = xvfs(n,ntpr,ns)
          enddo
        enddo
c
        call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $  narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $  xvec,yvec)
c
        ux24 = 'Volume'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
c
        do ntpr = 1,ntprt
          nt = narxt(ntpr)
          write (ndata1) (apr(i,ntpr), i = 1,nt)
          write (ndat1f,1890) (apr(i,ntpr), i = 1,nt)
        enddo
c
c       Fit interpolating polynomials to the grids for its
c       pressure derivatives.
c
        do ipc = 1,ipcv
c
          do ntpr = 1,ntprt
            do n = 1,narxt(ntpr)
              avgrid(n,ntpr) = dvfs(n,ntpr,ipc,ns)
            enddo
          enddo
c
          call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $    narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $    xvec,yvec)
c
          ux24 = 'dvfs( )'
          write (ux24(6:6),'(i1)') ipc
          j2 = ilnobl(ux24)
          write (ndata1) ux24
          write (ndat1f,1010) ux24(1:j2)
c
          do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1890) (apr(i,ntpr), i = 1,nt)
          enddo
        enddo
      endif
c
      write (ndat1f,1010) utermc(1:72)
c
c     Look for another species block or terminator line.
c
      read (ndat0s,1000,end=990,err=995) uline
      nn = nn + 1
      go to 100
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 write (noutpt,2000)
      write (nttyo,2000)
 2000 format(/' * Error - (EQPT/pcraq) Unexpectedly encountered',
     $ /7x,'end-of-file while reading the DATA0 file.')
c
      write (noutpt,2010) nn
      write (nttyo,2010) nn
 2010 format(7x,'The value of the species block counter is ',i3,'.')
      if (ns .gt. 0) then
        j2 = ilnobl(uspec(ns))
        if (j2 .gt. 0) then
          write (noutpt,2020) uspec(ns)(1:j2)
          write (nttyo,2020) uspec(ns)(1:j2)
 2020     format(7x,'The last species name read was "',a,'".')
        else
          nsm1 = ns - 1
          if ((nsm1) .gt. 0) then
            j2 = ilnobl(uspec(nsm1))
            if (j2 .gt. 0) then
              write (noutpt,2020) uspec(nsm1)(1:j2)
              write (nttyo,2020) uspec(nsm1)(1:j2)
            endif
          endif
        endif
      endif
      j2 = ilnobl(uline)
      if (j2 .gt. 0) then
        j2 = min(j2,70)
        write (noutpt,2030) uline(1:j2)
        write (nttyo,2030) uline(1:j2)
 2030   format(7x,'The last line read was the following:',
     $  /7x,'"',a,'"')
      endif
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  995 write (noutpt,2040)
      write (nttyo,2040)
 2040 format(/' * Error - (EQPT/pcraq) Encountered a read format',
     $ /7x,'error while reading the DATA0 file.')
c
      write (noutpt,2010) nn
      write (nttyo,2010) nn
      if (ns .gt. 0) then
        j2 = ilnobl(uspec(ns))
        if (j2 .gt. 0) then
          write (noutpt,2020) uspec(ns)(1:j2)
          write (nttyo,2020) uspec(ns)(1:j2)
        else
          nsm1 = ns - 1
          if ((nsm1) .gt. 0) then
            j2 = ilnobl(uspec(nsm1))
            if (j2 .gt. 0) then
              write (noutpt,2020) uspec(nsm1)(1:j2)
              write (nttyo,2020) uspec(nsm1)(1:j2)
            endif
          endif
        endif
      endif
      j2 = ilnobl(uline)
      if (j2 .gt. 0) then
        j2 = min(j2,67)
        write (noutpt,2030) uline(1:j2)
        write (nttyo,2030) uline(1:j2)
      endif
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
