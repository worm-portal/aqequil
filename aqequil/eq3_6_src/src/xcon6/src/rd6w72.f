      subroutine rd6w72(cdac,cesrb,cplim,csigma,dlzmx1,dlzmx2,dlzidp,
     $ dzpllg,dzplot,dzprlg,dzprnt,eact,electr,fk,iact,ifile,iktbt,
     $ iktmax,imchmx,imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,
     $ jreac,jtemp,jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksplmx,ksppmx,
     $ kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,
     $ ndact,ndctmx,nesrbt,nffg,nffgmx,ninpts,nmodl1,nmodl2,nodbmx,
     $ nopgmx,noprmx,noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,
     $ nrctmx,nrk,nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,
     $ ntrymx,nttkmx,nttyo,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,
     $ nxrtmx,qend,qrderr,rk0,rxbarb,sscrew,sk,tempci,tempc0,timemx,
     $ tolbt,toldl,tolsat,tolsst,tolx,trk0,tstrt,ttk,udac,uelemb,
     $ uendb,uesrb,uffg,undms,unrms,uprs,ureac,utitl1,utitl2,uxct16,
     $ uxmd24,uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,zistrt,zkfac,
     $ zklogl,zklogu,zvclgi)
c
c     This subroutine reads the EQ6 input file in compact ("W") format
c     for version 7.2.
c
c     This subroutine is called by:
c
c       XCON6/xcon6.f
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer iktmax,imchmx,kmax,nctmax,ndctmx,nffgmx,nodbmx,nopgmx,
     $ noprmx,noptmx,nprsmx,nrctmx,nsrtmx,nsscmx,ntitmx,nttkmx,nxmdmx,
     $ nxopmx,nxpemx,nxrtmx
c
      integer iact(imchmx,2,nrctmx),iktbt(nxrtmx),imech(2,nrctmx),
     $ iodb(nodbmx),iopg(nopgmx),iopr(noprmx),iopt(noptmx),
     $ jcode(nrctmx),jreac(nrctmx),jxmod(nxmdmx),kxmod(nxmdmx),
     $ ndact(imchmx,2,nrctmx),nesrbt(nsrtmx),nrk(2,nrctmx),nsk(nrctmx)
c
      integer ifile,ioscan,itermx,jtemp,kct,kdim,ksq,kmt,kprs,ksplmx,
     $ ksppmx,kstpmx,kxt,nffg,ninpts,nmodl1,nmodl2,nordlm,npslmx,nprmn,
     $ nprmx,nrct,nsslmx,ntitl1,ntitl2,ntrymx,nttyo,nxmod,nxopex,nxopt
c
      logical qend,qrderr
c
      character*80 utitl1(ntitmx),utitl2(ntitmx)
      character*48 uprs(nprsmx)
      character*24 udac(ndctmx,imchmx,2,nrctmx),uendb(iktmax,nxrtmx),
     $ uffg(nffgmx),undms(kmax),unrms(kmax),ureac(nrctmx),
     $ uxmd24(nxmdmx),uxopex(nxpemx)
      character*16 uxct16(nxopmx)
      character*8 uelemb(nctmax),uesrb(nctmax,nsrtmx),uxopt(nxopmx)
c
      real*8 cdac(ndctmx,imchmx,2,nrctmx),cesrb(nctmax,nsrtmx),
     $ csigma(imchmx,2,nrctmx),eact(imchmx,2,nrctmx),fk(nrctmx),
     $ hact(imchmx,2,nrctmx),modr(nrctmx),moffg(nffgmx),morr(nrctmx),
     $ mprs(nprsmx),mteaqb(nctmax),mteb(nctmax),rk0(imchmx,2,nrctmx),
     $ rxbarb(iktmax,nxrtmx),sscrew(nsscmx),sk(nrctmx),ttk(nttkmx),
     $ trk0(imchmx,2,nrctmx),vreac(nrctmx),xlkffg(nffgmx),
     $ xlkmod(nxmdmx),zvclgi(kmax)
c
      real*8 cplim,dlzmx1,dlzmx2,dlzidp,dzpllg,dzplot,dzprlg,dzprnt,
     $ electr,tempci,tempc0,timemx,tolbt,toldl,tolsat,tolsst,tolx,
     $ tstrt,zkfac,zklogl,zklogu,zimax,zistrt
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,iktb,j,kc,kcol,n,ncb,nrc,nsrt,nxrt
c
      character*80 uline
      character*24 ux,uxx
      character*8 uendit
c
      real*8 xx
c
c-----------------------------------------------------------------------
c
      data uendit /'endit.  '/
c
c-----------------------------------------------------------------------
c
      qrderr = .false.
c
c-----------------------------------------------------------------------
c
c     Main title.
c
c     Note: if there are exactly ntitpa lines in the title, the title
c     need not be terminated by an 'endit.'. The 'endit.' if present
c     is here not considered to be part of the title. The secondary
c     title is treated in the same manner.
c
      qend = .false.
      read (ninpts,1100,end=100,err=990) uline
 1100 format(a80)
      go to 105
c
  100 qend = .true.
      go to 999
c
  105 n = 1
      if (uline(1:8) .eq. uendit(1:8)) go to 120
      utitl1(1) = uline
c
      do 110 n = 2,ntitmx
        read (ninpts,1100,err=990) uline
        if (uline(1:8) .eq. uendit(1:8)) go to 120
        utitl1(n) = uline
  110 continue
      n = n + 1
c
      read (ninpts,1100,err=990) uline
      call locase(uline)
      j = index(uline,'nmodl1=')
      if (j .gt. 0) then
        backspace(ninpts)
      else
        write (nttyo,1015) ntitmx
 1015   format(/' * Error - (XCON6/rd6w72) Have too many lines in the',
     $  /7x,'main title. The code is only dimensioned for ',i4,
     $  /7x,'lines. Reduce the size of the title or increase the',
     $  /7x,'dimensioning parameter ntitpa.')
        go to 990
      endif
c
  120 ntitl1 = n - 1
c
c-----------------------------------------------------------------------
c
c     Nmodl option switches.
c
      read (ninpts,1110,err=990) nmodl1,nmodl2
 1110 format(12x,i2,22x,i2)
c
c     Temperature parameters.
c     Note: ttk(1) = tk1, etc.
c
      read (ninpts,1130,err=990) tempc0,jtemp,(ttk(i), i = 1,3)
 1130 format(12x,e12.5,12x,i2,/3(12x,e12.5))
c
c      Zi and time parameters.
c
      read (ninpts,1150,err=990) zistrt,zimax,tstrt,timemx,
     $ kstpmx,cplim
 1150 format(2(12x,e12.5),/2(12x,e12.5),/12x,i12,12x,e12.5)
c
c     Print interval parameters.
c
      read (ninpts,1170,err=990) dzprnt,dzprlg,ksppmx
 1170 format (2(12x,e12.5),12x,i5)
c
c     Plot interval parameters.
c
      read (ninpts,1190,err=990) dzplot,dzpllg,ksplmx
 1190 format (2(12x,e12.5),12x,i5)
c
c     Ifile.
c
      read (ninpts,1210,err=990) ifile
 1210 format(12x,i2,22x,i2)
c
c-----------------------------------------------------------------------
c
c     Iopt option switches.
c     Note: iopt(1) = iopt1, etc.
c
      read (ninpts,1240,err=990) (iopt(i), i = 1,20)
 1240 format(12x,10i5)
c
c     Iopr option switches.
c     Note: iopr(1) = iopr1, etc.
c
      read (ninpts,1260,err=990) (iopr(i), i = 1,20)
 1260 format(12x,10i5)
c
c     Iodb option switches.
c     Note: iodb(1) = iodb1, etc.
c
      read (ninpts,1280,err=990) (iodb(i), i = 1,20)
 1280 format(12x,10i5)
c
c-----------------------------------------------------------------------
c
c     Number of nxopt options.
c
      read (ninpts,1300,err=990) nxopt
 1300 format(12x,i2)
c
      if (nxopt .gt. nxopmx) then
        write (nttyo,1315) nxopmx
 1315   format(/' * Error - (XCON6/rd6w72) Have too many mineral',
     $  /7x,'subset-selection suppression options. The code is',
     $  /7x,'only dimensioned for ',i3,' such options. Reduce the',
     $  /7x,'number of options or increase the dimensioning',
     $  /7x,'parameter nxoppa.')
        go to 990
      endif
c
c     Nxopt options.
c
      if (nxopt .gt. 0) then
        do 140 n = 1,nxopt
          read (ninpts,1320,err=990) uxopt(n),uxct16(n)
 1320     format(12x,a6,1x,a8)
  140   continue
        read (ninpts,1340,err=990) nxopex
 1340   format(12x,i2)
c
        if (nxopex .gt. nxpemx) then
          write (nttyo,1345) nxpemx
 1345     format(/' * Error - (XCON6/rd6w72) Have too many',
     $    /7x,'exceptions specified to the mineral subset-selection',
     $    /7x,'suppression options. The code is only dimensioned',
     $    /7x,'for ',i3,'exceptions. Reduce the number of exceptions',
     $    /7x,'or increase the dimensioning parameter nxoppa.')
          go to 990
        endif
c
        if (nxopex .gt. 0) then
          do 150 n = 1,nxopex
            read (ninpts,1360,err=990) uxopex(n)
 1360       format(12x,a24)
  150     continue
        endif
      endif
c
c-----------------------------------------------------------------------
c
c     Number of nffg options.
c
      read (ninpts,1380,err=990) nffg
 1380 format(12x,i2)
c
      if (nffg .gt. nffgmx) then
        write (nttyo,1385) nffgmx
 1385   format(/' * Error - (XCON6/rd6w72) Have too many gases whose',
     $  /7x,'fugacities are to be fixed. The code is only dimensioned',
     $  /7x,'for ',i4,' such gases. Reduce the number of gases or',
     $  /7x,'increase the dimensioning parameter nffgpa.')
        go to 990
      endif
c
c     Nffg options.
c
      if (nffg .gt. 0) then
      do 160 n = 1,nffg
        read (ninpts,1400,err=990) uffg(n),moffg(n),xlkffg(n)
 1400   format(12x,a24,/12x,e12.5,12x,e12.5)
  160   continue
      endif
c
c-----------------------------------------------------------------------
c
c     Number of reactants.
c
      read (ninpts,1420,err=990) nrct
 1420 format(12x,i2)
c
      if (nrct .gt. nrctmx) then
        write (nttyo,1425) nrctmx
 1425   format(/' * Error - (XCON6/rd6w72) Have too many reactants',
     $  /7x,'The code is only dimensioned for ',i4,' reactants.',
     $  /7x,'Reduce the number of reactants or increase the',
     $  /7x,'dimensioning parameter nrctpa.')
        go to 990
      endif
c
c     Reactants.
c
      nsrt = 0
      nxrt = 0
      if (nrct .gt. 0) then
        do 390 nrc = 1,nrct
c
c         Name, flags, and masses.
c
          read (ninpts,1440,err=990) ureac(nrc),jcode(nrc),jreac(nrc),
     $    morr(nrc),modr(nrc)
 1440     format(12x,a24,/12x,i2,22x,i2,/2(12x,e12.5))
c
          if (jcode(nrc) .eq. 1) then
c
c           Solid solution compositions.
c
            nxrt = nxrt + 1
c
            if (nxrt .gt. nxrtmx) then
              write (nttyo,1445) nxrtmx
 1445         format(/' * Error - (XCON6/rd6w72) Have too many solid',
     $        /7x,'solution reactants. The code is only dimensioned',
     $        /7x,'for ',i4,' such reactants. Reduce the number of',
     $        /7x,'such reactants or increase the dimensioning',
     $        ' parameter nxrtpa.')
              go to 990
            endif
c
            iktb = 0
            do 170 i = 1,iktmax + 1
              read (ninpts,1460,err=990) ux,xx
 1460         format(3x,a24,1x,e12.5)
              if (ux(1:8) .eq. uendit(1:8)) go to 180
c
              iktb = iktb + 1
c
              if (iktb .gt. iktmax) then
                write (nttyo,1446) ureac(nrc),iktmax
 1446           format(/' * Error - (XCON6/rd6w72) Have too many',
     $          ' end-members',/7x,'in the solid solution reactant',
     $          ' "',a24,'".',/7x,'The code is only dimensioned for ',
     $          i4,' end-members per',/7x,'solid solution. Reduce',
     $          ' the number of end-members or',
     $          /7x,'increase the dimensioning parameter iktpar.')
                go to 990
              endif
c
              uendb(iktb,nxrt) = ux
              rxbarb(iktb,nxrt) = xx
  170       continue
  180       iktbt(nxrt) = iktb
          elseif (jcode(nrc) .eq. 2) then
c
c           Special reactant compositions.
c
            nsrt = nsrt + 1
c
            if (nsrt .gt. nsrtmx) then
              write (nttyo,1447) nsrtmx
 1447         format(/' * Error - (XCON6/rd6w72) Have too many special',
     $        /7x,'reactants. The code is only dimensioned for ',i4,
     $        /7x,'such reactants. Reduce the number of such reactants',
     $        /7x,'or increase the dimensioning parameter nsrtpa.')
              go to 990
            endif
c
            read (ninpts,1480,err=990) vreac(nrc)
 1480       format(12x,e12.5)
c
            ncb = 0
            do 190 n = 1,nctmax + 1
              read (ninpts,1500,err=990) ux,xx
 1500         format(3x,a8,1x,e22.15)
              if (ux(1:8) .eq. uendit(1:8)) go to 200
c
              ncb = ncb + 1
c
              if (ncb .gt. nctmax) then
                write (nttyo,1510) ureac(nrc),nctmax
 1510           format(/' * Error - (XCON6/rd6w72) Have too many',
     $          '  chemical',
     $          /7x,'elements in the special  reactant "',a24,'".',
     $          /7x,'The code is only dimensioned for ',i4,' elements.',
     $          /7x,'Reduce the number of elements or increase the',
     $          /7x,'dimensioning parameter nctpar.')
                go to 990
              endif
c
              uesrb(ncb,nsrt) = ux
              cesrb(ncb,nsrt) = xx
  190       continue
  200       nesrbt(nsrt) = ncb
          endif
c
c         Surface area parameters.
c
          read (ninpts,1520,err=990) nsk(nrc),sk(nrc),fk(nrc)
 1520     format(12x,i2,22x,e12.5,12x,e12.5)
c
c         Rate law codes.
c
          read (ninpts,1540,err=990) nrk(1,nrc),nrk(2,nrc)
 1540     format(12x,i2,22x,i2)
c
c         Rate law parameters, forward direction.
c
          if (nrk(1,nrc) .eq. 1) then
            imech(1,nrc) = 3
c
            if (imech(1,nrc) .gt. imchmx) then
              write (nttyo,1550) ureac(nrc),imchmx
 1550         format(/' * Error - (XCON6/rd6w72) Have too many rate',
     $        /7x,'constants in the forward direction rate law for',
     $        /7x,'reactant "',a24,'". The code is only',
     $        /7x,'dimensioned for ',i2,' rate constants per rate law.',
     $        /7x,'Reduce the number of rate constants or increase the',
     $        /7x,'dimensioning parameter imchpa.')
              go to 990
            endif
c
            read (ninpts,1560,err=990) (rk0(j,1,nrc), j = 1,3)
 1560       format (3(12x,e12.5))
          elseif (nrk(1,nrc) .eq. 2) then
            read (ninpts,1580,err=990) imech(1,nrc)
 1580       format(12x,i2)
c
            if (imech(1,nrc) .gt. imchmx) then
              write (nttyo,1550) ureac(nrc),imchmx
              go to 990
            endif
c
            do 220 i = 1,imech(1,nrc)
              read (ninpts,1600,err=990) rk0(i,1,nrc),trk0(i,1,nrc),
     $        iact(i,1,nrc)
 1600         format(12x,e12.5,12x,e12.5,12x,i2)
              read (ninpts,1620,err=990) eact(i,1,nrc),hact(i,1,nrc)
 1620         format(12x,e12.5,12x,e12.5)
              read (ninpts,1650,err=990) ndact(i,1,nrc),csigma(i,1,nrc)
 1650         format(12x,i2,22x,e12.5)
c
              if (ndact(i,1,nrc) .gt. ndctmx) then
                write (nttyo,1655) ureac(nrc),i,ndctmx
 1655           format(/' * Error - (XCON6/rd6w72) Have too many',
     $          /7x,'species in the activity product in term ',i2,
     $          /7x,'of the forward direction rate law for reactant',
     $          /7x,'"',a24,'". The code is only dimensioned for ',i3,
     $          /7x,'such species. Reduce the number of such species',
     $          /7x,'or increase the dimensioning parameter ndctpa.')
                go to 990
              endif
c
              if (ndact(i,1,nrc) .gt. 0) then
                do 210 n = 1,ndact(i,1,nrc)
                  read (ninpts,1670,err=990) udac(n,i,1,nrc),
     $            cdac(n,i,1,nrc)
 1670             format(12x,a24,12x,e12.5)
  210           continue
              endif
  220       continue
          elseif (nrk(1,nrc) .eq. 3) then
            imech(1,nrc) = 1
c
            if (imech(1,nrc) .gt. imchmx) then
              write (nttyo,1550) ureac(nrc),imchmx
              go to 990
            endif
c
            i = 1
            read (ninpts,1600,err=990) rk0(i,1,nrc),trk0(i,1,nrc),
     $      iact(i,1,nrc)
            read (ninpts,1620,err=990) eact(i,1,nrc),hact(i,1,nrc)
          elseif (nrk(1,nrc) .eq. 4) then
            read (ninpts,1580,err=990) imech(1,nrc)
c
            if (imech(1,nrc) .gt. imchmx) then
              write (nttyo,1550) ureac(nrc),imchmx
              go to 990
            endif
c
            do 240 i = 1,imech(1,nrc)
              read (ninpts,1600,err=990) rk0(i,1,nrc),trk0(i,1,nrc),
     $        iact(i,1,nrc)
              read (ninpts,1620,err=990) eact(i,1,nrc),hact(i,1,nrc)
              read (ninpts,1650,err=990) ndact(i,1,nrc)
c
              if (ndact(i,1,nrc) .gt. ndctmx) then
                write (nttyo,1655) ureac(nrc),i,ndctmx
                go to 990
              endif
c
              if (ndact(i,1,nrc) .gt. 0) then
                do 230 n = 1,ndact(i,1,nrc)
                  read (ninpts,1670,err=990) udac(n,i,1,nrc),
     $            cdac(n,i,1,nrc)
  230           continue
              endif
  240       continue
          endif
c
c         Rate law parameters, backward direction.
c
          if (nrk(2,nrc) .eq. 1) then
            imech(2,nrc) = 3
c
            if (imech(2,nrc) .gt. imchmx) then
              write (nttyo,1673) ureac(nrc),imchmx
 1673         format(/' * Error - (XCON6/rd6w72) Have too many rate',
     $        /7x,'constants in the backward direction rate law for',
     $        /7x,'reactant "',a24,'". The code is only',
     $        /7x,'dimensioned for ',i2,' rate constants per rate law.',
     $        /7x,'Reduce the number of rate constants or increase the',
     $        /7x,'dimensioning parameter imchpa.')
              go to 990
            endif
c
            read (ninpts,1560,err=990) (rk0(j,2,nrc), j = 1,3)
          elseif (nrk(2,nrc) .eq. 2) then
            read (ninpts,1580,err=990) imech(2,nrc)
c
            if (imech(2,nrc) .gt. imchmx) then
              write (nttyo,1673) ureac(nrc),imchmx
              go to 990
            endif
c
            do 260 i = 1,imech(2,nrc)
              read (ninpts,1600,err=990) rk0(i,2,nrc),trk0(i,2,nrc),
     $        iact(i,2,nrc)
              read (ninpts,1620,err=990) eact(i,2,nrc),hact(i,2,nrc)
              read (ninpts,1650,err=990) ndact(i,2,nrc),csigma(i,2,nrc)
c
              if (ndact(i,2,nrc) .gt. ndctmx) then
                write (nttyo,1675) ureac(nrc),i,ndctmx
 1675           format(/' * Error - (XCON6/rd6w72) Have too many',
     $          /7x,'species in the activity product in term ',i2,
     $          /7x,'of the backward direction rate law for reactant',
     $          /7x,'"',a24,'". The code is only dimensioned for ',i3,
     $          /7x,'such species. Reduce the number of such species',
     $          /7x,'or increase the dimensioning parameter ndctpa.')
                go to 990
              endif
c
              if (ndact(i,2,nrc) .gt. 0) then
                do 250 n = 1,ndact(i,2,nrc)
                  read (ninpts,1670,err=990) udac(n,i,2,nrc),
     $            cdac(n,i,2,nrc)
  250           continue
              endif
  260       continue
          elseif (nrk(2,nrc) .eq. 3) then
            imech(2,nrc) = 1
            if (imech(2,nrc) .gt. imchmx) then
              write (nttyo,1673) ureac(nrc),imchmx
              go to 990
            endif
c
            i = 1
            read (ninpts,1600,err=990) rk0(i,2,nrc),trk0(i,2,nrc),
     $      iact(i,2,nrc)
            read (ninpts,1620,err=990) eact(i,2,nrc),hact(i,2,nrc)
          elseif (nrk(2,nrc) .eq. 4) then
            read (ninpts,1580,err=990) imech(2,nrc)
c
            if (imech(2,nrc) .gt. imchmx) then
              write (nttyo,1673) ureac(nrc),imchmx
              go to 990
            endif
c
            do 280 i = 1,imech(2,nrc)
              read (ninpts,1600,err=990) rk0(i,2,nrc),trk0(i,2,nrc),
     $        iact(i,2,nrc)
              read (ninpts,1620,err=990) eact(i,2,nrc),hact(i,2,nrc)
              read (ninpts,1650,err=990) ndact(i,2,nrc)
c
              if (ndact(i,2,nrc) .gt. ndctmx) then
                write (nttyo,1675) ureac(nrc),i,ndctmx
                go to 990
              endif
c
              if (ndact(i,2,nrc) .gt. 0) then
                do 270 n = 1,ndact(i,2,nrc)
                  read (ninpts,1670,err=990) udac(n,i,2,nrc),
     $            cdac(n,i,2,nrc)
  270           continue
              endif
  280       continue
          endif
  390   continue
      endif
c
c-----------------------------------------------------------------------
c
c     Dump interval.
c
      read (ninpts,2000,err=990) dlzidp
 2000 format(3(12x,e12.5))
c
c     Tolerances.
c
      read (ninpts,2020,err=990) tolbt,toldl,tolx,tolsat,tolsst
 2020 format(3(12x,e12.5))
c
c     Setscrew parameters.
c     Note: sscrew(1) = screw1, etc.
c
      read (ninpts,2040,err=990) (sscrew(i), i = 1,6)
 2040 format(3(12x,e12.5))
c
c     Z parameters.
c
      read (ninpts,2060,err=990) zklogu,zklogl,zkfac
 2060 format(3(12x,e12.5))
c
c     Step size limits and maximum order.
c
      read (ninpts,2080,err=990) dlzmx1,dlzmx2,nordlm
 2080 format (2(12x,e12.5),12x,i2)
c
c     Maximum number of iterations and maximum number of phase
c     assemblage tries.
c
      read (ninpts,2100,err=990) itermx,ntrymx
 2100 format(12x,i2,22x,i2)
c
c     Slide and scan control parameters.
c
      read (ninpts,2120,err=990) npslmx,nsslmx,ioscan
 2120 format(12x,i2,22x,i2,22x,i2)
c
c-----------------------------------------------------------------------
c
c     Process the bottom half of the current output file.
c
c     Secondary title.
c
      do 400 n = 1,ntitmx
        read (ninpts,3000,err=990) uline
 3000   format(a80)
        if (uline(1:8) .eq. uendit(1:8)) go to 410
        utitl2(n) = uline
  400 continue
      n = n + 1
c
      read (ninpts,1100,err=990) uline
      call locase(uline)
      j = index(uline,'tempci=')
      if (j .gt. 0) then
        backspace(ninpts)
      else
        write (nttyo,3002) ntitmx
 3002   format(/' * Error - (XCON6/rd6w72) Have too many lines in the',
     $  /7x,'secondary title. The code is only dimensioned for ',i4,
     $  /7x,'lines. Reduce the size of the title or increase the',
     $  /7x,'dimensioning parameter ntitpa.')
        go to 990
      endif
c
  410 ntitl2 = n - 1
c
c-----------------------------------------------------------------------
c
c     Original temperature.
c
      read (ninpts,3010,err=990) tempci
 3010 format(12x,e12.5)
c
c-----------------------------------------------------------------------
c
c     Number of nxmod options.
c
      read (ninpts,3030,err=990) nxmod
 3030 format(12x,i2)
c
      if (nxmod .gt. nxmdmx) then
        write (nttyo,3035) nxmdmx
 3035   format(/' * Error - (XCON6/rd6w72) Have too many nxmod',
     $  /7x,'alter/suppress options. The code is only dimensioned',
     $  /7x,'for ',i3,' such options. Reduce the number of such',
     $  ' options',/7x,'or increase the dimensioning parameter',
     $  ' nxmdpa.')
        go to 990
      endif
c
c     Nxmod options.
c
      if (nxmod .gt. 0) then
        do 420 n = 1,nxmod
          read (ninpts,3050,err=990) uxmd24(n),jxmod(n),kxmod(n),
     $    xlkmod(n)
 3050     format(12x,a24,/12x,i2,22x,i2,22x,e12.5)
  420   continue
      endif
c
c-----------------------------------------------------------------------
c
c     Iopg options.
c     Note: iopg(1) = iopg1, etc.
c
      read (ninpts,3070,err=990) (iopg(i), i = 1,10)
 3070 format (12x,i2,22x,i2,22x,i2)
c
c-----------------------------------------------------------------------
c
c     Index limits.
c
      read (ninpts,3090,err=990) kct,ksq,kmt,kxt,kdim,kprs
 3090 format(3(12x,i2,10x),/3(12x,i2,10x))
c
      if (kct .gt. nctmax) then
        write (nttyo,3095) nctmax
 3095   format(/' * Error - (XCON6/rd6w72) Have too many chemical',
     $  /7x,'elements present. The code is only dimensioned',
     $  /7x,'for ',i3,' elements. Reduce the number of elements',
     $  /7x,'or increase the dimensioning parameter nctpar.')
        go to 990
      endif
c
      if (kdim .gt. kmax) then
        write (nttyo,3100) kmax
 3100   format(/' * Error - (XCON6/rd6w72) Have too many master',
     $  /7x,'variables. The code is only dimensioned for ',i3,
     $  /7x,'master variables. Reduce the number of such variables',
     $  /7x,'or increase the dimensioning parameter kpar.')
        go to 990
      endif
c
c-----------------------------------------------------------------------
c
c     Balance totals.
c
      do 430 kc = 1,kct
        read (ninpts,3120,err=990) uelemb(kc),mteb(kc),mteaqb(kc)
 3120   format(3x,a8,1x,e22.15,1x,e22.15)
  430 continue
      read (ninpts,3140,err=990) ux,electr
 3140 format(3x,a8,1x,e22.15)
c
c-----------------------------------------------------------------------
c
c     Basis variable data.
c
      do 440 kcol = 1,kdim
        read (ninpts,3160,err=990) unrms(kcol),undms(kcol),zvclgi(kcol)
 3160   format(3x,a24,1x,a24,1x,e22.15)
  440 continue
c
c-----------------------------------------------------------------------
c
c     Physically removed system data.
c
      nprmn = 0
      nprmx = 0
      if (kprs .gt. 0) then
c
        do 450 n = 1,nprsmx
          read (ninpts,3180,err=990) ux,xx
 3180     format(3x,a24,1x,e22.15)
          if (ux(1:8) .eq. uendit(1:8)) go to 460
          nprmn = nprmn + 1
c
          if (nprmn .gt. nprsmx) then
            write (nttyo,3190) nprsmx
 3190       format(/' * Error - (XCON6/rd6w72) Have too many mineral',
     $      /7x,'species in the physically removed system. The code is',
     $      /7x,'only dimensioned for ',i3,' such species. Reduce the',
     $      /7x,'number of such species or increase the dimensioning',
     $      /7x,'parameter nprspa.')
            go to 990
          endif
c
          uprs(nprmn)(1:24) = ux
          uprs(nprmn)(25:48) = ux
          mprs(nprmn) = xx
  450   continue
  460   continue
        nprmx = nprmn
c
        do 480 n = 1,nprsmx
          read (ninpts,3200,err=990) ux
 3200     format(3x,a24)
          if (ux(1:8) .eq. uendit(1:8)) go to 490
          do 470 i = 1,iktmax + 1
            read (ninpts,3220,err=990) uxx,xx
 3220       format(3x,a24,1x,e22.15)
            if (uxx(1:8) .eq. uendit(1:8)) go to 480
            nprmx = nprmx + 1
c
            if (nprmn .gt. nprsmx) then
              write (nttyo,3190) nprsmx
              go to 990
            endif
c
            uprs(nprmx)(1:24) = ux
            uprs(nprmx)(25:48) = uxx
            mprs(nprmx) = xx
  470     continue
  480   continue
  490   continue
      endif
c
      go to 999
c
  990 qrderr = .true.
c
  999 continue
      end
