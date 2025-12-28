      subroutine wr6d7(cdac,cesrb,csigma,dlzmx1,dlzmx2,dlzidp,
     $ dzprlg,dzprnt,eact,electr,fk,iact,iktbt,iktmax,imchmx,
     $ imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,jreac,jtemp,
     $ jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksppmx,kstpmx,kxmod,kxt,
     $ hact,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,ndact,ndctmx,
     $ nesrbt,newin,nffg,nffgmx,nmodl1,nmodl2,nodbmx,nopgmx,noprmx,
     $ noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,nrctmx,nrk,
     $ nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,ntrymx,nttkmx,
     $ nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrtmx,rk0,rxbarb,
     $ sscrew,sk,tempci,tempc0,timemx,tolbt,toldl,tolsat,tolsst,
     $ tolx,trk0,tstrt,ttk,ucode,udac,uelemb,uendb,uesrb,ueqlrn,
     $ ueqlst,uffg,undms,unrms,uprs,ureac,urelno,ustage,utitl1,
     $ utitl2,uxct16,uxmd24,uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,
     $ zistrt,zkfac,zklogl,zklogu,zvclgi)
c
c     This subroutine writes the EQ6 input file in menu-style ("D")
c     format for versions 7.0-7.1.
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
      integer ioscan,itermx,jtemp,kct,kdim,ksq,kmt,kprs,ksppmx,
     $ kstpmx,kxt,newin,nffg,nmodl1,nmodl2,nordlm,npslmx,nprmn,nprmx,
     $ nrct,nsslmx,ntitl1,ntitl2,ntrymx,nxmod,nxopex,nxopt
c
      character*80 utitl1(ntitmx),utitl2(ntitmx)
      character*48 uprs(nprsmx)
      character*24 udac(ndctmx,imchmx,2,nrctmx),uendb(iktmax,nxrtmx),
     $ uffg(nffgmx),undms(kmax),unrms(kmax),ureac(nrctmx),
     $ uxmd24(nxmdmx),uxopex(nxpemx)
      character*16 uxct16(nxopmx)
      character*8 uelemb(nctmax),uesrb(nctmax,nsrtmx),uxopt(nxopmx)
      character*8 ucode,urelno,ustage,ueqlrn,ueqlst
c
      real*8 cdac(ndctmx,imchmx,2,nrctmx),cesrb(nctmax,nsrtmx),
     $ csigma(imchmx,2,nrctmx),eact(imchmx,2,nrctmx),fk(nrctmx),
     $ hact(imchmx,2,nrctmx),modr(nrctmx),moffg(nffgmx),morr(nrctmx),
     $ mprs(nprsmx),mteaqb(nctmax),mteb(nctmax),rk0(imchmx,2,nrctmx),
     $ rxbarb(iktmax,nxrtmx),sscrew(nsscmx),sk(nrctmx),ttk(nttkmx),
     $ trk0(imchmx,2,nrctmx),vreac(nrctmx),xlkffg(nffgmx),
     $ xlkmod(nxmdmx),zvclgi(kmax)
c
      real*8 dlzmx1,dlzmx2,dlzidp,dzprlg,dzprnt,electr,tempci,tempc0,
     $ timemx,tolbt,toldl,tolsat,tolsst,tolx,tstrt,zkfac,zklogl,
     $ zklogu,zimax,zistrt
c
c-----------------------------------------------------------------------
c
      include 'xcon6/x6op7.h'
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,idesc,idescx,iktb,ival,j2,k,kc,kcol,n,ncb,npr,nrc,
     $ nsrt,nxrt
c
      integer ilnobl
c
      character*24 ux24(3)
      character*24 ublk24,unone,uxx
      character*16 udef(nto6pa),ujx(0:3),umint(0:4)
      character*8 usup(-1:2),ux8(3)
      character*8 ublk8
      character*1 ux1(3)
c
      real*8 zero
c
c-----------------------------------------------------------------------
c
      data ublk24 /'                        '/
      data unone  /'none                    '/
c
      data udef(1)  / '40      itermx  '/
      data udef(2)  / 'varies  dlzidp  '/
      data udef(3)  / '1.0e-06 tolbt   '/
      data udef(4)  / '1.0e-06 toldl   '/
      data udef(5)  / 'varies  tolx    '/
      data udef(6)  / 'varies  tolsat  '/
      data udef(7)  / 'varies  tolsst  '/
      data udef(8)  / '1.0e-04 screw1  '/
      data udef(9)  / 'n/a     screw2  '/
      data udef(10) / '1.0e-04 screw3  '/
      data udef(11) / '1.0e-04 screw4  '/
      data udef(12) / '4.0     screw5  '/
      data udef(13) / '4.0     screw6  '/
      data udef(14) / 'varies  zklogu  '/
      data udef(15) / '2.0     zklogl  '/
      data udef(16) / '.98     zkfac   '/
      data udef(17) / 'varies  dlzmx1  '/
      data udef(18) / 'varies  dlzmx2  '/
      data udef(19) / '6       nordlm  '/
      data udef(20) / '25      ntrymx  '/
      data udef(21) / '8       npslmx  '/
      data udef(22) / '3       nsslmx  '/
      data udef(23) / 'none    ioscan  '/
c
      data ujx(0)   / 'aqueous         '/
      data ujx(1)   / 'mineral         '/
      data ujx(2)   / 'gas             '/
      data ujx(3)   / 'solid solution  '/
c
      data umint(0) / 'mineral         '/
      data umint(1) / 'solid solution  '/
      data umint(2) / 'special         '/
      data umint(3) / 'aqueous         '/
      data umint(4) / 'gas             '/
c
      data usup(-1) / 'suppress'/
      data usup(0)  / 'replace '/
      data usup(1)  / 'augmentk'/
      data usup(2)  / 'augmentg'/
c
      data ublk8  /'        '/
c
      data zero   /0./
c
c-----------------------------------------------------------------------
c
c     Write the new input file.
c
c     Title.
c
      write (newin,1090)
 1090 format('|',70('-'),'|')
      do 110 n = 1,ntitl1
        write (newin,1100) utitl1(n)
 1100   format('|',a70,'|')
  110 continue
      write (newin,1090)
c
c-----------------------------------------------------------------------
c
c     Nmodl option switches.
c
      ux1(1) = ' '
      ux1(2) = ' '
      ux1(3) = ' '
      if (nmodl2 .eq. 0) then
        ux1(1) = '*'
      elseif (nmodl2 .eq. 1) then
        ux1(2) = '*'
      elseif (nmodl2.eq.2) then
        ux1(3) = '*'
      else
        ux1(1) = '*'
      endif
      write (newin,1120) ux1(1),ux1(2),ux1(3)
 1120 format('| calculational mode   |',a1,'normal',t37,'|',a1,
     $ 'economy',t57,'|',a1,'super economy',t72,'|')
c
      write (newin,1090)
      ux1(1) = ' '
      ux1(2) = ' '
      ux1(3) = ' '
      if (nmodl1 .eq. 1) then
        ux1(1) = '*'
      elseif (nmodl1.eq.2 .or. nmodl1.eq.0) then
        ux1(2) = '*'
      elseif (nmodl1 .eq. 3) then
        ux1(3) = '*'
      else
        ux1(1) = '*'
      endif
      write (newin,1130) ux1(1),ux1(2),ux1(3)
 1130 format('| model type',t24,'|',a1,'titration',t37,'|',a1,'closed',
     $ t57,'|',a1,'open',t72,'|')
c
c-----------------------------------------------------------------------
c
c     Temperature parameters.
c     Note: ttk(1) = tk1, etc.
c
      write (newin,1090)
      ux1(1) = ' '
      ux1(2) = ' '
      if (jtemp .eq. 0) then
         ux1(1) = '*'
      elseif (jtemp .eq. 1) then
         ux1(2) = '*'
      else
         ux1(1) = '*'
      endif
      write (newin,1140) ux1(1),ux1(2)
 1140 format('| temperature model',t24,'|',a1,'power',t37,'|',a1,
     $ 'fluid mixing',t72,'|')
c
      write (newin,1090)
      write (newin,1150)
 1150 format('c power model  -->   temp = tstart + tk1*zi +',
     $ ' tk2*zi**2 + tk3*zi**3',t72,'|',
     $ /'c mixing model -->   temp = (tstart * tk1 + zi*tk2) /',
     $ ' (zi + tk1)',t72,'|',/'c',70('-'),'|')
      write (newin,1155) tempc0,(ttk(i), i = 1,3)
 1155 format('| tstart(c)|',1pg11.4,'|tk1|',1pg11.4,'|tk2|',1pg11.4,
     $ '|tk3|',1pg11.4,t72,'|')
c
c-----------------------------------------------------------------------
c
c      Zi and time parameters.
c
      write (newin,1090)
      write (newin,1160) zistrt,zimax
 1160 format('| starting value of zi |',1pg12.5,'|max. value of zi',
     $ t58,'|',1pg12.5,t72,'|')
c
      write (newin,1090)
      write (newin,1170) tstrt,timemx
 1170 format('| starting time (sec)  |',1pg12.5,'|max. time (sec)',
     $ t58,'|',1pg12.5,t72,'|')
c
      write (newin,1090)
      write (newin,1175) kstpmx,ksppmx
 1175 format('| max. steps',t24,'|',i10,2x,'|max. steps w/o print|',
     $ i10,t72,'|')
c
c     Note: cplim does not appear on the "D" format input file.
c
c     Print interval parameters.
c
      write (newin,1090)
      write (newin,1180) dzprnt,dzprlg
 1180 format('| linear print interval|',1pg12.5,
     $ '|log print interval  |',1pg12.5,t72,'|')
c
c     Plot interval parameters.
c
c     Note: dzplot, dzpllg, and ksplmx do not appear on the "D"
c     format input file.
c
c     Note: ifile does not appear on the "D" format input file.
c
c-----------------------------------------------------------------------
c
c     Nxopt options.
c
      write (newin,1090)
      write (newin,1190)
 1190 format('| suppress mineral phases',t72,'|')
      write (newin,1090)
c
      ux8(1) = ublk8
      ux8(2) = ublk8
      if (nxopt .gt. 0) then
        k = 0
        do 140 n = 1,nxopt
          k = k + 1
          if (uxopt(n)(1:8) .eq. 'all     ') then
            ux8(k) = uxopt(n)
          elseif (uxopt(n)(1:8) .eq. 'alwith  ') then
            ux8(k) = uxct16(n)
          else
            ux8(k) = uxct16(n)
          endif
          if (k .eq. 2) then
            write (newin,1200)  ux8(1),ux8(2)
 1200       format('| phases w/ elements|  ',a8,t46,'|  ',a8,t72,'|n')
            k = 0
            ux8(1) = ublk8
            ux8(2) = ublk8
          endif
  140   continue
        if (k .eq. 1) write (newin,1200)  ux8(1),ux8(2)
      else
        write (newin,1200)  ux8(1),ux8(2)
      endif
c
      ux24(1) = ublk24
      ux24(2) = ublk24
      if (nxopex .gt. 0) then
        k = 0
        do 150 n = 1,nxopex
          k = k + 1
          ux24(k) = uxopex(n)
          if (k .eq. 2) then
            write (newin,1210)  ux24(1),ux24(2)
 1210       format('| phases except',t21,'|  ',a22,t46,'|  ',a23,
     $      t72,'|n')
            k = 0
            ux24(1) = ublk24
            ux24(2) = ublk24
          endif
  150   continue
        if (k .eq. 1) write (newin,1210)  ux24(1),ux24(2)
      else
        write (newin,1210)  ux24(1),ux24(2)
      endif
c
c-----------------------------------------------------------------------
c
c     Nffg options.
c
      write (newin,1090)
      write (newin,1220)
 1220 format('| fixed fugacity phases- species, moles(per kg h2o), ',
     $ 'log fugacity(bars)',t72,'|')
      write (newin,1090)
      if (nffg .gt. 0) then
        write (newin,1230) (uffg(i),moffg(i),xlkffg(i), i = 1,nffg)
 1230 format('| ',a24,t38,'|',1pg12.5,t55,'|',1pg12.5,t72,'|n')
      else
        write (newin,1240)
 1240 format('| none',t38,'|',t55,'|',t72,'|n')
      endif
c
c-----------------------------------------------------------------------
c
c     Reactants and rate laws.
c
      write (newin,1090)
      write (newin,1250)
 1250 format('c',t20,'R A T E     L A W S',/,
     $ 'c 1 = relative',t33,'rate = rk1 + rk2*zi + (1/2)rk3*zi*zi',
     $ /'c 2 = transition state theory',t33,
     $ 'rate = CHECK DOCUMENTATION',
     $ /'c 3 = specified rate',
     $ /'c 4 = activity term rate',t33,'rate = CHECK DOCUMENTATION',
     $ /'c',/'c ',t20,'R E A C T A N T     T Y P E S',
     $ /'c mineral     solid solution     special    aqueous    gas',
     $ /'c',/'c ',t20,'S U R F A C E    T Y P E',
     $ /'c 0 = fixed surface area     1 = fixed specific surface area',
     $ /'c',/,'c',t20,'N O T E S',
     $ /'c status and jreac are normally not set by the user')
      write (newin,1090)
      write (newin,1260)
 1260 format('| reactants   (ss) solid solution only',
     $       5x,'(sp) special reactant only  |')
      write (newin,1090)
c
      if (nrct .eq. 0) then
        write (newin,1270)
 1270   format('| REACTANT',t20,'| none',t45,'|status',t55,'|',t72,'|')
      else
        nsrt = 0
        nxrt = 0
        do 390 nrc = 1,nrct
          if (nrc .gt. 1) write (newin,1280)
 1280     format('c',70('-'),'|')
c
c         Name, flags, and masses.
c
          write (newin,1290) ureac(nrc),jreac(nrc)
 1290     format('| REACTANT',t20,'|  ',a22,t45,'|status',t55,'|',i2,
     $    t72,'|')
          write (newin,1300) morr(nrc),modr(nrc)
 1300     format('| moles remaining',t20,'|',1pg12.5,t45,'|destroyed',
     $    t55,'|',1pg12.5,t72,'|')
          write (newin,1310) umint(jcode(nrc)),sk(nrc)
 1310     format('| reactant type',t20,'|',2x,a14,t45,'|sk',t55,'|',
     $    1pg12.5,t72,'|')
          write (newin,1320) nsk(nrc),fk(nrc)
 1320     format('| surface type',t20,'|',2x,i2,t45,'|fk',t55,
     $    '|',1pg12.5,t72,'|')
c
          if (jcode(nrc) .eq. 1) then
c
c           Solid solution compositions.
c
            nxrt = nxrt + 1
            do 170 iktb = 1,iktbt(nxrt)
              write (newin,1330) uendb(iktb,nxrt),rxbarb(iktb,nxrt)
 1330         format('| end-member',t20,'|',a24,t45,
     $        '|mole fr',t55,'| ',1pg14.6,t72,'|ss,n')
  170       continue
  180       if (iktbt(nxrt) .eq. 0) write (newin,1330) unone,zero
          else
            write (newin,1340)
 1340       format('| end-member',t20,'|',t45,'|mole fr',t55,'|',
     $       t72,'|ss,n')
          endif
c
         if (jcode(nrc) .eq. 2) then
c
c           Special reactant compositions.
c
            nsrt = nsrt + 1
            write (newin,1350) vreac(nrc)
 1350       format('| volume',t20,'|',1pg12.5,t45,'|',t55,'|',t72,'|sp')
c
            do 190 ncb = 1,nesrbt(nsrt)
              write (newin,1360) uesrb(ncb,nsrt),cesrb(ncb,nsrt)
 1360         format('| element',t20,'|  ',a8,t35,'|moles',t45,'|',
     $        1x,1pe25.15,t72,'|sp,n')
  190       continue
  200       if (nesrbt(nsrt) .eq. 0) write (newin,1360) unone(1:8),zero
          else
            write (newin,1370)
 1370       format('| volume',t20,'|',t45,'|',t55,'|',t72,'|sp',
     $      /'| element',t20,'|',t45,'|moles',t55,'|',t72,'|sp,n')
          endif
c
c         Rate law parameters, forward direction.
c
          write (newin,1380) nrk(1,nrc)
 1380     format('| DISSOLUTION LAW  |',i5,t45,'|',t55,'|',t72,'|')
c
          if (nrk(1,nrc) .eq. 1) then
            do 210 i = 1,imech(1,nrc)
              write (newin,1390) i,rk0(i,1,nrc),i
 1390         format('| rate constant rk',i1,t20,'|',1pg12.5,t45,
     $        '|csigma',i1,t55,'|',t72,'|')
  210       continue
          elseif (nrk(1,nrc) .eq. 2) then
            do 230 i = 1,imech(1,nrc)
              write (newin,1395) i,rk0(i,1,nrc),i,csigma(i,1,nrc)
 1395         format('| rate constant rk',i1,t20,'|',1pg12.5,t45,
     $        '|csigma',i1,t55,'| ',1pg14.6,t72,'|')
              do 220 n = 1,ndact(i,1,nrc)
                write (newin,1400) udac(n,i,1,nrc),cdac(n,i,1,nrc)
 1400           format('| aqueous species  | ',a23,t45,'|cdac',t55,
     $          '| ',1pg14.5,t72,'|n')
  220         continue
              write (newin,1410) trk0(i,1,nrc)
 1410         format('| temperature (c)  | ',1pg12.5,t45,'|',t55,'|',
     $        t72,'|234')
              if (iact(i,1,nrc) .eq. 1) then
                write (newin,1420) eact(i,1,nrc)
 1420           format('| act. energy-kcal | ',1pg12.5,t45,'|enthalpy',
     $          t55,'|',t72,'|234')
              elseif (iact(i,1,nrc) .eq. 2) then
                write (newin,1430) hact(i,1,nrc)
 1430           format('| act. energy-kcal | ',t45,'|enthalpy',
     $          t55,'| ',1pg12.5,t72,'|234')
              endif
  230       continue
          elseif (nrk(1,nrc) .eq. 3) then
            do 250 i = 1,imech(1,nrc)
              write (newin,1390) i,rk0(i,1,nrc),i
              do 240 n = 1,ndact(i,1,nrc)
                write (newin,1400) udac(n,i,1,nrc),cdac(n,i,1,nrc)
  240         continue
              write (newin,1410) trk0(i,1,nrc)
              if (iact(i,1,nrc) .eq. 1) then
                write (newin,1420) eact(i,1,nrc)
              elseif (iact(i,1,nrc) .eq. 2) then
                write (newin,1430) hact(i,1,nrc)
              endif
  250       continue
          elseif (nrk(1,nrc) .eq. 4) then
            do 270 i = 1,imech(1,nrc)
              write (newin,1390) i,rk0(i,1,nrc),i
              do 260 n = 1,ndact(i,1,nrc)
                write (newin,1400) udac(n,i,1,nrc),cdac(n,i,1,nrc)
  260         continue
              write (newin,1410) trk0(i,1,nrc)
              if (iact(i,1,nrc) .eq. 1) then
                write (newin,1420) eact(i,1,nrc)
              elseif (iact(i,1,nrc) .eq. 2) then
                write (newin,1430) hact(i,1,nrc)
              endif
  270       continue
          endif
c
c         Rate law parameters, backward direction.
c
          write (newin,1440) nrk(2,nrc)
 1440     format('| PRECIPITATION LAW|',i5,t45,'|',t55,'|',t72,'|')
c
          if (nrk(2,nrc) .eq. 1) then
            do 280 i = 1,imech(2,nrc)
              write (newin,1390) i,rk0(i,2,nrc),i
  280       continue
          elseif (nrk(2,nrc) .eq. 2) then
            do 300 i = 1,imech(2,nrc)
              write (newin,1395) i,rk0(i,2,nrc),i,csigma(i,2,nrc)
              do 290 n = 1,ndact(i,2,nrc)
                write (newin,1400) udac(n,i,2,nrc),cdac(n,i,2,nrc)
  290         continue
              write (newin,1410) trk0(i,2,nrc)
              if (iact(i,2,nrc) .eq. 1) then
                write (newin,1420) eact(i,2,nrc)
              elseif (iact(i,2,nrc) .eq. 2) then
                write (newin,1430) hact(i,2,nrc)
              endif
  300       continue
          elseif (nrk(2,nrc) .eq. 3) then
            do 320 i = 1,imech(2,nrc)
              write (newin,1390) i,rk0(i,2,nrc),i
              do 310 n = 1,ndact(i,2,nrc)
                write (newin,1400) udac(n,i,2,nrc),cdac(n,i,2,nrc)
  310         continue
              write (newin,1410) trk0(i,2,nrc)
              if (iact(i,2,nrc) .eq. 1) then
                write (newin,1420) eact(i,2,nrc)
              elseif (iact(i,2,nrc) .eq. 2) then
                write (newin,1430) hact(i,2,nrc)
              endif
  320       continue
          elseif (nrk(2,nrc) .eq. 4) then
            do 340 i = 1,imech(2,nrc)
              write (newin,1390) i,rk0(i,2,nrc),i
              do 330 n = 1,ndact(i,2,nrc)
                write (newin,1400) udac(n,i,2,nrc),cdac(n,i,2,nrc)
  330         continue
              write (newin,1410) trk0(i,2,nrc)
              if (iact(i,2,nrc) .eq. 1) then
                write (newin,1420) eact(i,2,nrc)
              elseif (iact(i,2,nrc) .eq. 2) then
                write (newin,1430) hact(i,2,nrc)
              endif
  340       continue
          endif
c
  390   continue
      endif
c
c-----------------------------------------------------------------------
c
c     Options. Iopt, iopr, and iodb option switches, but not iopg
c     option switches. Skip any of the above which happen to be
c     development options.
c
c     Note: iopt(1) = iopt1, etc.
c
      write (newin,1090)
      write (newin,1460)
 1460 format('| options',t72,'|')
      write (newin,1090)
c
      do 410 idescx = 1,nop6pa
        if (uvar6(idescx)(1:4) .ne. 'iopg') then
          j2 = ilnobl(uopt6(idescx))
          write (newin,1465) uopt6(idescx)(1:j2)
 1465     format('| - ',a,' -',t72,'|')
c
          do 400 idesc = 1,nod6pa
            if (iopti6(idesc) .eq. idescx) then
              if (uvar6(idescx)(1:4) .eq. 'iopt') then
                if (iopt(index6(idescx)) .eq. ivalu6(idesc)) then
                  ux1(1) = '*'
                else
                  ux1(1) = ' '
                endif
              elseif (uvar6(idescx)(1:4) .eq. 'iopr') then
                if (iopr(index6(idescx)) .eq. ivalu6(idesc)) then
                  ux1(1) = '*'
                else
                  ux1(1) = ' '
                endif
              elseif (uvar6(idescx)(1:4) .eq. 'iodb') then
                if (iodb(index6(idescx)) .eq. ivalu6(idesc)) then
                  ux1(1) = '*'
                else
                  ux1(1) = ' '
                endif
              else
                ux1(1) = ' '
              endif
              j2 = ilnobl(udesc6(idesc))
              write (newin,1470) ux1(1),udesc6(idesc)(1:j2)
 1470         format('|',t5,a1,t7,a,t72,'|')
            endif
  400     continue
        endif
  410 continue
c
c-----------------------------------------------------------------------
c
c     Development options.
c
      write (newin,1090)
      write (newin,1480)
 1480 format('| development options  (used for code development)',
     $ t72,'|')
      write (newin,1090)
c
      do 420 idescx = 1,ndv6pa
        ival = 0
        if (udv6vr(idescx)(1:4) .eq. 'iopt') then
          ival = iopt(idev6i(idescx))
        elseif (udv6vr(idescx)(1:4) .eq. 'iopr') then
          ival = iopr(idev6i(idescx))
        elseif (udv6vr(idescx)(1:4) .eq. 'iodb') then
          ival = iodb(idev6i(idescx))
        elseif (udv6vr(idescx)(1:4) .eq. 'iopg') then
          ival = iopg(idev6i(idescx))
        endif
        j2 = ilnobl(udevl6(idescx))
        write (newin,1490) ival,udevl6(idescx)(1:j2)
 1490   format('|',i5,1x,a,t72,'|')
  420 continue
c
c-----------------------------------------------------------------------
c
c     Tolerances.
c
      write (newin,1090)
      write (newin,1500)
 1500 format('| tolerances',t36,'desired values - defaults info-only |')
      write (newin,1090)
c
c     Maximum number of iterations.
c
      if (itermx .eq. 0) then
        write (newin,1510) utol6(1),udef(1)
 1510   format('| ',a32,'|',t54,'|',2x,a14,t72,'|')
      else
        write (newin,1520) utol6(1),itermx,udef(1)
 1520   format('| ',a32,'|',i10,t54,'|',2x,a14,t72,'|')
      endif
c
c     Dump interval.
c
      if (dlzidp .eq. 0) then
        write (newin,1510) utol6(2),udef(2)
      else
        write (newin,1530) utol6(2),dlzidp,udef(2)
 1530   format('| ',a32,'|',1pg12.5,t54,'|',2x,a14,t72,'|')
      endif
c
c     True tolerances.
c
      if (tolbt .eq. 0) then
        write (newin,1510) utol6(3),udef(3)
      else
        write (newin,1530) utol6(3),tolbt,udef(3)
      endif
      if (toldl .eq. 0) then
        write (newin,1510) utol6(4),udef(4)
      else
        write (newin,1530) utol6(4),toldl,udef(4)
      endif
      if (tolx .eq. 0) then
        write (newin,1510) utol6(5),udef(5)
      else
        write (newin,1530) utol6(5),tolx,udef(5)
      endif
      if (tolsat .eq. 0) then
        write (newin,1510) utol6(6),udef(6)
      else
        write (newin,1530) utol6(6),tolsat,udef(6)
      endif
      if (tolsst .eq. 0) then
        write (newin,1510) utol6(7),udef(7)
      else
        write (newin,1530) utol6(7),tolsst,udef(7)
      endif
c
c     Setscrew parameters.
c     Note: sscrew(1) = screw1, etc.
c
      if (sscrew(1) .eq. 0) then
        write (newin,1510) utol6(8),udef(8)
      else
        write (newin,1530) utol6(8),sscrew(1),udef(8)
      endif
      if (sscrew(2) .eq. 0) then
        write (newin,1510) utol6(9),udef(9)
      else
        write (newin,1530) utol6(9),sscrew(2),udef(9)
      endif
      if (sscrew(3) .eq. 0) then
        write (newin,1510) utol6(10),udef(10)
      else
        write (newin,1530) utol6(10),sscrew(3),udef(10)
      endif
      if (sscrew(4) .eq. 0) then
        write (newin,1510) utol6(11),udef(11)
      else
        write (newin,1530) utol6(11),sscrew(4),udef(11)
      endif
      if (sscrew(5) .eq. 0) then
        write (newin,1510) utol6(12),udef(12)
      else
        write (newin,1530) utol6(12),sscrew(5),udef(12)
      endif
      if (sscrew(6) .eq. 0) then
        write (newin,1510) utol6(13),udef(13)
      else
        write (newin,1530) utol6(13),sscrew(6),udef(13)
      endif
c
c     Z parameters.
c
      if (zklogu .eq. 0) then
        write (newin,1510) utol6(14),udef(14)
      else
        write (newin,1530) utol6(14),zklogu,udef(14)
      endif
      if (zklogl .eq. 0) then
        write (newin,1510) utol6(15),udef(15)
      else
        write (newin,1530) utol6(15),zklogl,udef(15)
      endif
      if (zkfac .eq. 0) then
        write (newin,1510) utol6(16),udef(16)
      else
        write (newin,1530) utol6(16),zkfac,udef(16)
      endif
c
c     Step size limits and maximum order.
c
      if (dlzmx1 .eq. 0) then
        write (newin,1510) utol6(17),udef(17)
      else
        write (newin,1530) utol6(17),dlzmx1,udef(17)
      endif
      if (dlzmx2 .eq. 0) then
        write (newin,1510) utol6(18),udef(18)
      else
        write (newin,1530) utol6(18),dlzmx2,udef(18)
      endif
      if (nordlm .eq. 0) then
        write (newin,1510) utol6(19),udef(19)
      else
        write (newin,1520) utol6(19),nordlm,udef(19)
      endif
c
c     Maximum number of phase assemblage tries.
c
      if (ntrymx .eq. 0) then
        write (newin,1510) utol6(20),udef(20)
      else
        write (newin,1520) utol6(20),ntrymx,udef(20)
      endif
c
c     Slide and scan control parameters.
c
      if (npslmx .eq. 0) then
        write (newin,1510) utol6(21),udef(21)
      else
        write (newin,1520) utol6(21),npslmx,udef(21)
      endif
      if (nsslmx .eq. 0) then
        write (newin,1510) utol6(22),udef(22)
      else
        write (newin,1520) utol6(22),nsslmx,udef(22)
      endif
      if (ioscan .eq. 0) then
        write (newin,1510) utol6(23),udef(23)
      else
        write (newin,1520) utol6(23),ioscan,udef(23)
      endif
      write (newin,1090)
c
c-----------------------------------------------------------------------
c
c     Bottom half of the current output file.
c
c     Write recaptured data in comments on the originating code, if
c     known.
c
      if (ucode(1:3).eq.'EQ6' .or. ucode(1:3).eq.'eq6') then
        write (newin,2180) urelno,ustage,ueqlrn,ueqlst
 2180   format('c pickup file written by eq6.',a4,a6,t72,'|',
     $  /'c  supported by eqlib.',a4,a6,t72,'|')
      endif
c
      if (ucode(1:5).eq.'EQ3NR' .or. ucode(1:5).eq.'eq3nr') then
        write (newin,2190) urelno,ustage,ueqlrn,ueqlst
 2190   format('c pickup file written by eq3nr.',a4,a6,t72,'|',
     $  /'c  supported by eqlib.',a4,a6,t72,'|')
      endif
c
c-----------------------------------------------------------------------
c
c     Title.
c
      do 500 n = 1,ntitl2
        write (newin,3010) utitl2(n)
 3010   format('|',a70,'|')
  500 continue
      write (newin,3020)
 3020 format('|',70('-'),'|')
c
c-----------------------------------------------------------------------
c
c     Original temperature.
c
      write (newin,3030) tempci
 3030 format('| temperature (C)',t38,'|',3x,1pg12.5,t72,'|')
      write (newin,3020)
c
c-----------------------------------------------------------------------
c
c     Electrical imbalance
c
      write (newin,3040) electr
 3040 format('| electrical imbalance',t38,'|',2x,1pe25.15,t72,'|')
      write (newin,3020)
c
c-----------------------------------------------------------------------
c
c     Number of aqueous basis species.
c
      write (newin,3050) ksq
 3050 format('| number of aqueous master species',t38,'|',2x,i5,t72,
     $ '|')
      write (newin,3020)
c
c     Index of the last pure mineral.
c
      write (newin,3060) kmt
 3060 format('| position of last pure mineral',t38,'|',2x,i5,t72,'|')
      write (newin,3020)
c
c     Index of the last solid solution.
c
      write (newin,3070) kxt
 3070 format('| position of last solid solution',t38,'|',2x,i5,t72,'|')
      write (newin,3020)
c
c-----------------------------------------------------------------------
c
c     Nxmod options.
c
      write (newin,3080)
 3080 format('| suppressed species',5x,
     $ '(suppress,replace,augmentk,augmentg)  value',t72,'|')
      write (newin,3020)
      if (nxmod .gt. 0) then
        do 520 n = 1,nxmod
          write (newin,3090) uxmd24(n),ujx(jxmod(n)),usup(kxmod(n)),
     $    xlkmod(n)
 3090     format('| ',a24,t26,'| ',a14,' | ',a8,' | ',1pg12.5,t72,'|')
  520   continue
      else
        write(newin,3095) unone
 3095     format('| ',a8,t26,'|',t43,'|',t54,'|',t72,'|')
      endif
      write (newin,3020)
c
c-----------------------------------------------------------------------
c
c     Iopg options.
c     Note: iopg(1) = iopg1, etc.
c
      write (newin,3100)
 3100 format('| iopg options',t72,'|')
c
      write (newin,3020)
      do 540 idescx = 1,nop6pa
        if (uvar6(idescx)(1:4) .eq. 'iopg') then
          j2 = ilnobl(uopt6(idescx))
          write (newin,3110) uopt6(idescx)(1:j2)
 3110     format('| - ',a,' -',t72,'|')
          do 530 idesc = 1,nod6pa
            if (iopti6(idesc) .eq. idescx) then
              if (iopg(index6(idescx)) .eq. ivalu6(idesc)) then
                ux1(1) = '*'
              else
                ux1(1) = ' '
              endif
              j2 = ilnobl(udesc6(idesc))
              write (newin,3120) ux1(1),udesc6(idesc)(1:j2)
 3120         format('|    ',a1,1x,a,t72,'|')
            endif
  530     continue
        endif
  540 continue
      write (newin,3020)
c
c-----------------------------------------------------------------------
c
c     Balance totals.
c
      write (newin,3130)
 3130 format('| elements, moles and moles aqueous',t72,'|')
      write (newin,3020)
c
      do 550 kc = 1,kct
        write (newin,3140) uelemb(kc),mteb(kc),mteaqb(kc)
 3140   format('|',10x,a8,t20,'|',1pe25.15,'|',1pe25.15,'|')
  550 continue
      write (newin,3020)
c
c-----------------------------------------------------------------------
c
c     Basis variable data.
c
      write (newin,3150)
 3150 format('| master species and logarithmic basis variables',t72,'|')
      write (newin,3020)
c
      do 560 kcol = 1,kdim
        write (newin,3160) unrms(kcol),undms(kcol),zvclgi(kcol)
 3160   format('|',a19,'|',a24,t46,'|',1pe25.15,t72,'|')
  560 continue
      write (newin,3020)
c
c-----------------------------------------------------------------------
c
c     Physically removed system data.
c
      write (newin,3170)
 3170 format('| physically removed subsystem ',
     $ ' (solid solution, mineral, moles)',t72,'|')
      write (newin,3020)
c
      if (kprs .gt. 0) then
        do 570 npr = 1,nprmn
          write (newin,3180) uprs(npr),mprs(npr)
 3180     format('|',t21,'|',a24,'|',1pe25.15,'|')
  570   continue
c
        uxx = '        '
        do 590 npr = nprmn + 1,nprmx
          if (uprs(npr)(1:24) .ne. uxx(1:24)) then
            write (newin,3190) uprs(npr)(1:24)
 3190       format('|',a24,'|',t46,'|',t72,'|')
            uxx = uprs(npr)(1:24)
          endif
          write (newin,3200) uprs(npr)(25:48),mprs(npr)
 3200     format('|',t21,'|',a24,'|',1pe25.15,'|')
  590   continue
      else
        write (newin,3210) unone
 3210     format('| ',a8,t21,'|',t46,'|',t72,'|')
      endif
      write (newin,3020)
c
c-----------------------------------------------------------------------
c
  999 continue
      end
