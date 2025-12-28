      subroutine wr6w72(cdac,cesrb,cplim,csigma,dlzmx1,dlzmx2,dlzidp,
     $ dzpllg,dzplot,dzprlg,dzprnt,eact,electr,fk,iact,ifile,iktbt,
     $ iktmax,imchmx,imech,iodb,iopg,iopr,iopt,ioscan,itermx,jcode,
     $ jreac,jtemp,jxmod,kct,kdim,kmax,kmt,kprs,ksq,ksplmx,ksppmx,
     $ kstpmx,kxmod,kxt,hact,modr,moffg,morr,mprs,mteaqb,mteb,nctmax,
     $ ndact,ndctmx,nesrbt,newin,nffg,nffgmx,nmodl1,nmodl2,nodbmx,
     $ nopgmx,noprmx,noptmx,nordlm,npslmx,nprmn,nprmx,nprsmx,nrct,
     $ nrctmx,nrk,nsk,nsrtmx,nsscmx,nsslmx,ntitl1,ntitl2,ntitmx,
     $ ntrymx,nttkmx,nxmdmx,nxmod,nxopex,nxopmx,nxopt,nxpemx,nxrtmx,
     $ rk0,rxbarb,sscrew,sk,tempci,tempc0,timemx,tolbt,toldl,tolsat,
     $ tolsst,tolx,trk0,tstrt,ttk,ucode,udac,uelemb,uendb,ueqlrn,
     $ ueqlst,uesrb,uffg,undms,unrms,uprs,ureac,urelno,ustage,utitl1,
     $ utitl2,uxct16,uxmd24,uxopex,uxopt,vreac,xlkffg,xlkmod,zimax,
     $ zistrt,zkfac,zklogl,zklogu,zvclgi)
c
c     This subroutine writes the EQ6 input file in compact ("W") format
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
     $ ksppmx,kstpmx,kxt,newin,nffg,nmodl1,nmodl2,nordlm,npslmx,nprmn,
     $ nprmx,nrct,nsslmx,ntitl1,ntitl2,ntrymx,nxmod,nxopex,nxopt
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
      real*8 cplim,dlzmx1,dlzmx2,dlzidp,dzpllg,dzplot,dzprlg,dzprnt,
     $ electr,tempci,tempc0,timemx,tolbt,toldl,tolsat,tolsst,tolx,
     $ tstrt,zkfac,zklogl,zklogu,zimax,zistrt
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,iktb,j,kc,kcol,n,ncb,npr,nrc,nsrt,nxrt
c
      character*24 uxx
      character*8 uendit
c
c-----------------------------------------------------------------------
c
      data uendit /'endit.  '/
c
c-----------------------------------------------------------------------
c
c     Write the new input file.
c
c     Title.
c
      do 110 n = 1,ntitl1
        write (newin,1100) utitl1(n)
 1100 format(a80)
  110 continue
      if (ntitl1 .lt. ntitmx) write (newin,1110)
 1110 format('endit.')
c
c-----------------------------------------------------------------------
c
c     Nmodl option switches.
c
      write (newin,1120) nmodl1,nmodl2
 1120 format(4x,'nmodl1= ',i2,14x,'nmodl2= ',i2)
c
c     Temperature parameters.
c     Note: ttk(1) = tk1, etc.
c
      write (newin,1140) tempc0,jtemp,(ttk(i), i = 1,3)
 1140 format(4x,'tempc0= ',1pe12.5,5x,'jtemp= ',i2,
     $/7x,'tk1= ',1pe12.5,7x,'tk2= ',1pe12.5,7x,'tk3= ',1pe12.5)
c
c      Zi and time parameters.
c
      write (newin,1160) zistrt,zimax,tstrt,timemx,kstpmx,cplim
 1160 format(4x,'zistrt= ',1pe12.5,5x,'zimax= ',1pe12.5,
     $ /5x,'tstrt= ',1pe12.5,4x,'timemx= ',1pe12.5,
     $ /4x,'kstpmx= ',i12,5x,'cplim= ',1pe12.5)
c
c     Print interval parameters.
c
      write (newin,1180) dzprnt,dzprlg,ksppmx
 1180 format(4x,'dzprnt= ',1pe12.5,4x,'dzprlg= ',1pe12.5,4x,
     $ 'ksppmx= ',i5)
c
c     Plot interval parameters.
c
      write (newin,1200) dzplot,dzpllg,ksplmx
 1200 format (4x,'dzplot= ',1pe12.5,4x,'dzpllg= ',1pe12.5,4x,
     $ 'ksplmx= ',i5)
c
c     Ifile.
c
      write (newin,1220) ifile
 1220 format(5x,'ifile= ',i2)
c
c-----------------------------------------------------------------------
c
c     Comment header for iopt, iopr, and iodb option switches.
c
      write (newin,1230)
 1230 format('*',15x,'1    2    3    4    5    6    7    8    9   10')
c
c     Iopt option switches.
c     Note: iopt(1) = iopt1, etc.
c
      write (newin,1250) (iopt(i), i = 1,20)
 1250 format(2x,'iopt1-10= ',10i5,/1x,'iopt11-20= ',10i5)
c
c     Iopr option switches.
c     Note: iopr(1) = iopY1, etc.
c
      write (newin,1270) (iopr(i), i = 1,20)
 1270 format(2x,'iopr1-10= ',10i5,/1x,'iopr11-20= ',10i5)
c
c     Iodb option switches.
c     Note: iodb(1) = iodb1, etc.
c
      write (newin,1290) (iodb(i), i = 1,20)
 1290 format(2x,'iodb1-10= ',10i5,/1x,'iodb11-20= ',10i5)
c
c-----------------------------------------------------------------------
c
c     Number of nxopt options.
c
      write (newin,1310) nxopt
 1310 format(5x,'nxopt= ',i2)
c
c     Nxopt options.
c
      if (nxopt .gt. 0) then
        do 140 n = 1,nxopt
          write (newin,1330) uxopt(n),uxct16(n)
 1330     format(4x,'option= ',a6,1x,a8)
  140   continue
        write (newin,1350) nxopex
 1350   format(4x,'nxopex= ',i2)
        if (nxopex .gt. 0) then
          do 150 n = 1,nxopex
            write (newin,1370) uxopex(n)
 1370       format(1x,'exception= ',a8)
  150     continue
        endif
      endif
c
c-----------------------------------------------------------------------
c
c     Number of nffg options.
c
      write (newin,1390) nffg
 1390 format(6x,'nffg= ',i2)
c
c     Nffg options.
c
      if (nffg .gt. 0) then
      do 160 n = 1,nffg
        write (newin,1410) uffg(n),moffg(n),xlkffg(n)
 1410   format(3x,'species= ',a24,/5x,'moffg= ',1pe12.5,4x,
     $  'xlkffg= ',1pe12.5)
  160   continue
      endif
c
c-----------------------------------------------------------------------
c
c     Number of reactants.
c
      write (newin,1430) nrct
 1430 format(6x,'nrct= ',i2)
      write (newin,1435)
 1435 format('*-------------------------------------------------',
     $ '----------------------------')
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
          write (newin,1450) ureac(nrc),jcode(nrc),jreac(nrc),
     $    morr(nrc),modr(nrc)
 1450     format(2x,'reactant= ',a24,/5x,'jcode= ',i2,15x,'jreac= ',i2,
     $    /6x,'morr= ',1pe12.5,6x,'modr= ',1pe12.5)
c
          if (jcode(nrc) .eq. 1) then
c
c           Solid solution compositions.
c
            nxrt = nxrt + 1
            do 170 iktb = 1,iktbt(nxrt)
              write (newin,1470) uendb(iktb,nxrt),rxbarb(iktb,nxrt)
 1470         format(3x,a24,1x,1pe12.5)
  170       continue
            write (newin,1475) uendit
 1475       format(3x,a8)
          elseif (jcode(nrc) .eq. 2) then
c
c           Special reactant compositions.
c
            nsrt = nsrt + 1
            write (newin,1490) vreac(nrc)
 1490       format(5x,'vreac= ',1pe12.5)
            do 190 ncb = 1,nesrbt(nsrt)
              write (newin,1510) uesrb(ncb,nsrt),cesrb(ncb,nsrt)
 1510         format(3x,a8,1x,1pe22.15)
  190       continue
            write (newin,1475) uendit
          endif
c
c         Surface area parameters.
c
          write (newin,1530) nsk(nrc),sk(nrc),fk(nrc)
 1530     format(7x,'nsk= ',i2,18x,'sk= ',1pe12.5,8x,'fk= ',1pe12.5)
c
          write (newin,1550) nrk(1,nrc),nrk(2,nrc)
 1550     format(7x,'nrk= ',i2,16x,'nrpk= ',i2)
c
c         Rate law parameters, forward direction.
c
          if (nrk(1,nrc) .eq. 1) then
            i = imech(1,nrc)
            write (newin,1570) (rk0(j,1,nrc), j = 1,i)
 1570       format(7x,'rk1= ',1pe12.5,7x,'rk2= ',1pe12.5,7x,'rk3= ',
     $      1pe12.5)
          elseif (nrk(1,nrc) .eq. 2) then
            write (newin,1590) imech(1,nrc)
 1590       format(5x,'imech= ',i2)
            do 220 i = 1,imech(1,nrc)
              write (newin,1610) rk0(i,1,nrc),trk0(i,1,nrc),
     $        iact(i,1,nrc)
 1610         format(7x,'rk0= ',1pe12.5,6x,'trk0= ',1pe12.5,6x,
     $        'iact= ',i2)
              write (newin,1640) eact(i,1,nrc),hact(i,1,nrc)
 1640         format(6x,'eact= ',1pe12.5,6x,'hact= ',1pe12.5)
              write (newin,1660) ndact(i,1,nrc),csigma(i,1,nrc)
 1660         format(5x,'ndact= ',i2,14x,'csigma= ',1pe12.5)
              if (ndact(i,1,nrc) .gt. 0) then
                do 210 n = 1,ndact(i,1,nrc)
                  write (newin,1680) udac(n,i,1,nrc),cdac(n,i,1,nrc)
 1680             format(6x,'udac= ',a24,6x,'cdac= ',1pe12.5)
  210           continue
              endif
  220       continue
          elseif (nrk(1,nrc) .eq. 3) then
            i = imech(1,nrc)
            write (newin,1610) rk0(i,1,nrc),trk0(i,1,nrc),
     $      iact(i,1,nrc)
            write (newin,1640) eact(i,1,nrc),hact(i,1,nrc)
          elseif (nrk(1,nrc) .eq. 4) then
            write (newin,1590) imech(1,nrc)
            do 240 i = 1,imech(1,nrc)
              write (newin,1610) rk0(i,1,nrc),trk0(i,1,nrc),
     $        iact(i,1,nrc)
              write (newin,1640) eact(i,1,nrc),hact(i,1,nrc)
              write (newin,1800) ndact(i,1,nrc)
 1800         format(5x,'ndact= ',i2)
              if (ndact(i,1,nrc) .gt. 0) then
                do 230 n = 1,ndact(i,1,nrc)
                  write (newin,1680) udac(n,i,1,nrc),cdac(n,i,1,nrc)
  230           continue
              endif
  240       continue
          endif
c
c         Rate law parameters, backward direction.
c
          if (nrk(2,nrc) .eq. 1) then
            i = imech(2,nrc)
            write (newin,1830) (rk0(j,2,nrc), j = 1,i)
 1830       format(6x,'rpk1= ',1pe12.5,6x,'rpk2= ',1pe12.5,6x,'rpk3= ',
     $      1pe12.5)
          elseif (nrk(2,nrc) .eq. 2) then
            write (newin,1840) imech(2,nrc)
 1840       format(4x,'ipmech= ',i2)
            do 260 i = 1,imech(2,nrc)
              write (newin,1850) rk0(i,2,nrc),trk0(i,2,nrc),
     $        iact(i,2,nrc)
 1850         format(6x,'rpk0= ',1pe12.5,5x,'trpk0= ',1pe12.5,5x,
     $        'iactp= ',i2)
              write (newin,1860) eact(i,2,nrc),hact(i,2,nrc)
 1860         format(5x,'eactp= ',1pe12.5,5x,'hactp= ',1pe12.5)
              write (newin,1870) ndact(i,2,nrc),csigma(i,2,nrc)
 1870         format(4x,'npdact= ',i2,14x,'cpsigm= ',1pe12.5)
c
              if (ndact(i,2,nrc) .gt. 0) then
                do 250 n = 1,ndact(i,2,nrc)
                  write (newin,1880) udac(n,i,2,nrc),cdac(n,i,2,nrc)
 1880             format(5x,'updac= ',a24,5x,'cpdac= ',1pe12.5)
  250           continue
              endif
  260       continue
          elseif (nrk(2,nrc) .eq. 3) then
            i = imech(2,nrc)
            write (newin,1850) rk0(i,2,nrc),trk0(i,2,nrc),iact(i,2,nrc)
            write (newin,1860) eact(i,2,nrc),hact(i,2,nrc)
          elseif (nrk(2,nrc) .eq. 4) then
            write (newin,1840) imech(2,nrc)
            do 280 i = 1,imech(2,nrc)
              write (newin,1850) rk0(i,2,nrc),trk0(i,2,nrc),
     $        iact(i,2,nrc)
              write (newin,1860) eact(i,2,nrc),hact(i,2,nrc)
              write (newin,1940) ndact(i,2,nrc)
 1940         format(4x,'npdact= ',i2)
              if (ndact(i,2,nrc) .gt. 0) then
                do 270 n = 1,ndact(i,2,nrc)
                  write (newin,1880) udac(n,i,2,nrc),cdac(n,i,2,nrc)
  270           continue
              endif
  280       continue
          endif
c
          write (newin,1435)
  390   continue
      endif
c
c-----------------------------------------------------------------------
c
c     Dump interval.
c
      write (newin,2010) dlzidp
 2010 format(4x,'dlzidp= ',1pe12.5)
c
c     Tolerances.
c
      write (newin,2030) tolbt,toldl,tolx,tolsat,tolsst
 2030 format(5x,'tolbt= ',1pe12.5,5x,'toldl= ',1pe12.5,6x,'tolx= ',
     $ 1pe12.5,/4x,'tolsat= ',1pe12.5,4x,'tolsst= ',1pe12.5)
c
c     Setscrew parameters.
c     Note: sscrew(1) = screw1, etc.
c
      write (newin,2050) (sscrew(i), i = 1,6)
 2050 format(4x,'screw1= ',1pe12.5,4x,'screw2= ',1pe12.5,4x,'screw3= ',
     $ 1pe12.5,/4x,'screw4= ',1pe12.5,4x,'screw5= ',1pe12.5,
     $ 4x,'screw6= ',1pe12.5)
c
c     Z parameters.
c
      write (newin,2070) zklogu,zklogl,zkfac
 2070 format(4x,'zklogu= ',f12.3,4x,'zklogl= ',f12.3,5x,'zkfac= ',f12.3)
c
c     Step size limits and maximum order.
c
      write (newin,2090) dlzmx1,dlzmx2,nordlm
 2090 format(4x,'dlzmx1= ',1pe12.5,4x,'dlzmx2= ',1pe12.5,
     $ 4x,'nordlm= ',i2)
c
c     Maximum number of iterations and maximum number of phase
c     assemblage tries.
c
      write (newin,2110) itermx,ntrymx
 2110 format(4x,'itermx= ',i2,14x,'ntrymx= ',i2)
c
c     Slide and scan control parameters.
c
      write (newin,2130) npslmx,nsslmx,ioscan
 2130 format(4x,'npslmx= ',i2,14x,'nsslmx= ',i2,14x,'ioscan= ',i2)
c
      write (newin,1435)
c
c-----------------------------------------------------------------------
c
c     Process the bottom half of the current output file.
c
c     Write recaptured data in comments on the originating code, if
c     known.
c
      if (ucode(1:3).eq.'EQ6' .or. ucode(1:3).eq.'eq6') then
        write (newin,2180) urelno,ustage,ueqlrn,ueqlst
 2180   format('* pickup file written by EQ6, version ',a4,' (',a8,')',
     $  /'*  supported by EQLIB, version ',a4,' (',a8,')')
      endif
c
      if (ucode(1:5).eq.'EQ3NR' .or. ucode(1:5).eq.'eq3nr') then
        write (newin,2190) urelno,ustage,ueqlrn,ueqlst
 2190   format('* pickup file written by EQ3NR, version ',a4,' (',a8,
     $  ')',/'*  supported by EQLIB, version ',a4,' (',a8,')')
      endif
c
c-----------------------------------------------------------------------
c
c     Title.
c
      do 400 n = 1,ntitl2
        write (newin,3000) utitl2(n)
 3000   format(a80)
  400 continue
      if (ntitl2 .lt. ntitmx) write (newin,1110)
c
c-----------------------------------------------------------------------
c
c     Original temperature.
c
      write (newin,3020) tempci
 3020 format(4x,'tempci= ',1pe12.5)
c
c-----------------------------------------------------------------------
c
c     Number of nxmod options.
c
      write (newin,3040) nxmod
 3040 format(5x,'nxmod= ',i2)
c
c-----------------------------------------------------------------------
c
c     Nxmod options.
c
      if (nxmod .gt. 0) then
        do 420 n = 1,nxmod
          write (newin,3060) uxmd24(n),jxmod(n),kxmod(n),xlkmod(n)
 3060     format(3x,'species= ',a24,/6x,'type= ',i2,14x,'option= ',i2,
     $    14x,'xlkmod= ',1pe12.5)
  420   continue
      endif
c
c-----------------------------------------------------------------------
c
c     Iopg options.
c     Note: iopg(1) = iopg1, etc.
c
      write(newin,3080) (iopg(i), i = 1,10)
 3080 format(5x,'iopg1= ',i2,15x,'iopg2= ',i2,15x,'iopg3= ',i2,
     $       /5x,'iopg4= ',i2,15x,'iopg5= ',i2,15x,'iopg6= ',i2,
     $       /5x,'iopg7= ',i2,15x,'iopg8= ',i2,15x,'iopg9= ',i2,
     $       /4x,'iopg10= ',i2)
c
c-----------------------------------------------------------------------
c
c     Index limits.
c
      write (newin,3100) kct,ksq,kmt,kxt,kdim,kprs
 3100 format(7x,'kct= ',i2,17x,'ksq= ',i2,17x,'kmt= ',i2,
     $ /7x,'kxt= ',i2,16x,'kdim= ',i2,16x,'kprs= ',i2)
c
c-----------------------------------------------------------------------
c
c     Balance totals.
c
      write (newin,3110)
 3110 format('*  Component',4x,'Moles Total',12x,'Moles Aqueous')
      do 430 kc = 1,kct
        write (newin,3130) uelemb(kc),mteb(kc),mteaqb(kc)
 3130   format(3x,a8,1x,1pe22.15,1x,1pe22.15)
  430 continue
      write (newin,3150) electr
 3150 format(3x,'electr  ',1x,1pe22.15)
c
c-----------------------------------------------------------------------
c
c     Basis variable data.
c
      do 440 kcol = 1,kdim
        write (newin,3170) unrms(kcol),undms(kcol),zvclgi(kcol)
 3170   format(3x,a24,1x,a24,1x,1pe22.15)
  440 continue
c
c-----------------------------------------------------------------------
c
c     Physically removed system data.
c
      if (kprs .gt. 0) then
c
        do 450 npr = 1,nprmn
          write (newin,3190) uprs(npr),mprs(npr)
 3190     format(3x,a24,1x,1pe22.15)
  450   continue
        write (newin,3200) uendit
 3200   format(3x,a8)
c
        uxx  = '        '
        do 480 npr = nprmn + 1,nprmx
          if (uprs(npr)(1:8) .eq. '        ') go to 490
          if (uprs(npr)(1:24) .eq. uxx(1:24)) then
            write (newin,3200) uendit
            write (newin,3210) uprs(npr)(1:24)
 3210       format(3x,a24)
          endif
          write (newin,3230) uprs(npr)(25:48),mprs(npr)
 3230     format(3x,a24,1x,1pe22.15)
  480   continue
  490   continue
        write (newin,3200) uendit
        write (newin,3200) uendit
      endif
c
  999 continue
      end
