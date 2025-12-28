      subroutine wr3d72(cspb,fep,iktmax,iodb,iopg,iopr,iopt,itermx,
     $ jflagb,jxmod,kxmod,ncompb,newin,nodbmx,nopgmx,noprmx,noptmx,
     $ nsq,nsqmax,ntitl,ntitmx,nttyo,nxmdmx,nxmod,nxtb,nxtmax,
     $ rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat,ubasis,uebal,
     $ umemb,uphas1,uphas2,uredox,usolb,uspecb,utitl,xbarb,
     $ uxmd24,xlkmod)
c
c     This subroutine writes the EQ3NR input file in menu-style ("D")
c     format for version 7.0-7.1.
c
c     This subroutine is called by:
c
c       XCON3/xcon3.f
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer iktmax,nodbmx,nopgmx,noprmx,noptmx,nsqmax,ntitmx,nxmdmx,
     $ nxtmax
c
      integer iodb(nodbmx),iopg(nopgmx),iopr(noprmx),iopt(noptmx),
     $ jflagb(nsqmax),jxmod(nxmdmx),kxmod(nxmdmx),ncompb(nxtmax)
c
      integer itermx,newin,nsq,ntitl,nttyo,nxmod,nxtb
c
      character*80 utitl(ntitmx)
      character*24 ubasis(nsqmax),umemb(iktmax,nxtmax),uphas1(nsqmax),
     $ uphas2(nsqmax),usolb(nxtmax),uspecb(nsqmax),uxmd24(nxmdmx)
      character*24 uebal,uredox
c
      real*8 cspb(nsqmax),xbarb(iktmax,nxtmax),xlkmod(nxmdmx)
c
      real*8 fep,rho,tempc,tdspkg,tdspl,tolbt,toldl,tolsat
c
c-----------------------------------------------------------------------
c
      include 'xcon3/x3op7.h'
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,idesc,idescx,ik,iktb,j,jfl,jj,j2,j3,n,ncount,nredox,
     $ ns,nx
      integer ilnobl
c
      logical qlcase
c
      character*24 unone,ux24
      character*16 udef(nto3pa),ujx(0:3)
      character*16 ux16
      character*8 usup(-1:2)
      character*1 ux1(3)
c
      real*8 xx
c
c-----------------------------------------------------------------------
c
      data unone  /'none                    '/
c
      data udef(1)  / '1.0e-10  tolbt  '/
      data udef(2)  / '1.0e-10  toldl  '/
      data udef(3)  / '0.5      tolsat '/
      data udef(4)  / '30       itermx '/
c
      data ujx(0)   / 'aqueous         '/
      data ujx(1)   / 'mineral         '/
      data ujx(2)   / 'gas             '/
      data ujx(3)   / 'solid solution  '/
c
      data usup(-1) / 'suppress'/
      data usup(0)  / 'replace '/
      data usup(1)  / 'augmentk'/
      data usup(2)  / 'augmentg'/
c
c-----------------------------------------------------------------------
c
c     Write the new input file.
c
c     Title.
c
      write (newin,1090)
 1090 format('|',70('-'),'|')
      do 110 n = 1,ntitl
        write (newin,1100) utitl(n)
 1100   format('|',a70,'|')
  110 continue
      write (newin,1090)
c
c-----------------------------------------------------------------------
c
c     Temperature and density.
c
      write (newin,1050) tempc,rho
 1050 format('|Temperature (C)         |',f6.2,8x,'|Density(gm/cm3)|',
     $ f9.5,t72,'|')
      write (newin,1090)
c
c-----------------------------------------------------------------------
c
c     Total dissolved salts.
c
      ux1(1) = ' '
      ux1(2) = ' '
      ux1(3) = ' '
      if (tdspkg .gt. 0.) then
        ux1(1) = '*'
        xx = tdspkg
      elseif (tdspl .gt. 0.) then
        ux1(2) = '*'
        xx = tdspl
      else
        ux1(3) = '*'
        xx = 0.
      endif
c
      if (xx .gt. 0.) then
        write (newin,1060) xx,(ux1(i), i = 1,3)
 1060   format('|Total Dissolved Salts   | ',g12.5,t41,'|',a1,'mg/kg',
     $  t49,'|',a1,'mg/l',t57,'|',a1,'not used',t72,'|')
      else
        write (newin,1070) (ux1(i), i = 1,3)
 1070   format('|Total Dissolved Salts   | ',t41,'|',a1,'mg/kg',
     $  t49,'|',a1,'mg/l',t57,'|',a1,'not used',t72,'|')
      endif
      write (newin,1095)
 1095 format('|',78('-'),'|')
c
c-----------------------------------------------------------------------
c
c     Species for electrical balancing.
c
      ux1(1) = ' '
      ux1(2) = ' '
      if (uebal(1:8) .eq. 'pick1.  ') then
        ux1(1) = '*'
        uebal = ' '
      elseif (uebal(1:8) .eq. '        ') then
        ux1(2) = '*'
      endif
c
      write (newin,1080) uebal,ux1(1),ux1(2)
 1080 format('|Electrical Balancing on |',a24,'|',a1,'code selects|',
     $ a1,'not performed|')
      write (newin,1095)
c
c-----------------------------------------------------------------------
c
c     Basis species and associated constraints.
c
      write (newin,1110)
 1110 format('| SPECIES',16x,'| BASIS SWITCH/CONSTRAINT|',
     $ ' CONC/ETC  | UNITS OR TYPE  |')
      write (newin,1095)
c
      if (nsq .le. 0) then
        write (newin,1120)
 1120   format('| none',t26,'|',t51,'|',t63,'|',t80,'|')
        write (newin,1095)
      endif
c
      ux16 = ' '
      if (iopt(1) .eq. -2) then
        ux16 = 'pe'
      elseif (iopt(1) .eq. -1) then
        ux16 = 'Eh'
      elseif (iopt(1) .eq. 0) then
        ux16 = 'LogfO2'
      endif
c
      if (ux16(1:8) .ne. '        ') then
        write (newin,1130) fep,ux16
 1130   format('|redox',t26,'|',t51,'|',g11.5,'|',a16,'|')
      endif
c
      if (iopt(1) .eq. -3) then
        qlcase = .true.
        do 115 ns = 1, nsq
          ux24 = uspecb(ns)
          call locase(ux24)
          if (ux24(1:24). ne. uspecb(ns)(1:24)) then
            qlcase = .false.
            go to 117
          endif
  115   continue
  117   if (qlcase) then
          uredox = 'o2(g)'
        else
          uredox = 'O2(g)'
        endif
      endif
c
      nredox = 0
      if (iopt(1).eq.1 .or. iopt(1).eq.-3) then
        write (newin,1140)
1140    format('|redox',t26,'| Defined by next species',t51,'|',t63,
     $  '|Redox couple',t80,'|')
        do 120 ns = 1,nsq
          if (uredox(1:24) .eq. uspecb(ns)(1:24)) then
            nredox = ns
            go to 130
          endif
  120   continue
c
        j2 = ilnobl(uredox)
        write (nttyo,1145) uredox(1:j2)
        write (newin,1145) uredox(1:j2)
 1145   format(" * Error - (XCON3/wr3d72) Can't find the redox",
     $  /7x,'couple-defining  species "',a,'" in the set of',
     $  /7x,'basis species.')
        go to 999
c
  130   if (nredox .ne. 1) then
          ux24 = uspecb(1)
          uspecb(1) = uspecb(nredox)
          uspecb(nredox) = ux24
c
          xx = cspb(1)
          cspb(1) = cspb(nredox)
          cspb(nredox) = xx
c
          jfl = jflagb(1)
          jflagb(1) = jflagb(nredox)
          jflagb(nredox) = jfl
c
          ux24 = ubasis(1)
          ubasis(1) = ubasis(nredox)
          ubasis(nredox) = ux24
c
          ux24 = uphas1(1)
          uphas1(1) = uphas1(nredox)
          uphas1(nredox) = ux24
c
          ux24 = uphas2(1)
          uphas2(1) = uphas2(nredox)
          uphas2(nredox) = ux24
c
          nredox = 1
        endif
      endif
c
      do 180 ns = 1,nsq
        if (uphas1(ns)(1:8) .ne. '        ') then
          ux24 = uphas1(ns)
        else
          ux24 = ubasis(ns)
        endif
        jfl = jflagb(ns)
        write (newin,1150) uspecb(ns),ux24,cspb(ns),ujflg7(jfl)
 1150   format('|',a24,'|',a24,'|',g11.5,'|',a16,'|')
c
        if (jfl .eq. 20) then
          do 140 nx = 1,nxtb
            if (uphas1(ns)(1:24) .eq. usolb(nx)(1:24)) go to 150
  140     continue
c
          j2 = ilnobl(uphas1(ns))
          write (nttyo,1160) uphas1(ns)(1:j2)
          write (newin,1160) uphas1(ns)(1:j2)
 1160     format(" * Error - (XCON3/wr3d72) Can't find the solid",
     $    /7x,'solution "',a,'" in the set of solid solutions',
     $    /7x,'for which compositions are defined.')
          go to 999
c
  150     continue
          iktb = ncompb(nx)
          do 160 ik = 1,iktb
            if (uphas2(ns)(1:24) .eq. umemb(ik,nx)(1:24)) go to 170
  160     continue
c
          j3 = ilnobl(uphas2(ns))
          write (nttyo,1170) uphas2(ns)(1:j3),uphas1(ns)(1:j2)
          write (newin,1170) uphas2(ns)(1:j3),uphas1(ns)(1:j2)
 1170     format(" * Error - (XCON3/wr3d72) Can't find the end-member",
     $    /7x,'"',a,'" in the set of components used to,'
     $    /7x,'specify the composition of solid solution "',a,'".')
          go to 999
c
  170     write (newin,1180) umemb(ik,nx),xbarb(ik,nx)
 1180     format('|',24x,'|',a24,'|',g11.5,'|Mole fraction',t80,'|')
c
          do 175 i = 1,iktb
            if (i .ne. ik) write (newin,1180) umemb(i,nx),xbarb(i,nx)
  175     continue
        endif
  180 continue
c
  190 write (newin,1095)
c
c-----------------------------------------------------------------------
c
c     Input solid solutions, excluding those defined as part of
c     solubility equilibrium constraints.
c
      write (newin,1800)
 1800 format('|Input Solid Solutions',t80,'|')
      write (newin,1095)
c
      ncount = 0
      if (nxtb .gt. 0) then
        do 210 nx = 1,nxtb
          do 200 ns = 1,nsq
            if (jflagb(ns) .eq. 20) then
              if (usolb(nx)(1:24) .ne. uphas1(ns)(1:24)) then
                ncount = ncount + 1
                write (newin,1810) usolb(nx)
 1810           format('|',a24,'|',t51,'|',t63,'|',t80,'|')
                iktb = ncompb(nx)
                do 195 i = 1,iktb
                  write (newin,1815) umemb(i,nx),xbarb(i,nx)
 1815             format('|',t26,'|',a24,'|',g11.5,'|Mole fraction',
     $            t80,'|')
  195           continue
              endif
            endif
  200     continue
  210   continue
      endif
c
      if (ncount .le. 0) then
        write (newin,1820)
 1820   format('| none',t26,'|',t51,'|',t63,'|',t80,'|')
      endif
c
      write (newin,1095)
c
c-----------------------------------------------------------------------
c
c     Nxmod options.
c
      write (newin,3080)
 3080 format('|SUPPRESSED SPECIES',3x,
     $ '(suppress,replace,augmentk,augmentg)    value',t72,'|')
      write (newin,1090)
      if (nxmod .gt. 0) then
        do 520 n = 1,nxmod
          write (newin,3090) uxmd24(n),ujx(jxmod(n)),usup(kxmod(n)),
     $    xlkmod(n)
 3090     format('| ',a24,t26,'| ',a14,' | ',a8,' | ',g12.5,t72,'|')
  520   continue
      else
        write(newin,3095) unone
 3095     format('| ',a8,t26,'|',t43,'|',t54,'|',t72,'|')
      endif
      write (newin,1090)
c
c-----------------------------------------------------------------------
c
c     Options. Iopt, iopg, and iopr option switches, but not iodb
c     option switches. Skip any of the above which happen to be
c     development options.
c
c     Note: iopt(1) = iopt1, etc.
c
      write (newin,1460)
 1460 format('|OPTIONS',t72,'|')
      write (newin,1090)
c
      do 410 idescx = 1,nop3pa
        if (uvar3(idescx)(1:4) .ne. 'iodb') then
          j = ilnobl(uopt3(idescx))
          write (newin,1465) uopt3(idescx)(1:j)
 1465     format('| - ',a,' -',t72,'|')
c
          do 400 idesc = 1,nod3pa
            if (iopti3(idesc) .eq. idescx) then
              if (uvar3(idescx)(1:4) .eq. 'iopt') then
                if (iopt(index3(idescx)) .eq. ivalu3(idesc)) then
                  ux1(1) = '*'
                else
                  ux1(1) = ' '
                endif
              elseif (uvar3(idescx)(1:4) .eq. 'iopg') then
                if (iopg(index3(idescx)) .eq. ivalu3(idesc)) then
                  ux1(1) = '*'
                else
                  ux1(1) = ' '
                endif
              elseif (uvar3(idescx)(1:4) .eq. 'iopr') then
                if (iopr(index3(idescx)) .eq. ivalu3(idesc)) then
                  ux1(1) = '*'
                else
                  ux1(1) = ' '
                endif
              else
                ux1(1) = ' '
              endif
              j = ilnobl(udesc3(idesc))
              write (newin,1470) ux1(1),udesc3(idesc)(1:j)
 1470         format('|',t5,a1,t7,a,t72,'|')
            endif
  400     continue
        endif
  410 continue
      write (newin,1090)
c
c-----------------------------------------------------------------------
c
c     Iodb switches.
c
      write (newin,1480)
 1480 format('|DEBUGGING SWITCHES (0 = off, 1,2 = on',t72,'|')
      write (newin,1090)
c
      do 420 idescx = 1,ndb3pa
        i = idbugi(idescx)
        j2 = ilnobl(udebug(idescx))
        j = index(udebug(idescx),'generic debugging')
        if (j .eq. 0) j = index(udebug(idescx),'pre-Newton-Raphson')
        if (j .eq. 0) then
          j = index(udebug(idescx),'stoichiometric factors')
          if (j .gt. 0) then
            jj = index(udebug(idescx),'factors calculation')
            if (jj .gt. 0) j = 0
          endif
        endif
        if (j .eq. 0) then
          write (newin,1490) iodb(i),udebug(idescx)(1:j2)
 1490     format('|',i1,2x,a,t72,'|')
        else
          write (newin,1495) iodb(i),udebug(idescx)(1:j2)
 1495     format('|',i1,2x,a,t72,'|2')
        endif
  420 continue
      write (newin,1090)
c
c-----------------------------------------------------------------------
c
c     Development options (there are none).
c
      write (newin,1700)
 1700 format('|DEVELOPMENT OPTIONS  (used for code development)',
     $ t72,'|')
      write (newin,1090)
      write (newin,1710)
 1710 format('| none',t72,'|')
      write (newin,1090)
c
c-----------------------------------------------------------------------
c
c     Tolerances.
c
      write (newin,1500)
 1500 format('|TOLERANCES',19x,'(desired values)',8x,'(defaults)',
     $ 7x,'|')
      write (newin,1090)
c
      if (tolbt .eq. 0) then
        write (newin,1510) utol3(1),udef(1)
 1510   format('|',a27,t28,'|',t51,'|',a16,t72,'|')
      else
        write (newin,1520) utol3(1),tolbt,udef(1)
 1520   format('|',a27,t28,'|',i10,t51,'|',a16,t72,'|')
      endif
c
      if (toldl .eq. 0) then
        write (newin,1510) utol3(2),udef(2)
      else
        write (newin,1520) utol3(2),toldl,udef(2)
      endif
c
      if (tolsat .eq. 0) then
        write (newin,1510) utol3(3),udef(3)
      else
        write (newin,1520) utol3(3),tolsat,udef(3)
      endif
c
      if (itermx .eq. 0) then
        write (newin,1510) utol3(4),udef(4)
      else
        write (newin,1520) utol3(4),itermx,udef(4)
      endif
c
      write (newin,1090)
c
c-----------------------------------------------------------------------
c
  999 continue
      end
