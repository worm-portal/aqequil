      subroutine wrtabc(acflg,actlg,actw,afrc1,aft1,alk,conclg,
     $ cteaq,ctb,dvoso,dwoso,eh,fje,fo2lg,fugac,fxi,iktmax,iopt,
     $ jflag,jsflag,kmax,kstep,kx1,kxt,mrmlra,modr,mosp,mospt,
     $ moph,mopht,narn1,narn2,nbasp,nbaspd,nbt,nbtmax,ncmpr,nct,
     $ nctmax,nelect,ngrn1,ngrn2,ngtmax,nhydr,nhydx,nllnmx,
     $ no2gaq,noptmx,noutpt,npt,nptmax,nrct,nrctmx,nstmax,ntabx,
     $ ntidmx,ntitl2,ntitld,ntitmx,nttyo,nxrn1,nxrn2,nxtmax,pe,ph,
     $ phmes,ppmwb,ppmwe,prcinf,press,prminf,qrho,qriinf,rho,rhowc,
     $ sidrph,sigmam,tdsgks,tdsglw,tempc,time1,uelem,ulinex,
     $ uphase,uplatm,ureac,uspec,usteq6,utitl2,utitld,uveeq6,
     $ vodrt,vosoct,wkgh2o,wodrt,wosoct,xbar,xbarlg,xi1)
c
c     This subroutine writes to TABX (the scrambled TAB file) using
c     a csv (comma separated value) format. A .csv file can be opened
c     by spreadsheets and other software that support data analysis
c     and plotting.
c
c     This subroutine is called by:
c
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer iktmax,kmax,nbtmax,nctmax,ngtmax,nllnmx,noptmx,nptmax,
     $ nrctmx,nstmax,ntidmx,ntitmx,nxtmax
c
      integer noutpt,ntabx,nttyo
c
      integer iopt(noptmx),jflag(nstmax),jsflag(nstmax),nbasp(nbtmax),
     $ nbaspd(nbtmax),ncmpr(2,nptmax)
c
      integer kstep,kx1,kxt,narn1,narn2,nbt,nct,nelect,ngrn1,ngrn2,
     $ nhydr,nhydx,no2gaq,npt,nrct,ntitl2,ntitld,nxrn1,nxrn2
c
      logical qrho,qriinf
c
      character(len=nllnmx) ulinex
      character(len=80) utitl2(ntitmx),utitld(ntidmx)
      character(len=48) uspec(nstmax)
      character(len=24) uphase(nptmax),ureac(nrctmx)
      character(len=8) uelem(nctmax)
      character(len=8) uplatm,usteq6,uveeq6
c
      real(8) acflg(nstmax),actlg(nstmax),afrc1(nrctmx),
     $ conclg(nstmax),ctb(nbtmax),cteaq(nctmax),fugac(ngtmax),
     $ modr(nrctmx),moph(nptmax),mopht(nptmax),mosp(nstmax),
     $ mospt(nstmax),ppmwb(nbtmax),ppmwe(nctmax),sidrph(nptmax),
     $ xbar(nstmax),xbarlg(nstmax)
c
      real(8) actw,aft1,alk,dvoso,dwoso,eh,fje,fo2lg,fxi,mrmlra,
     $ pe,ph,phmes,prcinf,press,prminf,rho,rhowc,sigmam,tdsgks,
     $ tdsglw,tempc,time1,vodrt,vosoct,wkgh2o,wodrt,wosoct,xi1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with variable global dimensioning.
c
      integer ngfmax,npmmax,npsmax,nxcmax
c
      SAVE ngfmax,npmmax,npsmax,nxcmax
c
      integer, dimension(:), allocatable :: ngasfx
c
      SAVE ngasfx
c
      character(len=24), dimension(:,:), allocatable :: ussphx
      character(len=24), dimension(:), allocatable :: ubaspx,
     $ ubasqx,usidrx,uphasx
      character(len=8), dimension(:), allocatable :: uelacx
c
      SAVE ubaspx,ubasqx,uelacx,uphasx,usidrx,ussphx
c
      real(8), dimension(:), allocatable :: acflbx,actlbx,afrcx,
     $ conlbx,ctbaqx,cteaqx,molrbx,molrex,mophx,ppmvbx,ppmvex,ppmwbx,
     $ ppmwex,sidrx,xfrac
c
      SAVE acflbx,actlbx,afrcx,conlbx,ctbaqx,cteaqx,molrbx,molrex,
     $ mophx,ppmvbx,ppmvex,ppmwbx,ppmwex,sidrx,xfrac
c
c-----------------------------------------------------------------------
c
c     Local variable declarations required to be saved between calls.
c
      integer ngft,npmt,npst,nxct
      integer nbtpr,nctpr,nqtpr
c
      SAVE ngft,npmt,npst,nxct
      SAVE nbtpr,nctpr,nqtpr
c
      integer ilstmx
c
      SAVE ilstmx
c
      logical qgftag,qpstag
c
      SAVE qgftag,qpstag
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,jj,j1,j2,j3,j4,jc,jlead,jline,jq,jsum,jt1,jt2,jt3,
     $ jthv1,jthv2,jthv3,jts,n,nb,nc,ng,ngt,nlim,np,nrc,nr1,nr2,ns,
     $ ns1,ns2
c
      integer ilnobl
c
      logical qnewph,qnewss,qpsext,qwrthd
c
      character(len=48) ux48
      character(len=16) ux16
      character(len=8) ux8,ux8a,ux8b,ultag1,ultag2,ultag3,ultags,
     $ ultgv1,ultgv2,ultgv3,unotav
c
      real(8) mx,sidrcv,tdays
c
      real(8) tlg
c
c-----------------------------------------------------------------------
c
      data unotav/'N/A'/
c
      data sidrcv/-10./
c
c-----------------------------------------------------------------------
c
c
c     Allocate or reallocate local work arrays as needed.
c
      if (.not.ALLOCATED(uphasx)) then
c
c       Local work arrays are not allocated. Note that only one array
c       is tested to see if all are allocated or not.
c
        ilstmx = (nllnmx/25) - 3
c
        npsmax = min(ilstmx,nptmax)
c
        ngfmax = min(ilstmx,ngtmax)
c
        nxcmax = nxtmax*iktmax
        nxcmax = min(ilstmx,nxcmax)
c
        npmmax = 8*kmax
        npmmax = min(ilstmx,npmmax)
c
        npmt = 0
        npst = 0
        ngft = 0
        nxct = 0
c
        qgftag = .false.
        qpstag = .false.
c
        nctpr = 0
        nqtpr = 0
        nbtpr = 0
c
        ALLOCATE(uphasx(npmmax))
        ALLOCATE(mophx(npmmax))
c
        ALLOCATE(uelacx(nctmax))
        ALLOCATE(cteaqx(nctmax),ppmwex(nctmax))
c
        ALLOCATE(ubasqx(nbtmax))
        ALLOCATE(ctbaqx(nbtmax),ppmwbx(nbtmax))
c
        if (qrho) then
          ALLOCATE(molrex(nctmax),ppmvex(nctmax))
          ALLOCATE(molrbx(nbtmax),ppmvbx(nbtmax))
        endif
c
        ALLOCATE(ubaspx(nbtmax))
        ALLOCATE(acflbx(nbtmax),actlbx(nbtmax),conlbx(nbtmax))
c
        ALLOCATE(afrcx(nrctmx))
c
        ALLOCATE(usidrx(npsmax))
        ALLOCATE(sidrx(npsmax))
c
        ALLOCATE(ngasfx(ngfmax))
c
        ALLOCATE(ussphx(2,nxcmax))
        ALLOCATE(xfrac(nxcmax))
c
c       Zero the contents of the local work arrays.
c
        do n = 1,nctmax
          uelacx(n) = ' '
          cteaqx(n) = 0.
          ppmwex(n) = 0.
        enddo
c
        do n = 1,nbtmax
          ubasqx(n) = ' '
          ctbaqx(n) = 0.
          ppmwbx(n) = 0.
        enddo
c
        if (qrho) then
          do n = 1,nctmax
            molrex(n) = 0.
            ppmvex(n) = 0.
          enddo
          do n = 1,nbtmax
            molrbx(n) = 0.
            ppmvbx(n) = 0.
          enddo
        endif
c
        do n = 1,nbtmax
          ubaspx(n) = ' '
          acflbx(n) = 0.
          actlbx(n) = 0.
          conlbx(n) = 0.
        enddo
c
        do n = 1,nrctmx
          afrcx(n) = 0.
        enddo
c
        do n = 1,npmmax
          uphasx(n) = ' '
          mophx(n) = 0.
        enddo
c
        do n = 1,npsmax
          usidrx(n) = ' '
          sidrx(n) = 0.
        enddo
c
        do n = 1,ngfmax
          ngasfx(n) = 0
        enddo
c
        do n = 1,nxcmax
          ussphx(1,n) = ' '
          ussphx(2,n) = ' '
          xfrac(n) = 0.
        enddo
c
c       Construct the local element list.
c
        nctpr = 0
        do nc = 1,nct
          if (uelem(nc)(1:2).ne.'H ' .and. uelem(nc)(1:2) .ne. 'O ')
     $      then
            nctpr = nctpr + 1
            uelacx(nctpr) = uelem(nc)
          endif
        enddo
c
c       Construct the local solute basis species list.
c
        nqtpr = 0
        do nb = 1,nbt
          ns1 = nbaspd(nb)
          ns2 = nbasp(nb)
          if (ns1.ge.narn1 .and. ns1.le.narn2) then
            if (jsflag(ns2).lt.2 .and. jflag(ns2).ne.30) then
              if (ns1.ne.narn1 .and. ns1.ne.nelect .and.
     $          ns1.ne.no2gaq) then
                nqtpr = nqtpr + 1
                ubasqx(nqtpr) = uspec(ns1)
              endif
            endif
          endif
        enddo
c
c       Construct the local complete basis species list.
c
        nbtpr = 0
        do nb = 1,nbt
          ns1 = nbaspd(nb)
          ns2 = nbasp(nb)
          if (ns1.ge.narn1 .and. ns1.le.narn2) then
            if (jsflag(ns2).lt.2 .and. jflag(ns2).ne.30) then
              nbtpr = nbtpr + 1
              ubaspx(nbtpr) = uspec(ns1)
            endif
          endif
        enddo
c
c       Construct the local gas species index map list.
c
        ngft = 0
        do ns = ngrn1,ngrn2
          if (jsflag(ns) .lt. 2) then
            if (ngft .ge. ngfmax) then
c
c             Prepare to write a "More omitted" label
c
              qgftag = .true.
              go to 100
            endif
            ngft = ngft + 1
            ngasfx(ngft) = ns
          endif
        enddo
  100   continue
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write output tables of some selected parameters. These are
c     identified by short strings, commonly letters or letters and
c     numbers. They are written on the TABX file in scrambled order
c     (i.e., the lines of any given tables are interspersed with those
c     of the other tables). Each line is marked at the beginning with
c     an id string marking the table to which it belongs. The lines
c     will later be descrambled so that all lines belonging to a given
c     table are contiguous.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write Table A:
c
c       Run title (first 15 lines max)
c       Code version identification
c       Data file title (first 10 lines max)
c
c     Note: the original primary title (utitl1) has already been rolled
c     over to the current secondary title (utitl2).
c
      qwrthd = .false.
      if (kstep.eq.0 .and. iopt(18).eq.0) then
c
        ultag1 = 'A       '
        jt1 = len_trim(ultag1)
        ultags = ultag1
        jts = jt1
c
c       Write the table header.
c
        write (ntabx,1000) ultag1(1:jt1)
 1000   format(a,',Table A,Title information (filtered),')
c
        qwrthd = .true.
c
        write (ntabx,'(a,",")') ultag1(1:jt1)
        j2 = ilnobl(uveeq6)
        j3 = ilnobl(usteq6)
        j4 = ilnobl(uplatm)
        write (ntabx,1020) ultag1(1:jt1),uveeq6(1:j2),usteq6(1:j3),
     $  uplatm(1:j4)
 1020   format(a,',Running EQ6 from EQ3/6-V',a,' (',a,') for ',a,',')
        write (ntabx,'(a,",---,")') ultag1(1:jt1)
        write (ntabx,'(a,",")') ultag1(1:jt1)
c
c       Write the table data.
c
        nlim = min(ntitl2,15)
        do n = 1,nlim
          jj = min(80,nllnmx)
          ulinex = utitl2(n)(1:jj)
          j1 = len_trim(ulinex)
c
c         Check for double quotes and commas.
c
          jq = index(ulinex,'"')
          jc = index(ulinex,',')
          jsum = jq + jc
          if (jsum .gt. 0) then
c
c           Replace any double quotes in the title text by
c           single quotes.
c
            if (jq .gt. 0) then
c
  110         ulinex(jq:jq) = "'"
              jq = index(ulinex,'"')
              if (jq .gt. 0) go to 110
c
            endif
c
c           Replace any commas in the title text by semicolons.
c
            jc = index(ulinex,',')
            if (jc .gt. 0) then
c
  120         ulinex(jc:jc) = ";"
              jc = index(ulinex,',')
              if (jc .gt. 0) go to 120
c
            endif
          endif
c
          write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
        enddo
        write (ntabx,'(a,",---,")') ultag1(1:jt1)
        write (ntabx,'(a,",")') ultag1(1:jt1)
c
        nlim = min(ntitld,10)
        do n = 1,nlim
          jj = min(80,nllnmx)
          ulinex = utitld(n)(1:jj)
          j1 = len_trim(ulinex)
c
c         Check for double quotes and commas.
c
          jq = index(ulinex,'"')
          jc = index(ulinex,',')
          jsum = jq + jc
          if (jsum .gt. 0) then
c
c           Replace any double quotes in the title text by
c           single quotes.
c
            if (jq .gt. 0) then
c
  130         ulinex(jq:jq) = "'"
              jq = index(ulinex,'"')
              if (jq .gt. 0) go to 130
c
            endif
c
c           Replace any commas in the title text by semicolons.
c
            jc = index(ulinex,',')
            if (jc .gt. 0) then
c
  140         ulinex(jc:jc) = ";"
              jc = index(ulinex,',')
              if (jc .gt. 0) go to 140
c
            endif
          endif
c
          write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
        enddo
        write (ntabx,'(a,",---,")') ultag1(1:jt1)
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write Table B1:
c
c       Xi
c       t(days)
c       Temperature(C)
c       Pressure(bars)
c       pH
c       pmH
c       log fO2
c       Eh(v)
c       pe
c       Activity of water
c
c     Prepare the data.
c
      if (.not.qriinf) then
        tdays = time1/86400.
      endif
c
      if (eh .le. -99999.) eh = -99999.
      if (pe .le. -99999.) pe = -99999.
      if (aft1 .ge. 99999.) aft1 = 99999.
c
      ultag1 = 'B1      '
      jt1 = len_trim(ultag1)
c
      if (qwrthd) then
c
c       Write an end marker for the previous table.
c
        write (ntabx,'(a,",EndTable:,",a,",")')
     $  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1
c
c       Write two empty lines before the table header.
c
        do i = 1,2
          write (ntabx,'(a,",")') ultag1(1:jt1)
        enddo
c
c       Write the table header.
c
        write (ntabx,1040) ultag1(1:jt1)
 1040   format(a,',Table B1,Miscellaneous parameters I,')
c
c       Write the column header.
c
        ulinex = 'Xi,t(days),Temp(C),Press(bars),pH,pmH,'
        j1 = len_trim(ulinex)
        ux48 = 'log fO2,Eh(v),pe,aw'
        j2 = len_trim(ux48)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux48(1:j2)
        j1 = jsum
c
        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
      endif
c
c     Write the table data.
c
      ulinex = ''
      j1 = 0
c
      write (ux16,'(1pg11.4)') xi1
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') tdays
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(f11.4)') tempc
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(f11.4)') press
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(f11.4)') ph
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(f11.4)') phmes
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(f11.4)') fo2lg
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(f11.4)') eh
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(f11.4)') pe
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(f11.4)') actw
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      j1 = jsum
c
      write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write Table B2:
c
c       Xi
c       t(days)
c       Mass of solvent H2O, kg
c       Alkalinity(eq/kg H2O)
c       Sigma m
c       Ionic strength
c       Ionic asymmetry
c       TDS, g/kg.sol
c       TDS, g/L
c       Aqueous solution density, g/L
c       Molarity/molality ratio
c
c     Prepare the data. None to prepare.
c
      ultag1 = 'B2      '
      jt1 = len_trim(ultag1)
c
      if (qwrthd) then
c
c       Write an end marker for the previous table.
c
        write (ntabx,'(a,",EndTable:,",a,",")')
     $  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1
c
c       Write two empty lines before the table header.
c
        do i = 1,2
          write (ntabx,'(a,",")') ultag1(1:jt1)
        enddo
c
c       Write the table header.
c
        write (ntabx,1060) ultag1(1:jt1)
 1060   format(a,',Table B2,Miscellaneous parameters II,')
c
c       Write the column header.
c
        ulinex = 'Xi,t(days),H2O(kg),Alk(eq/kg.H2O),Signam(m),I(m),'
        j1 = len_trim(ulinex)
        ux48 = 'J(m),TDS(g/kg.sol),TDS(g/L),density(g/L),M/m'
        j2 = len_trim(ux48)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux48(1:j2)
        j1 = jsum
c
        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
      endif
c
c     Write the table data.
c
      ulinex = ''
      j1 = 0
c
      write (ux16,'(1pg11.4)') xi1
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') tdays
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') wkgh2o
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') alk
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') sigmam
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') fxi
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') fje
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(f11.4)') tdsgks
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(f11.4)') tdsglw
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(f11.4)') rhowc
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(f11.4)') mrmlra
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      j1 = jsum
c
      write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write Table C1:
c
c       Xi
c       t(days)
c       Molalities of dissolved elements
c
c     and Table C2:
c
c       Xi
c       t(days)
c       ppm (mg/kg.sol) of dissolved elements
c
c     Prepare the data.
c
      nctpr = 0
      do nc = 1,nct
        if (uelem(nc)(1:2).ne.'H ' .and. uelem(nc)(1:2) .ne. 'O ') then
          nctpr = nctpr + 1
          uelacx(nctpr) = uelem(nc)
          cteaqx(nctpr) = cteaq(nc)
          ppmwex(nctpr) = ppmwe(nc)
        endif
      enddo
c
      ultag1 = 'C1      '
      ultag2 = 'C2      '
      jt1 = len_trim(ultag1)
      jt2 = len_trim(ultag2)
c
      if (qwrthd) then
c
c       Write an end marker for the table previous to C1.
c
        write (ntabx,'(a,",EndTable:,",a,",")')
     $  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1
c
c       Write an end marker for the table previous to C2.
c
        write (ntabx,'(a,",EndTable:,",a,",")')
     $  ultag2(1:jt2),ultags(1:jts)
        ultags = ultag2
        jts = jt2
c
c       Write two empty lines before the table header.
c
        do i = 1,2
          write (ntabx,'(a,",")') ultag1(1:jt1)
          write (ntabx,'(a,",")') ultag2(1:jt2)
        enddo
c
c       Write the table header.
c
        write (ntabx,1100) ultag1(1:jt1)
 1100   format(a,',Table C1,Dissolved elements(molality),')
        write (ntabx,1110) ultag2(1:jt2)
 1110   format(a,',Table C2,Dissolved elements(ppm: mg/kg.sol),')
c
c       Write the column headers.
c
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)
        do n = 1,nctpr
          j1 = j1 + 1
          ulinex(j1:j1) = ","
          j2 = len_trim(uelacx(n))
          jsum = j1 + j2
          ulinex(j1 + 1:jsum) = uelacx(n)(1:j2)
          j1 = jsum
        enddo
c
        jline = j1
        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
        write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
      endif
c
c     Write the table data. First, the common data (Xi, t( days)).
c
      ulinex = ''
      j1 = 0
c
      write (ux16,'(1pg11.4)') xi1
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') tdays
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      j1 = jsum
      jlead = j1
c
c     Write the data specific to Table C1.
c
      do n = 1,nctpr
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(1pg11.4)') cteaqx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
      enddo
c
      write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
c
c     Write the data specific to Table C2.
c
      ulinex(jlead + 1:j1) = ' '
      j1 = jlead
      do n = 1,nctpr
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(1pg11.4)') ppmwex(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
      enddo
c
      write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qrho) then
c
c       Write Table C3:
c
c         Xi
c         t(days)
c         Molarities of dissolved elements
c
c       and Table C4:
c
c         Xi
c         t(days)
c         ppm (mg/L) of dissolved elements
c
c       Prepare the data.
c
        nctpr = 0
        do nc = 1,nct
          if (uelem(nc)(1:2).ne.'H ' .and. uelem(nc)(1:2) .ne. 'O ')
     $      then
            nctpr = nctpr + 1
            uelacx(nctpr) = uelem(nc)
            molrex(nctpr) = cteaq(nc)*mrmlra
            ppmvex(nctpr) = ppmwe(nc)*rho
          endif
        enddo
c
        ultag1 = 'C3       '
        ultag2 = 'C4       '
        jt1 = len_trim(ultag1)
        jt2 = len_trim(ultag2)
c
        if (qwrthd) then
c
c         Write an end marker for the table previous to C3.
c
          write (ntabx,'(a,",EndTable:,",a,",")')
     $    ultag1(1:jt1),ultags(1:jts)
          ultags = ultag1
          jts = jt1
c
c         Write an end marker for the table previous to C4.
c
          write (ntabx,'(a,",EndTable:,",a,",")')
     $    ultag2(1:jt2),ultags(1:jts)
          ultags = ultag2
          jts = jt2
c
c         Write two empty lines before the table header.
c
          do i = 1,2
            write (ntabx,'(a,",")') ultag1(1:jt1)
            write (ntabx,'(a,",")') ultag2(1:jt2)
          enddo
c
c         Write the table header.
c
          write (ntabx,1120) ultag1(1:jt1)
 1120     format(a,',Table C3,Dissolved elements(molarity),')
          write (ntabx,1130) ultag2(1:jt2)
 1130     format(a,',Table C4,Dissolved elements(ppm: mg/L),')
c
c         Write the column headers.
c
          ulinex = 'Xi,t(days)'
          j1 = len_trim(ulinex)
          do n = 1,nctpr
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            j2 = len_trim(uelacx(n))
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = uelacx(n)(1:j2)
            j1 = jsum
          enddo
c
          jline = j1
          write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
          write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
        endif
c
c       Write the table data. First, the common data (Xi, t( days)).
c
        ulinex = ''
        j1 = 0
c
        write (ux16,'(1pg11.4)') xi1
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        jsum = jsum + 1
        ulinex(jsum:jsum) = ','
        j1 = jsum
c
        write (ux16,'(1pg11.4)') tdays
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
        jlead = j1
c
c       Write the data specific to Table C3.
c
        do n = 1,nctpr
          j1 = j1 + 1
          ulinex(j1:j1) = ","
          write (ux16,'(1pg11.4)') molrex(n)
          call lejust(ux16)
          j2 = len_trim(ux16)
          jsum = j1 + j2
          ulinex(j1 + 1:jsum) = ux16(1:j2)
          j1 = jsum
        enddo
c
        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
c
c       Write the data specific to Table C4.
c
        ulinex(jlead + 1:j1) = ' '
        j1 = jlead
        do n = 1,nctpr
          j1 = j1 + 1
          ulinex(j1:j1) = ","
          write (ux16,'(1pg11.4)') ppmvex(n)
          call lejust(ux16)
          j2 = len_trim(ux16)
          jsum = j1 + j2
          ulinex(j1 + 1:jsum) = ux16(1:j2)
          j1 = jsum
        enddo
c
        write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write Table D1:
c
c       Xi
c       t(days)
c       Total molalities of solute basis species
c
c     and Table D2:
c
c       Xi
c       t(days)
c       Total ppm (mg/kg.sol) of solute basis species
c
c     Prepare the data.
c
      n = 0
      do nb = 1,nbt
        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)
        if (ns1.ge.narn1 .and. ns1.le.narn2) then
          if (jsflag(ns2).lt.2 .and. jflag(ns2).ne.30) then
            if (ns1.ne.narn1 .and. ns1.ne.nelect .and.
     $        ns1.ne.no2gaq) then
              n = n + 1
              ctbaqx(n) = ctb(nb)
              ppmwbx(n) = ppmwb(nb)
            endif
          endif
        endif
      enddo
      if (n .ne. nqtpr) then
c
c       Programming error trap.
c
        write (ux8a,'(i8)') n
        write (ux8b,'(i8)') nqtpr
        call lejust(ux8a)
        call lejust(ux8b)
        j1 = len_trim(ux8a)
        j2 = len_trim(ux8b)
        write (noutpt,1140) ux8a(1:j1),ux8b(1:j2)
        write (nttyo,1140) ux8a(1:j1),ux8b(1:j2)
 1140   format(/' * Error - (EQ6/wrtabc) Programming error trap: ',
     $  'The solute',/7x,'basis species counts for Tables D1 and D2',
     $  ' do not match',/7x,'(',a,' versus ',a,').')
        stop
      endif
c
      ultag1 = 'D1      '
      ultag2 = 'D2      '
      jt1 = len_trim(ultag1)
      jt2 = len_trim(ultag2)
c
      if (qwrthd) then
c
c       Write an end marker for the table previous to D1.
c
        write (ntabx,'(a,",EndTable:,",a,",")')
     $  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1
c
c       Write an end marker for the table previous to D2.
c
        write (ntabx,'(a,",EndTable:,",a,",")')
     $  ultag2(1:jt2),ultags(1:jts)
        ultags = ultag2
        jts = jt2
c
c       Write two empty lines before the table header.
c
        do i = 1,2
          write (ntabx,'(a,",")') ultag1(1:jt1)
          write (ntabx,'(a,",")') ultag2(1:jt2)
        enddo
c
c       Write the table header.
c
        write (ntabx,1220) ultag1(1:jt1)
 1220   format(a,',Table D1,Solute basis species(total molality),')
        write (ntabx,1230) ultag2(1:jt2)
 1230   format(a,',Table D2,Solute basis species',
     $  '(total ppm: mg/kg.sol),')
c
c       Write the column headers.
c
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)
        do n = 1,nqtpr
          j1 = j1 + 1
          ulinex(j1:j1) = ","
          j2 = len_trim(ubasqx(n))
          jsum = j1 + j2
          ulinex(j1 + 1:jsum) = ubasqx(n)(1:j2)
          j1 = jsum
        enddo
c
        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
        write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
      endif
c
c     Write the table data. First, the common data (Xi, t( days)).
c
      ulinex = ''
      j1 = 0
c
      write (ux16,'(1pg11.4)') xi1
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') tdays
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      j1 = jsum
      jlead = j1
c
c     Write the data specific to Table D1.
c
      do n = 1,nqtpr
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(1pg11.4)') ctbaqx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
      enddo
c
      write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
c
c     Write the data specific to Table D2.
c
      ulinex(jlead + 1:j1) = ' '
      j1 = jlead
      do n = 1,nqtpr
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(f11.4)') ppmwbx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
      enddo
c
      write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qrho) then
c
c       Write Table D3:
c
c         Xi
c         t(days)
c         Total molarities of solute basis species
c
c       and Table D4:
c
c         Xi
c         t(days)
c         Total ppm (mg/L) solute file basis species
c
c       Prepare the data.
c
        n = 0
        do nb = 1,nbt
          ns1 = nbaspd(nb)
          ns2 = nbasp(nb)
          if (ns1.ge.narn1 .and. ns1.le.narn2) then
            if (jsflag(ns2).lt.2 .and. jflag(ns2).ne.30) then
              if (ns1.ne.narn1 .and. ns1.ne.nelect .and.
     $          ns1.ne.no2gaq) then
                n = n + 1
                molrbx(n) = ctb(nb)*mrmlra
                ppmvbx(n) = ppmwb(nb)*rho
              endif
            endif
          endif
        enddo
        if (n .ne. nqtpr) then
c
c         Programming error trap.
c
          write (ux8a,'(i8)') n
          write (ux8b,'(i8)') nqtpr
          call lejust(ux8a)
          call lejust(ux8b)
          j1 = len_trim(ux8a)
          j2 = len_trim(ux8b)
          write (noutpt,1235) ux8a(1:j1),ux8b(1:j2)
          write (nttyo,1235) ux8a(1:j1),ux8b(1:j2)
 1235     format(/' * Error - (EQ6/wrtabc) Programming error trap: ',
     $    'The solute',/7x,'basis species counts for Tables D3 and D4',
     $    ' do not match',/7x,'(',a,' versus ',a,').')
          stop
        endif
c
        ultag1 = 'D3      '
        ultag2 = 'D4      '
        jt1 = len_trim(ultag1)
        jt2 = len_trim(ultag2)
c
        if (qwrthd) then
c
c         Write an end marker for the table previous to D3.
c
          write (ntabx,'(a,",EndTable:,",a,",")')
     $    ultag1(1:jt1),ultags(1:jts)
          ultags = ultag1
          jts = jt1
c
c         Write an end marker for the table previous to D4.
c
          write (ntabx,'(a,",EndTable:,",a,",")')
     $    ultag2(1:jt2),ultags(1:jts)
          ultags = ultag2
          jts = jt2
c
c         Write two empty lines before the table header.
c
          do i = 1,2
            write (ntabx,'(a,",")') ultag1(1:jt1)
            write (ntabx,'(a,",")') ultag2(1:jt2)
          enddo
c
c         Write the table header.
c
          write (ntabx,1240) ultag1(1:jt1)
 1240     format(a,',Table D3,Solute basis species(total molarity),')
          write (ntabx,1250) ultag2(1:jt2)
 1250     format(a,',Table D4,Solute basis species(total ppm: mg/L),')
c
c         Write the column headers.
c
          ulinex = 'Xi,t(days)'
          j1 = len_trim(ulinex)
          do n = 1,nqtpr
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            j2 = len_trim(ubasqx(n))
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = ubasqx(n)(1:j2)
            j1 = jsum
          enddo
c
          write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
          write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
        endif
c
c       Write the table data. First, the common data (Xi, t( days)).
c
        ulinex = ''
        j1 = 0
c
        write (ux16,'(1pg11.4)') xi1
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        jsum = jsum + 1
        ulinex(jsum:jsum) = ','
        j1 = jsum
c
        write (ux16,'(1pg11.4)') tdays
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
        jlead = j1
c
c       Write the data specific to Table D3.
c
        do n = 1,nqtpr
          j1 = j1 + 1
          ulinex(j1:j1) = ","
          write (ux16,'(1pg11.4)') molrbx(n)
          call lejust(ux16)
          j2 = len_trim(ux16)
          jsum = j1 + j2
          ulinex(j1 + 1:jsum) = ux16(1:j2)
          j1 = jsum
        enddo
c
        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
c
c       Write the data specific to Table D4.
c
        ulinex(jlead + 1:j1) = ' '
        j1 = jlead
        do n = 1,nqtpr
          j1 = j1 + 1
          ulinex(j1:j1) = ","
          write (ux16,'(f11.4)') ppmvbx(n)
          call lejust(ux16)
          j2 = len_trim(ux16)
          jsum = j1 + j2
          ulinex(j1 + 1:jsum) = ux16(1:j2)
          j1 = jsum
        enddo
c
        write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write Table E1:
c
c       Xi
c       t(days)
c       Log true molalities of basis species (log mole fraction for H2O,
c       (log fugacity for aqueous O2(g))
c
c     Table E2:
c
c       Xi
c       t(days)
c       Log activities of basis species (log fugacity for aqueous O2(g))
c
c     and Table E3:
c
c       Xi
c       t(days)
c       Log true activity coefficients of basis species
c
c     Prepare the data.
c
      n = 0
      do nb = 1,nbt
        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)
        if (ns1.ge.narn1 .and. ns1.le.narn2) then
          if (jsflag(ns2).lt.2 .and. jflag(ns2).ne.30) then
            n = n + 1
            if (ns1 .eq. narn1) then
              conlbx(n) = xbarlg(ns1)
              actlbx(n) = actlg(ns1)
              acflbx(n) = acflg(ns1)
            elseif (ns1 .eq. no2gaq) then
              conlbx(n) = fo2lg
              actlbx(n) = fo2lg
              acflbx(n) = acflg(ns1)
            else
              conlbx(n) = conclg(ns1)
              actlbx(n) = actlg(ns1)
              acflbx(n) = acflg(ns1)
            endif
          endif
        endif
      enddo
      if (n .ne. nbtpr) then
c
c       Programming error trap.
c
        write (ux8a,'(i8)') n
        write (ux8b,'(i8)') nbtpr
        call lejust(ux8a)
        call lejust(ux8b)
        j1 = len_trim(ux8a)
        j2 = len_trim(ux8b)
        write (noutpt,1260) ux8a(1:j1),ux8b(1:j2)
        write (nttyo,1260) ux8a(1:j1),ux8b(1:j2)
 1260   format(/' * Error - (EQ6/wrtabc) Programming error trap: ',
     $  'The local complete',/7x,'basis species counts for Tables E1,',
     $  ' E2, and E3 do not match',/7x,'(',a,' versus ',a,').')
        stop
      endif
c
      ultag1 = 'E1      '
      ultag2 = 'E2      '
      ultag3 = 'E3      '
      jt1 = len_trim(ultag1)
      jt2 = len_trim(ultag2)
      jt3 = len_trim(ultag3)
c
      if (qwrthd) then
c
c       Write an end marker for the table previous to E1.
c
        write (ntabx,'(a,",EndTable:,",a,",")')
     $  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1
c
c       Write an end marker for the table previous to E2.
c
        write (ntabx,'(a,",EndTable:,",a,",")')
     $  ultag2(1:jt2),ultags(1:jts)
        ultags = ultag2
        jts = jt2
c
c       Write an end marker for the table previous to E3.
c
        write (ntabx,'(a,",EndTable:,",a,",")')
     $  ultag3(1:jt3),ultags(1:jts)
        ultags = ultag3
        jts = jt3
c
c       Write two empty lines before the table header.
c
        do i = 1,2
          write (ntabx,'(a,",")') ultag1(1:jt1)
          write (ntabx,'(a,",")') ultag2(1:jt2)
          write (ntabx,'(a,",")') ultag3(1:jt3)
        enddo
c
c       Write the table header.
c
        write (ntabx,1270) ultag1(1:jt1)
 1270   format(a,',Table E1,Basis species(log true molality;',
     $  ' log mole fraction for H2O, log fugacity for O2(g)),')
        write (ntabx,1280) ultag2(1:jt2)
 1280   format(a,',Table E2,Basis species(log activity;',
     $  ' log fugacity for O2(g)),')
        write (ntabx,1290) ultag3(1:jt3)
 1290   format(a,
     $  ',Table E3,Basis species(log activity coefficient),')
c
c       Write the column headers.
c
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)
        do n = 1,nbtpr
          j1 = j1 + 1
          ulinex(j1:j1) = ","
          j2 = len_trim(ubaspx(n))
          jsum = j1 + j2
          ulinex(j1 + 1:jsum) = ubaspx(n)(1:j2)
          j1 = jsum
        enddo
c
        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
        write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
        write (ntabx,'(a,",",a,",")') ultag3(1:jt3),ulinex(1:j1)
      endif
c
c     Write the table data. First, the common data (Xi, t( days)).
c
      ulinex = ''
      j1 = 0
c
      write (ux16,'(1pg11.4)') xi1
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') tdays
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      j1 = jsum
      jlead = j1
c
c     Write the data specific to Table E1.
c
      do n = 1,nbtpr
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(1pg11.4)') conlbx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
      enddo
c
      write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
c
c     Write the data specific to Table E2.
c
      ulinex(jlead + 1:j1) = ' '
      j1 = jlead
      do n = 1,nbtpr
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(f11.4)') actlbx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
      enddo
c
      write (ntabx,'(a,",",a,",")') ultag2(1:jt2),ulinex(1:j1)
c
c     Write the data specific to Table E3.
c
      ulinex(jlead + 1:j1) = ' '
      j1 = jlead
      do n = 1,nbtpr
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(f11.4)') acflbx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
      enddo
c
      write (ntabx,'(a,",",a,",")') ultag3(1:jt3),ulinex(1:j1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write Table J:
c
c       Xi
c       t(days)
c       Moles of reactants destroyed/created
c
c     Prepare the data. None to prepare.
c
      ultag1 = 'J       '
      jt1 = len_trim(ultag1)
c
      if (qwrthd) then
c
c       Write an end marker for the previous table.
c
        write (ntabx,'(a,",EndTable:,",a,",")')
     $  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1
c
c       Write two empty lines before the table header.
c
        do i = 1,2
          write (ntabx,'(a,",")') ultag1(1:jt1)
        enddo
c
c       Write the table header.
c
        write (ntabx,1340) ultag1(1:jt1)
 1340   format(a,',Table J,Moles of reactants destroyed/created,')
c
c       Write the column header.
c
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)
        do n = 1,nrct
          j1 = j1 + 1
          ulinex(j1:j1) = ","
          j2 = len_trim(ureac(n))
          jsum = j1 + j2
          ulinex(j1 + 1:jsum) = ureac(n)(1:j2)
          j1 = jsum
        enddo
c
        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
      endif
c
c     Write the table data.
c
      ulinex = ''
      j1 = 0
c
      write (ux16,'(1pg11.4)') xi1
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') tdays
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      j1 = jsum
c
      do n = 1,nrct
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(1pg11.4)') modr(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
      enddo
c
      write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write Table K:
c
c       Xi
c       t(days)
c       Total Affinity, kcal
c       Affinities of reactants, kcal
c
c     Prepare the data.
c
      do n = 1,nrct
        afrcx(n) = min(99999.,afrc1(n))
      enddo
c
      ultag1 = 'K       '
      jt1 = len_trim(ultag1)
c
      if (qwrthd) then
c
c       Write an end marker for the previous table.
c
        write (ntabx,'(a,",EndTable:,",a,",")')
     $  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1
c
c       Write two empty lines before the table header.
c
        do i = 1,2
          write (ntabx,'(a,",")') ultag1(1:jt1)
        enddo
c
c       Write the table header.
c
        write (ntabx,1360) ultag1(1:jt1)
 1360   format(a,',Table K,Affinities of reactants (kcal),')
c
c       Write the column header.
c
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)
        ux48 = ',Total'
        j2 = len_trim(ux48)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux48(1:j2)
        j1 = jsum
        do n = 1,nrct
          j1 = j1 + 1
          ulinex(j1:j1) = ","
          j2 = len_trim(ureac(n))
          jsum = j1 + j2
          ulinex(j1 + 1:jsum) = ureac(n)(1:j2)
          j1 = jsum
        enddo
c
        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
      endif
c
c     Write the table data.
c
      ulinex = ''
      j1 = 0
c
      write (ux16,'(1pg11.4)') xi1
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') tdays
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(f11.4)') aft1
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      j1 = jsum
c
      do n = 1,nrct
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(f11.4)') afrcx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
      enddo
c
      write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write Table P:
c
c       Xi
c       t(days)
c       Moles of product phases
c
c     Note: in these tables, moles of product phase refers to what
c     is present in the Equilibrium System for Titration and Closed
c     system models, and to cumulative totals for the Fluid-Centered
c     Flow-Through system.
c
c     Prepare the data.
c
      do n = 1,npmt
        mophx(n) = 0.
      enddo
c
      qnewph = .false.
      if (iopt(1) .ne. 2) then
c
c       Closed or titration system.
c
        do np = 2,npt
c**       if (moph(np).gt.0. .and. uphase(np)(1:5).ne."fix_f") then
          if (moph(np) .gt. 0.) then
            do n = 1,npmt
              if (uphase(np) .eq. uphasx(n)) then
c
c               Have found a product mineral in the existing list.
c
                mophx(n) = moph(np)
                go to 310
              endif
            enddo
c
c           Did not find a current product mineral in the existing
c           list. Add it to the list.
c
            qnewph = .true.
            npmt = npmt + 1
            n = npmt
            uphasx(n) = uphase(np)
            mophx(n) = moph(np)
          endif
  310     continue
        enddo
      else
c
c       Fluid-centered flow-through system.
c
        do np = 2,npt
c**       if (moph(np).gt.0. .and. uphase(np)(1:5).ne."fix_f") then
          if (moph(np) .gt. 0.) then
            do n = 1,npmt
              if (uphase(np) .eq. uphasx(n)) then
c
c               Have found a product mineral in the existing list.
c
                mophx(n) = mopht(np)
                go to 320
              endif
            enddo
c
c           Did not find a current product mineral in the existing
c           list. Add it to the list.
c
            qnewph = .true.
            npmt = npmt + 1
            n = npmt
            uphasx(n) = uphase(np)
            mophx(n) = mopht(np)
          endif
  320     continue
        enddo
      endif
c
      ultag1 = 'P       '
      jt1 = len_trim(ultag1)
c
      if (qwrthd) then
c
c       Write an end marker for the previous table.
c
        write (ntabx,'(a,",EndTable:,",a,",")')
     $  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1
c
c       Write two empty lines before the table header.
c
        do i = 1,2
          write (ntabx,'(a,",")') ultag1(1:jt1)
        enddo
c
c       Write the table header.
c
        if (iopt(1) .ne. 2) then
          write (ntabx,1370) ultag1(1:jt1)
 1370     format(a,',Table P,Moles of product minerals,')
        else
          write (ntabx,1380) ultag1(1:jt1)
 1380     format(a,',Table P,Moles of product minerals (cumulative),')
        endif
      endif
c
      if (qnewph .or. kstep.le.0) then
c
c       Write the column header. This is an example of a variable
c       header.
c
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)
        do n = 1,npmt
          j1 = j1 + 1
          ulinex(j1:j1) = ","
          j2 = len_trim(uphasx(n))
          jsum = j1 + j2
          ulinex(j1 + 1:jsum) = uphasx(n)(1:j2)
          j1 = jsum
        enddo
c
        ultgv1(1:6) = ultag1(1:6)
        ultgv1(7:8) = 'vh'
        jthv1 = len_trim(ultgv1)
        write (ntabx,'(a,",",a,",")') ultgv1(1:jthv1),ulinex(1:j1)
      endif
c
c     Write the table data.
c
      ulinex = ''
      j1 = 0
c
      write (ux16,'(1pg11.4)') xi1
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') tdays
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      j1 = jsum
c
      do n = 1,npmt
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(1pg11.4)') mophx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
      enddo
c
      write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write Table Q:
c
c       Xi
c       t(days)
c       Saturation indices of potential product phases
c
c     Note: in these tables, the saturation index (SI) is represented
c     by the function log(Q/K). Values are contained in the sidrph
c     array. A potential product phase is considered to be one which
c     at some point in the run has a saturation index above a certain
c     threshhold (sidrcv).
c
c     Prepare the data.
c
      do n = 1,npst
        sidrx(n) = 0.
      enddo
c
      qpsext = .false.
c
      do np = 2,npt
        if (uphase(np)(1:5).ne."fix_f") then
          do n = 1,npst
            if (uphase(np) .eq. usidrx(n)) then
c
c             Have found a potential product phase in the existing list.
c
              sidrx(n) = sidrph(np)
              go to 350
            endif
          enddo
c
          if (sidrph(np).gt.sidrcv) then
c
c           Did not find a potential product phase in the existing
c           list. Add a phase, write a "More omitted" tag, or skip if
c           a tag has been previously written.
c
            if (qpstag) then
c
c             The list has previously maxed out and a "More omitted"
c             tag has been written. Do nothing.
c
              continue
c
            elseif (npst .eq. npsmax) then
c
c             The list has previously maxed out, but a tag has
c             not been written. Prepare to write a tag.
c
              qpsext = .true.
              qpstag = .true.
c
            else
c
c             Add a new potential product phase to the list.
c
              qpsext = .true.
              npst = npst + 1
              n = npst
              usidrx(n) = uphase(np)
              sidrx(n) = sidrph(np)
            endif
          endif
        endif
  350   continue
      enddo
c
      ultag1 = 'Q       '
      jt1 = len_trim(ultag1)
c
      if (qwrthd) then
c
c       Write an end marker for the previous table.
c
        write (ntabx,'(a,",EndTable:,",a,",")')
     $  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1
c
c       Write two empty lines before the table header.
c
        do i = 1,2
          write (ntabx,'(a,",")') ultag1(1:jt1)
        enddo
c
c       Write the table header.
c
        write (ntabx,1392) ultag1(1:jt1)
 1392   format(a,',Table Q,Saturation indices of potential',
     $  ' product phases,')
      endif
c
      if (qpsext .or. kstep.le.0) then
c
c       Write the column header. This is an example of a variable
c       header.
c
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)
        do n = 1,npst
          j1 = j1 + 1
          ulinex(j1:j1) = ","
          j2 = len_trim(usidrx(n))
          jsum = j1 + j2
          ulinex(j1 + 1:jsum) = usidrx(n)(1:j2)
          j1 = jsum
        enddo
        if (qpstag) then
c
c         Have more phases than can fit.
c
          j1 = j1 + 1
          ulinex(j1:j1) = ","
          ux16 = 'More omitted'
          j2 = len_trim(ux16)
          jsum = j1 + j2
          ulinex(j1 + 1:jsum) = ux16(1:j2)
          j1 = jsum
        endif
c
        ultgv1(1:6) = ultag1(1:6)
        ultgv1(7:8) = 'vh'
        jthv1 = len_trim(ultgv1)
        write (ntabx,'(a,",",a,",")') ultgv1(1:jthv1),ulinex(1:j1)
      endif
c
c     Write the table data.
c
      ulinex = ''
      j1 = 0
c
      write (ux16,'(1pg11.4)') xi1
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') tdays
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      j1 = jsum
c
      do n = 1,npst
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        write (ux16,'(f11.4)') sidrx(n)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
      enddo
c
      write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write Table T:
c
c       Xi
c       t(days)
c       Fugacities of gases, bars
c
c     Prepare the data. Nothing to prepare.
c
      ultag1 = 'T       '
      jt1 = len_trim(ultag1)
c
      if (qwrthd) then
c
c       Write an end marker for the previous table.
c
        write (ntabx,'(a,",EndTable:,",a,",")')
     $  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1
c
c       Write two empty lines before the table header.
c
        do i = 1,2
          write (ntabx,'(a,",")') ultag1(1:jt1)
        enddo
c
c       Write the table header.
c
        write (ntabx,1400) ultag1(1:jt1)
 1400   format(a,',Table T,Fugacities (bars),')
c
c       Write the column header.
c
        ulinex = 'Xi,t(days)'
        j1 = len_trim(ulinex)
        do n = 1,ngft
          j1 = j1 + 1
          ulinex(j1:j1) = ","
          ns = ngasfx(n)
          j2 = len_trim(uspec(ns)(1:24))
          jsum = j1 + j2
          ulinex(j1 + 1:jsum) = uspec(ns)(1:j2)
          j1 = jsum
        enddo
        if (qgftag) then
c
c         Have more gas species than can fit.
c
          j1 = j1 + 1
          ulinex(j1:j1) = ","
          ux16 = 'More omitted'
          j2 = len_trim(ux16)
          jsum = j1 + j2
          ulinex(j1 + 1:jsum) = ux16(1:j2)
          j1 = jsum
        endif
c
        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
      endif
c
c     Write the table data.
c
      ulinex = ''
      j1 = 0
c
      write (ux16,'(1pg11.4)') xi1
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') tdays
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      j1 = jsum
c
      do n = 1,ngft
        j1 = j1 + 1
        ulinex(j1:j1) = ","
        ns = ngasfx(n)
        ng = ns - ngrn1 + 1
        write (ux16,'(1pg11.4)') fugac(ng)
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
      enddo
c
      write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write Table W:
c
c       Xi
c       t(days)
c       Total mass of solids destroyed, grams
c       Total mass of solids created, grams
c       Net mass of solids created, grams
c       Total volume of solids destroyed, cc
c       Total volume of solids created, cc
c       Net volume of solids created, cc
c
      ultag1 = 'W       '
      jt1 = len_trim(ultag1)
c
      if (qwrthd) then
c
c       Write an end marker for the previous table.
c
        write (ntabx,'(a,",EndTable:,",a,",")')
     $  ultag1(1:jt1),ultags(1:jts)
        ultags = ultag1
        jts = jt1
c
c       Write two empty lines before the table header.
c
        do i = 1,2
          write (ntabx,'(a,",")') ultag1(1:jt1)
        enddo
c
c       Write the table header.
c
        write (ntabx,1410) ultag1(1:jt1)
 1410   format(a,',Table W,Overall mass and volume changes,')
c
c       Write the column header.
c
        ulinex = 'Xi,t(days),g destroyed,g created,g net,'
        j1 = len_trim(ulinex)
        ux48 = 'cc destroyed,cc created,cc net'
        j2 = len_trim(ux48)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux48(1:j2)
        j1 = jsum
c
        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
      endif
c
c     Write the table data.
c
      ulinex = ''
      j1 = 0
c
      write (ux16,'(1pg11.4)') xi1
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') tdays
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') wodrt
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') wosoct
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') dwoso
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') vodrt
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') vosoct
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      jsum = jsum + 1
      ulinex(jsum:jsum) = ','
      j1 = jsum
c
      write (ux16,'(1pg11.4)') dvoso
      call lejust(ux16)
      j2 = len_trim(ux16)
      jsum = j1 + j2
      ulinex(j1 + 1:jsum) = ux16(1:j2)
      j1 = jsum
c
      write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(4) .ge. 1) then
c
c       Write Table X:
c
c         Xi
c         t(days)
c         Mole fractions of solid solutions.
c
c       Note: cumulative averages are presented for the case of the
c       Fluid-Centered Flow-Through system.
c
c       Prepare the data.
c
        do n = 1,nxct
          xfrac(n) = 0.
        enddo
c
        qnewss = .false.
        if (iopt(1) .ne. 2) then
c
c         Closed or titration system.
c
          do np = nxrn1,nxrn2
            if (moph(np) .gt. 0.) then
              do n = 1,nxct
                if (uphase(np) .eq. ussphx(1,n)) then
c
c                 Have found a solid solution in the existing list.
c
                  nr1 = ncmpr(1,np)
                  nr2 = ncmpr(2,np)
                  do ns = nr1,nr2
                    if (uspec(ns) .eq. ussphx(2,n)) then
                      xfrac(n) = xbar(ns)
                      go to 410
                    endif
                  enddo
                endif
              enddo
c
c             Did not find a current solid solution in the existing
c             list. Add it to the list.
c
              qnewss = .true.
              nr1 = ncmpr(1,np)
              nr2 = ncmpr(2,np)
              do ns = nr1,nr2
                nxct = nxct + 1
                n = nxct
                ussphx(1,n) = uphase(np)
                ussphx(2,n) = uspec(ns)
                xfrac(n) = xbar(ns)
              enddo
            endif
  410       continue
          enddo
        else
c
c         Fluid-centered flow-through system.
c
          do np = nxrn1,nxrn2
            if (mopht(np) .gt. 0.) then
              do n = 1,nxct
                if (uphase(np) .eq. ussphx(1,n)) then
c
c                 Have found a solid solution in the existing list.
c
                  nr1 = ncmpr(1,np)
                  nr2 = ncmpr(2,np)
                  do ns = nr1,nr2
                    if (uspec(ns) .eq. ussphx(2,n)) then
                      xfrac(n) = mospt(ns)/mopht(np)
                      go to 420
                    endif
                  enddo
                endif
              enddo
c
c             Did not find a current solid solution in the existing
c             list. Add it to the list.
c
              qnewss = .true.
              nr1 = ncmpr(1,np)
              nr2 = ncmpr(2,np)
              do ns = nr1,nr2
                nxct = nxct + 1
                n = nxct
                ussphx(1,n) = uphase(np)
                ussphx(2,n) = uspec(ns)
                xfrac(n) = mospt(ns)/mopht(np)
              enddo
            endif
  420       continue
          enddo
        endif
c
        ultag1 = 'X       '
        jt1 = len_trim(ultag1)
c
        if (qwrthd) then
c
c         Write an end marker for the previous table.
c
          write (ntabx,'(a,",EndTable:,",a,",")')
     $    ultag1(1:jt1),ultags(1:jts)
          ultags = ultag1
          jts = jt1
c
c         Write two empty lines before the table header.
c
          do i = 1,2
            write (ntabx,'(a,",")') ultag1(1:jt1)
          enddo
c
c         Write the table header.
c
          write (ntabx,1440) ultag1(1:jt1)
 1440     format(a,',Table X,Solid solution mole fractions,')
        endif
c
        if (qnewss .or. kstep.le.0) then
c
c         Write the column header.
c
          ulinex = 'Xi,t(days)'
          j1 = len_trim(ulinex)
          do n = 1,nxct
            j1 = j1 + 1
            ulinex(j1:j1) = ","
            j2 = len_trim(ussphx(1,n))
            jsum = j1 + j2
            ulinex(j1 + 1:jsum) = ussphx(1,n)(1:j2)
            j1 = jsum
          enddo
c
          write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
        endif
c
c       Write the table data.
c
        ulinex = ''
        j1 = 0
c
        write (ux16,'(1pg11.4)') xi1
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        jsum = jsum + 1
        ulinex(jsum:jsum) = ','
        j1 = jsum
c
        write (ux16,'(1pg11.4)') tdays
        call lejust(ux16)
        j2 = len_trim(ux16)
        jsum = j1 + j2
        ulinex(j1 + 1:jsum) = ux16(1:j2)
        j1 = jsum
c
        do n = 1,nxct
          j1 = j1 + 1
          ulinex(j1:j1) = ","
          write (ux16,'(1pg11.4)') xfrac(n)
          call lejust(ux16)
          j2 = len_trim(ux16)
          jsum = j1 + j2
          ulinex(j1 + 1:jsum) = ux16(1:j2)
          j1 = jsum
        enddo
c
        write (ntabx,'(a,",",a,",")') ultag1(1:jt1),ulinex(1:j1)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write Table Z:
c
c       End of data marker
c
      if (qwrthd) then
c
c       Write an end marker for the previous table.
c
        write (ntabx,'("Z,EndTable:,",a,",")') ultags(1:jts)
c
        write (ntabx,'("Z,")')
        write (ntabx,'("Z,")')
c
c       Write the table header.
c
        write (ntabx,'("Z,Table Z,Endfile,")')
c
        write (ntabx,'("Z,Endtable:,Z,")')
c
      endif
c
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
