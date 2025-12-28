      subroutine wrtabx(actlg,afrc1,aft1,alk,cteaq,dvoso,dwoso,
     $ eh,fo2lg,iindx1,iktmax,iopt,ipndx1,kmax,km1,kmt,kstep,kx1,
     $ kxt,loph,ncmpr,modr,mopht,narn1,mosp,nct,nctmax,noptmx,
     $ nptmax,nrct,nrctmx,nstmax,ntabx,ntidmx,ntitl2,ntitld,ntitmx,
     $ nxtmax,pe,ph,ppmwe,prcinf,press,prminf,qbye,qmod,qriinf,
     $ tempc,time1,uelem,uphase,uplatm,ureac,uspec,usteq6,utitl2,
     $ utitld,uveeq6,vodrt,vosoct,wodrt,woh2o,wosoct,xbar,xi1)
c
c     This subroutine writes the scratch tab file tabx. The length of
c     any line should not exceed 129 characters.
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
      integer iktmax,kmax,nctmax,noptmx,nptmax,nrctmx,nstmax,
     $ ntidmx,ntitmx,nxtmax
c
      integer ntabx
c
      integer iindx1(kmax),iopt(noptmx),ipndx1(kmax),ncmpr(2,nptmax)
c
      integer km1,kmt,kstep,kx1,kxt,narn1,nct,nrct,ntitl2,ntitld
c
      logical qbye,qmod,qriinf
c
      character*80 utitl2(ntitmx),utitld(ntidmx)
      character*48 uspec(nstmax)
      character*24 uphase(nptmax),ureac(nrctmx)
      character*8 uelem(nctmax)
      character*8 uplatm,usteq6,uveeq6
c
      real*8 actlg(nstmax),afrc1(nrctmx),cteaq(nctmax),loph(nptmax),
     $ modr(nrctmx),mopht(nptmax),mosp(nstmax),ppmwe(nctmax),
     $ xbar(nstmax)
c
      real*8 aft1,alk,dvoso,dwoso,eh,fo2lg,wodrt,wosoct,pe,ph,
     $ prcinf,press,prminf,tempc,time1,vodrt,vosoct,woh2o,xi1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with global dimensioning.
c
      character(len=16), dimension(:,:), allocatable :: uprss
      character(len=8), dimension(:,:), allocatable :: uprmn
      character(len=8), dimension(:), allocatable :: uelac
c
      real(8), dimension(:), allocatable :: lcteaq,ppmaq,xfrac,zvclgx
c
c     Saved values of local array sizes.
c
      integer isv_kmax,isv_nctmax,isv_nxcmax
c
      SAVE isv_kmax,isv_nctmax,isv_nxcmax
c
      SAVE lcteaq,ppmaq,uelac,uprmn,uprss,xfrac,zvclgx
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nxcmax
c
      integer i,j2,j3,j4,k,kcol,ktot,ktotlm,n,nc,nctpr,nlim,nlim1,nlim2,
     $ np,np1,nrc,ns
c
      integer ilnobl
c
      logical qphasc,qstart,qtitl
c
      character*1 ua,ub,uc,ud,ue,uf,ug,uh,uj,uk,ul,um,un,uo,up,uq,ur,
     $ us,ut,uu,uv,ux,uy
c
      real*8 cx,dx1,dx2,dx3,lalk,mx,tdays,tlogd,wkgh2o,xilog
c
      real*8 tlg
c
c-----------------------------------------------------------------------
c
      data ua,ub,uc,ud,ue,uf,ug,uh,uj,uk,ul,um,un,uo,up,uq,ur,us,
     $ ut,uu,uv,ux,uy/'a','b','c','d','e','f','g','h','j',
     $ 'k','l','m','n','o','p','q','r','s','t','u','v','x','y'/
c
c-----------------------------------------------------------------------
c
      nxcmax = nxtmax*iktmax
c
c     Allocate or reallocate local work arrays as needed.
c
      if (.not.ALLOCATED(uprss)) then
c
c       Local work arrays are not allocated. Zero the saved
c       array size variables. Note that only one array is tested
c       to see if it is allocated. It is assumed that all local
c       work arrays are either allocated or not.
c
        isv_kmax = 0
        isv_nctmax = 0
        isv_nxcmax = 0
      else
c
c       Local work arrays are allocated. Check to see if any of the
c       array size variables have changed. If so, deallocate
c       the corresponding local work arrays and zero the corresponding
c       saved size variables.
c
        if (kmax .ne. isv_kmax) then
          DEALLOCATE(uprmn)
          DEALLOCATE(uprss)
          DEALLOCATE(zvclgx)
          isv_kmax = 0
        endif
c
        if (nctmax .ne. isv_nctmax) then
          DEALLOCATE(uelac)
          DEALLOCATE(lcteaq,ppmaq)
          isv_nctmax = 0
        endif
c
        if (nxcmax .ne. isv_nxcmax) then
          DEALLOCATE(xfrac)
          isv_nxcmax = 0
        endif
      endif
c
c     At this point, the saved array size values are zero if the
c     corresponding arrays need to be allocated.
c
      if (isv_kmax .eq. 0) then
        ALLOCATE(uprmn(3,kmax))
        ALLOCATE(uprss(2,kmax))
        ALLOCATE(zvclgx(kmax))
        isv_kmax = kmax
      endif
c
      if (isv_nctmax .eq. 0) then
        ALLOCATE(uelac(nctmax))
        ALLOCATE(lcteaq(nctmax),ppmaq(nctmax))
        isv_nctmax = nctmax
      endif
c
      if (isv_nxcmax .eq. 0) then
        ALLOCATE(xfrac(nxcmax))
        isv_nxcmax = nxcmax
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Zero the contents of the local work arrays.
c
      do k = 1,kmax
        uprmn(1,k) = ' '
        uprmn(2,k) = ' '
        uprmn(3,k) = ' '
        uprss(1,k) = ' '
        uprss(2,k) = ' '
      enddo
c
      do k = 1,kmax
        zvclgx(k) = -99999.
      enddo
c
      do n = 1,nctmax
        uelac(n) = ' '
      enddo
c
      do n = 1,nctmax
        lcteaq(n) = -99999.
        ppmaq(n) = 0.
      enddo
c
      do n = 1,nxcmax
        xfrac(n) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write output tables of some selected parameters. These are
c     identified by the letters a, b, c, etc. They are written on
c     the file tabx in scrambled order (i.e., the lines of any given
c     tables are interspersed with those of the other tables). Each
c     line is marked in column one with the letter of the table to
c     which it belongs. When the run is complete, EQLIBU/dscram.f
c     writes the lines in descrambled form on file tab. File tabx
c     is preserved. The maximum line length for the tabx file is 129
c     characters. The maximum line length for the tab file is 128
c     characters.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c      Write table a. Header data.
c
c               Run title
c               Code version identification
c               Data file title
c
c     Note: the original primary title (utitl1) has been rolled over
c     to the current secondary title (utitl2).
c
      qstart = kstep .eq. 0
      if (qstart) then
        do i = 1,3
          write (ntabx,1000) ua
 1000     format(a1,1x)
        enddo
c
        do n = 1,ntitl2
          j2 = ilnobl(utitl2(n))
          write (ntabx,1005) ua,utitl2(n)(1:j2)
 1005     format(a1,1x,a)
        enddo
c
        write (ntabx,1000) ua
        j2 = ilnobl(uveeq6)
        j3 = ilnobl(usteq6)
        j4 = ilnobl(uplatm)
        write (ntabx,1010) ua,uveeq6(1:j2),usteq6(1:j3),uplatm(1:j4)
 1010   format(a1,1x,'Running EQ3/6-V',a,'-EQ6-EXE-',a,'-',a)
        write (ntabx,1015) ua
 1015   format(a1)
        write (ntabx,1000) ua
c
        do n = 1,ntitld
          j2 = ilnobl(utitld(n))
          write (ntabx,1005) ua,utitld(n)(1:j2)
        enddo
      endif
c
      qtitl = qstart
      if (iopt(18) .ge. 1) qtitl = .false.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write table b. Major run parameters.
c
c               Overall reaction progress
c               Log of overall reaction progress
c               Time, days
c               Log days
c               Temperature, C
c               Pressure, bars
c               pH
c               log fO2
c               Eh, volts
c               pe
c               Mass of solvent, kg
c               Total affinity, kcal
c
      if (qtitl) then
        do i = 1,3
          write (ntabx,1000) ub
        enddo
c
        write (ntabx,1020) ub
 1020   format(a1,'       xi     log xi    time, d  log days',
     $  '    tempc  ','   press       ph     log fo2      eh',
     $  '        pe      kg h2o   tot aff')
        write (ntabx,1000) ub
      endif
c
cXX  Will eventually have to revise the tab file (-99999. vs. -999.).
cXX   xilog = -99999.
      xilog = -999.
      if (xi1 .gt. 0.) xilog = tlg(xi1)
      wkgh2o = 0.001*woh2o
      if (qriinf) then
        tdays = prcinf
cXX     tlogd = +99999.
        tlogd = +999.
      else
        tdays = time1/86400.
        if (time1 .gt. prminf) then
          tlogd = tlg(tdays)
        else
cXX       tlogd = -99999.
          tlogd = -999.
        endif
      endif
c
cXX
      if (eh .le. -99999.) eh = -999.
      if (pe .le. -99999.) pe = -999.
      if (aft1 .ge. 99999.) aft1 = 999.
cXX
c
      write (ntabx,1025) ub,xi1,xilog,tdays,tlogd,tempc,press,ph,fo2lg,
     $ eh,pe,wkgh2o,aft1
 1025 format(a1,1x,1pe10.3,0pf10.4,1pe10.3,9(0pf10.4))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write table c. Miscellaneous aqueous solution composition
c     parameters.
c
c               Log of activity of water
c               Alkalinity, eq/kg.H2O (not defined for
c                 temperatures greater than 50 C)
c
      if (qtitl) then
        do i = 1,3
          write (ntabx,1000) uc
        enddo
c
        write (ntabx,1030) uc
 1030   format(a1,1x,'   log xi   time, d   log days  log alk',
     $  '  log tot   log tot   log tot   log a h2o')
        write (ntabx,1035) uc
 1035   format(a1,41x,'  co3--     so4--     s--')
        write (ntabx,1000) uc
      endif
c
cXX   lalk = -99999.
      lalk = -999.
      if (alk .gt. 0.) lalk = tlg(alk)
      dx1 = 0.
      dx2 = 0.
      dx3 = 0.
      write (ntabx,1040) uc,xilog,tdays,tlogd,lalk,dx1,dx2,dx3,
     $ actlg(narn1)
 1040 format(a1,1x,0pf10.4,1pe10.3,10(0pf10.4))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write tables d, e, f, g, and h. Aqueous solution composition
c     in terms of molalities of chemical elements. Write
c     corresponding tables j, k, l, m, and n (aqueous solution
c     composition in terms of ppm of chemical elements).
c
      nctpr = 0
      do nc = 1,nct
        if (uelem(nc)(1:2).ne.'H ' .and. uelem(nc)(1:2) .ne. 'O ') then
          nctpr = nctpr + 1
          uelac(nctpr) = uelem(nc)
          cx = cteaq(nc)
          lcteaq(nctpr) = tlg(cx)
cXX  Will eventually have to revise the tab file (-99999. vs. -999.).
          if (lcteaq(nctpr) .lt. -999.) lcteaq(nctpr) = -999.
          ppmaq(nctpr) = ppmwe(nc)
        endif
      enddo
c
      if (nctpr .le. 0) go to 100
      nlim2 = min(9,nctpr)
c
      if (qtitl) then
        do i = 1,3
          write (ntabx,1000) ud
          write (ntabx,1000) uj
        enddo
c
        write (ntabx,1045) ud
 1045   format(a1,20x,'log molality of dissolved elements')
        write (ntabx,1050) uj
 1050   format(a1,17x,'ppm (mg/kg.sol) of dissolved elements')
        write (ntabx,1000) ud
        write (ntabx,1000) uj
        write (ntabx,1055) ud,(uelac(n), n = 1,nlim2)
        write (ntabx,1055) uj,(uelac(n), n = 1,nlim2)
 1055   format(a1,1x,'   log xi   time, d   log days ',10(4x,a3,3x))
        write (ntabx,1000) ud
        write (ntabx,1000) uj
      endif
c
      write (ntabx,1040) ud,xilog,tdays,tlogd,(lcteaq(n), n = 1,nlim2)
      write (ntabx,1065) uj,xilog,tdays,tlogd,(ppmaq(n), n = 1,nlim2)
 1065 format(a1,1x,0pf10.4,1pe10.3,0pf10.4,9(1x,g8.3,1x))
c
      if (nctpr .le. nlim2) go to 100
      nlim1 = nlim2 + 1
      nlim2 = min(18,nctpr)
c
      if (qtitl) then
        do i = 1,3
          write (ntabx,1000) ue
          write (ntabx,1000) uk
        enddo
c
        write (ntabx,1055) ue,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1055) uk,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1000) ue
        write (ntabx,1000) uk
      endif
c
      write (ntabx,1040) ue,xilog,tdays,tlogd,
     $ (lcteaq(n), n = nlim1,nlim2)
      write (ntabx,1065) uk,xilog,tdays,tlogd,
     $ (ppmaq(n), n = nlim1,nlim2)
c
      if (nctpr .le. nlim2) go to  100
      nlim1 = nlim2 + 1
      nlim2 = min(27,nctpr)
c
      if (qtitl) then
        do i = 1,3
          write (ntabx,1000) uf
          write (ntabx,1000) ul
        enddo
c
        write (ntabx,1055) uf,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1055) ul,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1000) uf
        write (ntabx,1000) ul
      endif
c
      write (ntabx,1040) uf,xilog,tdays,tlogd,
     $ (lcteaq(n), n = nlim1,nlim2)
      write (ntabx,1065) ul,xilog,tdays,tlogd,
     $ (ppmaq(n), n = nlim1,nlim2)
c
      if (nctpr .le. nlim2) go to  100
      nlim1 = nlim2 + 1
      nlim2 = min(36,nctpr)
c
      if (qtitl) then
        do i = 1,3
          write (ntabx,1000) ug
          write (ntabx,1000) um
        enddo
c
        write (ntabx,1055) ug,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1055) um,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1000) ug
        write (ntabx,1000) um
      endif
c
      write (ntabx,1040) ug,xilog,tdays,tlogd,
     $ (lcteaq(n), n = nlim1,nlim2)
      write (ntabx,1065) um,xilog,tdays,tlogd,
     $ (ppmaq(n), n = nlim1,nlim2)
c
      if (nctpr .le. nlim2) go to  100
      nlim1 = nlim2 + 1
      nlim2 = min(45,nctpr)
c
      if (qtitl) then
        do i = 1,3
          write (ntabx,1000) uh
          write (ntabx,1000) un
        enddo
c
        write (ntabx,1055) uh,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1055) un,(uelac(n), n = nlim1,nlim2)
        write (ntabx,1000) uh
        write (ntabx,1000) un
      endif
c
      write (ntabx,1040) uh,xilog,tdays,tlogd,
     $ (lcteaq(n), n = nlim1,nlim2)
      write (ntabx,1065) un,xilog,tdays,tlogd,
     $ (ppmaq(n), n = nlim1,nlim2)
c
  100 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write table o. Product solid solution compositions.
c
c                  Endmember compositions of product phases.
c
c     Write header on first pass through.
c
      if (qtitl) then
        do i = 1,3
          write (ntabx,1000) uo
        enddo
c
        write (ntabx,2100) uo
 2100 format(a1,20x,'solid solution product compositions')
        write (ntabx,1000) uo
      endif
c
c     First check if solid solution phases are present in matrix.
c
      ktot = kxt - kx1 + 1
c
      if (ktot .le. 0) go to 400
c
c     Go through solid solutions in matrix.
c
      n = 0
      do kcol = kx1,kxt
        ns = iindx1(kcol)
        np = ipndx1(kcol)
        n = n + 1
        uprss(1,n) = uspec(ns)(1:10)
        uprss(2,n) = uspec(ns)(11:20)
        xfrac(n) = xbar(ns)
      enddo
c
c     Write out solid solution name and composition.
c
      nlim = min(6,n)
      write (ntabx,2006) uo,(uprss(1,i),i = 1,nlim)
      write (ntabx,2011) uo,(uprss(2,i),i = 1,nlim)
 2006 format(a1,1x,'   log xi ',8x,6(1x,a10,1x))
 2011 format(a1,19x,6(1x,a10,1x))
      write (ntabx,2007) uo,xilog,(xfrac(i),i = 1,nlim)
 2007 format(a1,4x,f9.4,6x,6(1x,f10.8,1x))
      if (n .gt. nlim) then
        nlim = min(12,n)
        write (ntabx,2008) uo,(uprss(1,i),i = 7,nlim)
        write (ntabx,2008) uo,(uprss(2,i),i = 7,nlim)
        write (ntabx,2009) uo,(xfrac(i),i = 7,nlim)
 2008   format(/,a1,19x,5(1x,a10,1x),1x,a10)
 2009   format(a1,19x,6(1x,f10.8,1x))
      endif
      if (n .gt. nlim) then
        nlim = min(18,n)
        write (ntabx,2008) uo,(uprss(1,i),i = 13,nlim)
        write (ntabx,2008) uo,(uprss(2,i),i = 13,nlim)
        write (ntabx,2009) uo,(xfrac(i),i = 13,nlim)
      endif
      if (n .gt. nlim) then
        nlim = min(24,n)
        write (ntabx,2008) uo,(uprss(1,i),i = 19,nlim)
        write (ntabx,2008) uo,(uprss(2,i),i = 19,nlim)
        write (ntabx,2009) uo,(xfrac(i),i = 19,nlim)
      endif
      if (n .gt. 24) write (ntabx,2051) uo
 2051 format(a1,1x,'number of solid solution product phases > 24')
  400 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write tables p, q, r, and s. Product minerals.
c
c               Log of mass (moles) for closed system
c               Log of cumulative mass (moles) for flow-through
c                  system
c
      if (qtitl) then
        do i = 1,3
          write (ntabx,1000) up
        enddo
c
        if (iopt(1) .eq. 2) then
          write (ntabx,2000) up
 2000     format(a1,20x,'log of moles of product minerals (cumulative)')
        else
          write (ntabx,2002) up
 2002     format(a1,20x,'log of moles of product minerals')
        endif
        write (ntabx,1000) up
      endif
c
      ktot = kxt - km1 + 1
c
      if (ktot .le. 0) go to 200
c
c     ktotlm = no. of minerals that can be printed in the tables.
c
      ktotlm = 36
c
      if (kmt .ge. km1) then
        do kcol = km1,kmt
          n = kcol - km1 + 1
          if (n .gt. ktotlm) go to 125
          np = ipndx1(kcol)
          uprmn(1,n) = uphase(np)(1:8)
          uprmn(2,n) = uphase(np)(9:16)
          uprmn(3,n) = uphase(np)(17:24)
          if (iopt(1) .eq. 2) then
            zvclgx(n) = tlg(mopht(np))
          else
            zvclgx(n) = loph(np)
          endif
cXX  Will eventually have to revise the tab file (-99999. vs. -999.).
          if (zvclgx(n) .lt. -999.) zvclgx(n) = -999.
        enddo
      else
c
c     If all products present are solid solutions, reset counter to 0.
c
        n = 0
      endif
c
      if (kxt .ge. kx1) then
        np1 = 0
        do kcol = kx1,kxt
          ns = iindx1(kcol)
          np = ipndx1(kcol)
c
c         The values of kx1 - kxt correspond to all the endmembers of
c         solid solutions in the PRS. The following line checks to
c         see if the next solid solution phase has been reached.
c
          if (np .ne. np1) then
            np1 = np
            n = n + 1
            if (n .gt. ktotlm) go to 125
            uprmn(1,n) = uphase(np)(1:8)
            uprmn(2,n) = uphase(np)(9:16)
            uprmn(3,n) = uphase(np)(17:24)
            if (iopt(1) .eq. 2) then
              zvclgx(n) = tlg(mopht(np))
            else
              zvclgx(n) = loph(np)
            endif
            if (zvclgx(n) .lt. -999.) zvclgx(n) = -999.
          endif
        enddo
      endif
  125 ktot = n
c
      nlim2 = min(9,ktot)
c
      qphasc = qmod .or. qbye .or. qtitl
      if (qphasc) then
c
        do i = 1,3
          write (ntabx,1000) up
        enddo
c
        write (ntabx,2005) up,(uprmn(1,n), n = 1,nlim2)
 2005   format(a1,1x,'   log xi   time, d   log days ',9(1x,a8,1x))
        write (ntabx,2010) up,(uprmn(2,n), n = 1,nlim2)
 2010   format(a1,32x,9(1x,a8,1x),1x,a8)
        write (ntabx,2010) up,(uprmn(3,n), n = 1,nlim2)
        write (ntabx,1000) up
      endif
c
      write (ntabx,1040) up,xilog,tdays,tlogd,(zvclgx(n), n = 1,nlim2)
c
      if (ktot .le. nlim2) go to 200
      nlim1 = nlim2 + 1
      nlim2 = min(18,ktot)
c
      if (qphasc) then
        do i = 1,3
          write (ntabx,1000) uq
        enddo
c
        write (ntabx,2005) uq,(uprmn(1,n), n = nlim1,nlim2)
        write (ntabx,2010) uq,(uprmn(2,n), n = nlim1,nlim2)
        write (ntabx,2010) uq,(uprmn(3,n), n = nlim1,nlim2)
        write (ntabx,1000) uq
      endif
c
      write (ntabx,1040) uq,xilog,tdays,tlogd,
     $ (zvclgx(n), n = nlim1,nlim2)
c
      if (ktot .le. nlim2) go to 200
      nlim1 = nlim2 + 1
      nlim2 = min(27,ktot)
c
      if (qphasc) then
        do i = 1,3
          write (ntabx,1000) ur
        enddo
c
        write (ntabx,2005) ur,(uprmn(1,n), n = nlim1,nlim2)
        write (ntabx,2010) ur,(uprmn(2,n), n = nlim1,nlim2)
        write (ntabx,2010) ur,(uprmn(3,n), n = nlim1,nlim2)
        write (ntabx,1000) ur
      endif
c
      write (ntabx,1040) ur,xilog,tdays,tlogd,
     $ (zvclgx(n), n = nlim1,nlim2)
c
      if (ktot .le. nlim2) go to 200
      nlim1 = nlim2 + 1
      nlim2 = min(36,ktot)
c
      if (qphasc) then
        do i = 1,3
          write (ntabx,1000) us
        enddo
c
        write (ntabx,2005) us,(uprmn(1,n), n = nlim1,nlim2)
        write (ntabx,2010) us,(uprmn(2,n), n = nlim1,nlim2)
        write (ntabx,2010) us,(uprmn(3,n), n = nlim1,nlim2)
        write (ntabx,1000) us
      endif
c
      write (ntabx,1040) us,xilog,tdays,tlogd,
     $ (zvclgx(n), n = nlim1,nlim2)
c
  200 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write tables t and u. Reactants.
c
c               Log of mass (moles) of reactant destroyed.
c
      if (nrct .le. 0) go to 215
c
      do nrc = 1,nrct
        mx = abs(modr(nrc))
cXX     zvclgx(nrc) = -99999.
        zvclgx(nrc) = -999.
        if (mx .gt. 0.) zvclgx(nrc) = tlg(mx)
      enddo
c
      nlim2 = min(9,nrct)
c
      if (qtitl) then
        do i = 1,3
          write (ntabx,1000) ut
        enddo
c
        write (ntabx,2020) ut
 2020 format(a1,20x,'log of destroyed moles of reactants')
        write (ntabx,1000) ut
        write (ntabx,2005) ut,(ureac(nrc)(1:8), nrc = 1,nlim2)
        write (ntabx,2010) ut,(ureac(nrc)(9:16), nrc = 1,nlim2)
        write (ntabx,2010) ut,(ureac(nrc)(17:24), nrc = 1,nlim2)
        write (ntabx,1000) ut
      endif
c
      write (ntabx,1040) ut,xilog,tdays,tlogd,
     $ (zvclgx(nrc), nrc = 1,nlim2)
c
      if (nrct .le. nlim2) go to 215
      nlim1 = nlim2 + 1
      nlim2 = min(18,nrct)
c
      if (qtitl) then
        do i = 1,3
          write (ntabx,1000) uu
        enddo
c
        write (ntabx,2005) uu,(ureac(nrc)(1:8), nrc = nlim1,nlim2)
        write (ntabx,2010) uu,(ureac(nrc)(9:16), nrc = nlim1,nlim2)
        write (ntabx,2010) uu,(ureac(nrc)(17:24), nrc = nlim1,nlim2)
        write (ntabx,1000) uu
      endif
c
      write (ntabx,1040) uu,xilog,tdays,tlogd,
     $ (zvclgx(nrc), nrc = nlim1,nlim2)
c
  215 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write table v. Overall masses and volumes.
c
c               Total mass of solids destroyed, grams
c               Total mass of solids created, grams
c               Net mass of solids created, grams
c               Total volume of solids destroyed, cc
c               Total volume of solids created, cc
c               Net volume of solids created, cc
c
      if (qtitl) then
        do i = 1,3
          write (ntabx,1000) uv
        enddo
c
        write (ntabx,2030) uv
 2030   format(a1,1x,'      xi     log xi   time, d   log days   g des',
     $  '     g cre     g net     cc des    cc cre    cc net')
        write (ntabx,1000) uv
      endif
c
      write (ntabx,2035) uv,xi1,xilog,tdays,tlogd,wodrt,wosoct,dwoso,
     $ vodrt,vosoct,dvoso
 2035 format(a1,1x,1pe10.3,0pf10.4,1pe10.3,0pf10.4,8(1pe10.3))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write tables x and y. Affinities of irreversible reactions.
c
c               Affinities of individual reactants, kcal.
c
      if (nrct .le. 0) go to 300
      nlim2 = min(9,nrct)
c
      if (qtitl) then
        do i = 1,3
          write (ntabx,1000) ux
        enddo
c
        write (ntabx,2040) ux
 2040   format(a1,20x,'affinities of irreversible reactions')
        write (ntabx,1000) ux
        write (ntabx,2005) ux,(ureac(nrc)(1:8), nrc = 1,nlim2)
        write (ntabx,2010) ux,(ureac(nrc)(9:16), nrc = 1,nlim2)
        write (ntabx,2010) ux,(ureac(nrc)(17:24), nrc = 1,nlim2)
        write (ntabx,1000) ux
      endif
c
      write (ntabx,1040) ux,xilog,tdays,tlogd,
     $ (afrc1(nrc), nrc = 1,nlim2)
c
      if (nrct .le. nlim2) go to 300
      nlim1 = nlim2 + 1
      nlim2 = min(18,nrct)
      if (qtitl) then
        do i = 1,3
          write (ntabx,1000) uy
        enddo
c
        write (ntabx,2005) uy,(ureac(nrc)(1:8), nrc = nlim1,nlim2)
        write (ntabx,2010) uy,(ureac(nrc)(9:16), nrc = nlim1,nlim2)
        write (ntabx,2010) uy,(ureac(nrc)(17:24), nrc = nlim1,nlim2)
        write (ntabx,1000) uy
      endif
c
      write (ntabx,1040) uy,xilog,tdays,tlogd,
     $ (afrc1(nrc), nrc = nlim1,nlim2)
c
  300 continue
c
      end
