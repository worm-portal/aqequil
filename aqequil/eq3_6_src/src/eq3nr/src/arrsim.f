      subroutine arrsim(aamatr,acflg,actlg,bbig,cdrs,cjbasp,cnufac,
     $ conc,conclg,coval,delvec,dlogxw,eh,ehfac,eps100,gmmatr,iction,
     $ iindx1,iodb,ipivot,irdxc3,ixbasp,jcsort,jflag,jjndex,kbt,ker,
     $ khydr,kkndex,kmax,kwater,narn1,narn2,nbasp,nbt,nbti,nbtmax,
     $ nbw,ncosp,ndecsp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,nhydr,
     $ nodbmx,no2gaq,noutpt,npass,nredox,nstmax,nttyo,omega,qawfix,
     $ rhsvec,ucospi,uspec,xbar,xbarlg,xbarw,xbrwlg,xlke,xlks,
     $ zchar,zvclg1)
c
c     This subroutine computes starting estimates of species
c     concentrations that must be evaluated simultaneously. These
c     include cases of mean activity constraints, cases of equilibrium
c     constraints, and cases in which the log fO2 is constrained by Eh,
c     pe, or a redox couple.
c
c     This subroutine is called by:
c
c       EQ3NR/arrset.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
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
      integer kmax,nbtmax,ndrsmx,nodbmx,nstmax
c
      integer noutpt,nttyo
c
      integer iindx1(kmax),iodb(nodbmx),iction(nbtmax),ipivot(kmax),
     $ ixbasp(nbtmax),jcsort(nstmax),jjndex(nbtmax),jflag(nstmax),
     $ kkndex(nbtmax),nbasp(nbtmax),ncosp(nbtmax),ndecsp(nbtmax),
     $ ndrs(ndrsmx),ndrsr(2,nstmax)
c
      integer kbt,ker,khydr,kwater,narn1,narn2,nbt,nbti,nbw,nelect,
     $ nern1,nern2,nhydr,no2gaq,npass,nredox
c
      logical qawfix
c
      character(len=48) ucospi(nbtmax),uspec(nstmax)
c
      real(8) aamatr(kmax,kmax),acflg(nstmax),actlg(nstmax),
     $ cdrs(ndrsmx),cjbasp(nbtmax),cnufac(nstmax),conc(nstmax),
     $ conclg(nstmax),coval(nbtmax),delvec(kmax),dlogxw(nbtmax),
     $ gmmatr(kmax,kmax),rhsvec(kmax),xbar(nstmax),xbarlg(nstmax),
     $ xlks(nstmax),zchar(nstmax),zvclg1(kmax)
c
      real(8) bbig,eh,ehfac,eps100,omega,xbarw,xbrwlg,xlke
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ib,ibt,icol,ielect,ier,ihydr,io2gaq,irdxc3,irow,jfl,
     $ j2,j3,krow,n,nb,nbh,nb1,nb2,nerr,nr1,nr2,ns,nsc,ns1,ns2
c
      integer ilnobl,nbasis
c
      logical qpr,qx
c
      character(len=24) uaqeq,ublk24,ueh,umolfr,ured,ust1,ust2,ust3
      character(len=8) uptgas
c
      real(8) azns,azns1,azsum,cx,cxx,dzx,rx,xbrwlo,zp,zx,zxo
c
      real(8) coefdr,texp
c
c-----------------------------------------------------------------------
c
      data ublk24 /'                        '/
      data umolfr /'Mole fraction           '/
      data ueh    /'Eh                      '/
      data ured   /'Aqueous redox reaction  '/
      data uaqeq  /'Aqueous equilibrium     '/
      data uptgas /'Gas     '/
c
      data qpr    /.false./
c
c-----------------------------------------------------------------------
c
c     Save the current value of the mole fraction of water.
c
      xbrwlo = xbrwlg
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Build a structure for a matrix equation. Note that water
c     and its mole fraction relation may or may not be included
c     in this structure. Consequently, the kkndex array is redefined
c     from that built by EQ3NR/dawfix.f.
c
      ihydr = 0
      ielect = 0
      io2gaq = 0
c
      ib = 0
      do krow = 1,kbt
        nb = iindx1(krow)
        kkndex(nb) = 0
        ns = nbasp(nb)
        jfl = jflag(ns)
        qx = .false.
        if (ns .eq. narn1) then
          qx = jfl.eq.0 .and. npass.gt.1 .and. bbig.le.0.1
          qx = qx .or. qawfix
        elseif (ns.eq.nelect .or. ns.eq.no2gaq) then
          qx = irdxc3 .ne. 0
        else
          qx = jfl.eq.17 .or. jfl.eq.18 .or. jfl.eq.21 .or.
     $    jfl.eq.25 .or. jfl.eq.27
        endif
        if (qx) then
          ib = ib + 1
          jjndex(ib) = nb
          kkndex(nb) = 1
          if (ns .eq. nhydr) ihydr = ib
          if (ns .eq. nelect) ielect = ib
          if (ns .eq. no2gaq) io2gaq = ib
        endif
      enddo
      ibt = ib
c
c     Quit if no species concentrations are to be solved simultaneously.
c
      if (ibt. le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the iction array, which supports the jflag = 17, 18, and 21
c     options.
c
      do nb = 1,nbt
        iction(nb) = 0
      enddo
c
      do irow = 1,ibt
        nb = jjndex(irow)
        ns = nbasp(nb)
        jfl = jflag(ns)
        ns1 = 0
        if (jfl.eq.17 .or. jfl.eq.18 .or. jfl.eq.21) then
          ns1 = ncosp(nb)
c
c         Calling sequence substitutions:
c           ns1 for ns
c
          nb1 = nbasis(nbasp,nbt,nbtmax,ns1)
          if (nb1 .eq. 0) then
            j3 = ilnobl(uspec(ns)(1:24))
            j2 = ilnobl(uspec(ns1)(1:24))
            write (noutpt,1000) uspec(ns1)(1:j2),uspec(ns)(1:j3)
            write (nttyo,1000) uspec(ns1)(1:j2),uspec(ns)(1:j3)
 1000       format(/' * Error - (EQ3NR/arrsim) The specified',
     $      ' counterion',/7x,a,' in the',/7x,'jflag = 17, 18, or',
     $      ' 21 option  for the basis species',/7x,a," isn't in",
     $      ' the basis set.')
            stop
          endif
          if (kkndex(nb1) .ge. 1) then
            do icol = 1,ibt
              nb2 = jjndex(icol)
              if (nb2 .eq. nb1) then
                iction(nb) = icol
                go to 110
              endif
            enddo
          endif
        endif
  110   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Optional print of the matrix structure.
c
      if (iodb(3) .ge. 3) then
        write (noutpt,1020)
 1020   format(/7x,'--- Structure for Simultaneous Evaluation of',
     $  ' Starting Estimates ---',/)
c
        do irow = 1,ibt
          nb1 = jjndex(irow)
          ns1 = nbasp(nb1)
          jfl = jflag(ns1)
          ust1 = uspec(ns1)
          ust3 = ublk24
c
          if (ns1 .eq. narn1) ust2 = umolfr
c
          if (ns1.eq.nelect .or. ns1.eq.no2gaq) then
            if (ncosp(nb1) .eq. 0) then
              if (irdxc3 .lt. 0) then
                ust2 = ueh
              else
                ust2 = ured
              endif
              go to 120
            endif
          endif
c
          if (jfl .eq. 17) then
            ns = ncosp(nb1)
            ust2 = 'Log activity combination'
            ust3 = uspec(ns)
          elseif (jfl .eq. 18) then
            ns = ncosp(nb1)
            ust2 = 'Mean log activity'
            ust3 = uspec(ns)
          elseif (jfl .eq. 21) then
            ns = ncosp(nb1)
            ust2 = 'pHCl'
            ust3 = uspec(ns)
          elseif (jfl .eq. 27) then
            ust2 = uaqeq
          else
            do nb = 1,nbti
              nb2 = ndecsp(nb)
              if (nb1 .eq. nb2) then
                ust2 = ucospi(nb)(1:24)
                ust3 = ucospi(nb)(25:48)
                go to 120
              endif
            enddo
          endif
c
  120     j2 = ilnobl(ust2)
          write (noutpt,1030) irow,ust1,ust2(1:j2)
 1030     format(2x,i5,2x,a24,2x,a)
          j3 = ilnobl(ust3)
          if (ust3(1:24) .ne. ublk24(1:24))
     $    write (noutpt,1040) ust3(1:j3)
 1040     format(37x,'Constraint species= ',a)
        enddo
c
        write (noutpt,1050)
 1050   format(//10x,'--- The zvclg1 and acflg Values ---',/)
c
        do krow = 1,kbt
          nb = iindx1(krow)
          ns = nbasp(nb)
          if (kkndex(nb) .le. 0) then
            write (noutpt,1060) krow,uspec(ns),zvclg1(krow),acflg(ns)
 1060       format(2x,i5,2x,a24,2x,f10.4,2x,f10.4)
          endif
        enddo
c
        write (noutpt,1070)
 1070   format(/1x)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Set up the cnufac array.
c
      do nsc = narn1,narn2
        cnufac(nsc) = 0.
      enddo
c
      nr1 = ndrsr(1,narn1)
      if (jflag(narn1) .eq. 30) cnufac(narn1) = xbarlg(narn1)/cdrs(nr1)
c
      do nsc = narn1 + 1,narn2
        nr1 = ndrsr(1,nsc)
        if (jflag(nsc) .eq. 30) cnufac(nsc) = conc(nsc)/cdrs(nr1)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Build the matrix.
c
      do icol = 1,ibt
        do irow = 1,ibt
          aamatr(irow,icol) = 0.
        enddo
      enddo
c
      do irow = 1,ibt
        nb = jjndex(irow)
        ns = nbasp(nb)
        jfl = jflag(ns)
c
        if (ns .eq. narn1) then
c
c         Mole fraction of water, based on a first-order expansion:
c           log x(w)(new) = log x(w)(old)
c             + sum over s' in irow (d log x(w)/d log m(s'))
c               * (log m(s')(new) - log m(s')(old))
c
c         Compute the dlogxw array (d log xw/d log ms').
c
          call gdlgxw(cdrs,cjbasp,cnufac,conc,dlogxw,eps100,
     $    ixbasp,jcsort,jflag,narn1,narn2,nbasp,nbt,nbtmax,nbw,
     $    ndrs,ndrsmx,ndrsr,nern1,nern2,noutpt,nstmax,nttyo,omega,
     $    xbar,xbarw)
c
          do icol = 1,ibt
            nb1 = jjndex(icol)
            aamatr(irow,icol) = -dlogxw(nb1)
          enddo
          aamatr(irow,irow) = 1.0
        elseif (jfl .eq. 17) then
          ns1 = ncosp(nb)
          azns = abs(zchar(ns))
          azns1 = abs(zchar(ns1))
          aamatr(irow,irow) = azns1
          icol = iction(nb)
          if (icol .gt. 0) then
            zp = zchar(ns)*zchar(ns1)
            if (zp .lt. 0) then
              aamatr(irow,icol) = azns
            else
              aamatr(irow,icol) = -azns
            endif
          endif
        elseif (jfl .eq. 18) then
          ns1 = ncosp(nb)
          azns = abs(zchar(ns))
          azns1 = abs(zchar(ns1))
          azsum = azns + azns1
          aamatr(irow,irow) = azns1/azsum
          icol = iction(nb)
          if (icol .gt. 0) aamatr(irow,icol) = azns/azsum
        elseif (jfl .eq. 21) then
          ns1 = ncosp(nb)
          azns = abs(zchar(ns))
          azns1 = abs(zchar(ns1))
          aamatr(irow,irow) = -azns1
          icol = iction(nb)
          if (icol .gt. 0) then
            zp = zchar(ns)*zchar(ns1)
            if (zp .lt. 0) then
              aamatr(irow,icol) = -azns
            else
              aamatr(irow,icol) = azns
            endif
          endif
        elseif (jfl .eq. 27) then
          do icol = 1,ibt
            nb1 = jjndex(icol)
            ns1 = nbasp(nb1)
c
c           Calling sequence substitutions:
c             ns1 for nse
c
            aamatr(irow,icol) =
     $      coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,ns,nstmax)
          enddo
c
c         Calling sequence substitutions:
c           ns for nse
c
          aamatr(irow,irow)
     $    = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns,ns,nstmax)
        elseif (jfl .eq. 25) then
          ns1 = ncosp(nb)
          do icol = 1,ibt
            nb2 = jjndex(icol)
            ns2 = nbasp(nb2)
c
c           Calling sequence substitutions:
c             ns2 for nse
c             ns1 for ns
c
            aamatr(irow,icol)
     $      = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns2,ns1,nstmax)
          enddo
        elseif (irdxc3 .lt. 0) then
          aamatr(irow,io2gaq) = 1.
          if (ihydr .gt. 0) aamatr(irow,ihydr) = 4.
        else
          do icol = 1,ibt
            nb1 = jjndex(icol)
            ns1 = nbasp(nb1)
            if (kkndex(nb1).gt.0) then
c
c             Calling sequence substitutions:
c               ns1 for nse
c               nredox for ns
c
              aamatr(irow,icol)
     $        = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,nredox,nstmax)
c
            endif
          enddo
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 3) then
        write (noutpt,1100)
 1100   format(/10x,'--- Matrix ---',/)
        do irow = 1,ibt
          write (noutpt,1110) (aamatr(irow,icol), icol = 1,ibt)
 1110     format(2x,10(f7.2,2x))
        enddo
        write (noutpt,1070)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Build the right-hand-side vector.
c
      do irow = 1,ibt
        nb = jjndex(irow)
        ns = nbasp(nb)
        jfl = jflag(ns)
c
        if (ns .eq. narn1) then
c
c         Mole fraction of water, based on a first-order expansion:
c           log x(w)(new) = log x(w)(old)
c             + sum over s' in irow (d log x(w)/d log m(s'))
c               * (log m(s')(new) - log m(s')(old))
c
          rx = xbrwlo
          do krow = 1,kbt
            if (krow .ne. kwater) then
              nb1 = iindx1(krow)
              if (kkndex(nb1) .ge. 1)
     $        rx = rx - dlogxw(nb1)*zvclg1(krow)
            endif
          enddo
          rhsvec(irow) = rx
        elseif (jfl .eq. 17) then
          ns1 = ncosp(nb)
          azns = abs(zchar(ns))
          azns1 = abs(zchar(ns1))
          rx = coval(nb) - azns1*acflg(ns)
          zp = zchar(ns)*zchar(ns1)
          if (zp .lt. 0.) then
            rx = rx - azns*acflg(ns1)
          else
            rx = rx + azns*acflg(ns1)
          endif
          icol = iction(nb)
          if (icol .eq. 0) then
            if (zp .lt. 0.) then
              rx = rx - azns*conclg(ns1)
            else
              rx = rx + azns*conclg(ns1)
            endif
          endif
          rhsvec(irow) = rx
        elseif (jfl .eq. 18) then
          ns1 = ncosp(nb)
          azns = abs(zchar(ns))
          azns1 = abs(zchar(ns1))
          azsum = azns + azns1
          rx = coval(nb) - ( azns1/azsum )*acflg(ns)
     $    - ( azns/azsum )*acflg(ns1)
          icol = iction(nb)
          if (icol .eq. 0) rx = rx - ( azns/azsum )*conclg(ns1)
          rhsvec(irow) = rx
        elseif (jfl .eq. 21) then
          ns1 = ncosp(nb)
          azns = abs(zchar(ns))
          azns1 = abs(zchar(ns1))
          rx = coval(nb) + azns1*acflg(ns)
          zp = zchar(ns)*zchar(ns1)
          if (zp .lt. 0.) then
            rx = rx + azns*acflg(ns1)
          else
            rx = rx - azns*acflg(ns1)
          endif
          icol = iction(nb)
          if (icol .eq. 0) then
            if (zp .lt. 0.) then
              rx = rx + azns*conclg(ns1)
            else
              rx = rx - azns*conclg(ns1)
            endif
          endif
          rhsvec(irow) = rx
        elseif (jfl .eq. 27) then
          if (iodb(3) .ge. 3) then
c
c           Calling sequence substitutions:
c             noutpt for nf
c
            call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns,nstmax,uspec)
          endif
          rx = xlks(ns)
          nr1 = ndrsr(1,ns)
          nr2 = ndrsr(2,ns)
          do n = nr1,nr2
            ns2 = ndrs(n)
            cx = cdrs(n)
            if (ns2 .eq. narn1) then
              rx = rx - cx*actlg(narn1)
            else
c
c             Calling sequence substitutions:
c               ns2 for ns
c
              nb2 = nbasis(nbasp,nbt,nbtmax,ns2)
              if (kkndex(nb2) .le. 0) then
                rx = rx - cx*actlg(ns2)
              else
                rx = rx - cx*acflg(ns2)
              endif
            endif
          enddo
          rhsvec(irow) = rx
        elseif (jfl .eq. 25) then
          ns1 = ncosp(nb)
          if (iodb(3) .ge. 3) then
c
c           Calling sequence substitutions:
c             noutpt for nf
c             ns1 for ns
c
            call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns1,nstmax,uspec)
          endif
          rx = xlks(ns1)
          nr1 = ndrsr(1,ns1)
          nr2 = ndrsr(2,ns1)
          do n = nr1,nr2
            ns2 = ndrs(n)
            cx = cdrs(n)
            if (ns2 .eq. ns1) then
              if (uspec(ns1)(25:32) .eq. uptgas(1:8)) then
                rx = rx - cx*coval(nb)
              else
                rx = rx - cx*actlg(ns1)
              endif
            elseif (ns2 .eq. narn1) then
              rx = rx - cx*actlg(narn1)
            else
c
c             Calling sequence substitutions:
c               ns2 for ns
c
              nb2 = nbasis(nbasp,nbt,nbtmax,ns2)
              if (kkndex(nb2) .le. 0) then
                rx = rx - cx*actlg(ns2)
              else
                rx = rx - cx*acflg(ns2)
              endif
            endif
          enddo
          rhsvec(irow) = rx
        elseif (irdxc3 .lt. 0) then
          rx = 4.*eh/ehfac
          rx = rx + xlke + 2.*actlg(narn1)
          rx = rx - 4.*acflg(nhydr)
          nbh = iindx1(khydr)
          if (kkndex(nbh) .le. 0) rx = rx - 4.*conclg(nhydr)
          rhsvec(irow) = rx
        elseif (irdxc3 .gt. 0) then
          rx = xlks(nredox)
          nr1 = ndrsr(1,nredox)
          nr2 = ndrsr(2,nredox)
          do n = nr1,nr2
            ns2 = ndrs(n)
            cx = cdrs(n)
            if (ns2 .eq. narn1) then
              rx = rx - cx*actlg(narn1)
            else
c
c             Calling sequence substitutions:
c               ns2 for ns
c
              nb2 = nbasis(nbasp,nbt,nbtmax,ns2)
              if (kkndex(nb2) .le. 0) then
                rx = rx - cx*actlg(ns2)
              else
                rx = rx - cx*acflg(ns2)
              endif
            endif
          enddo
          rhsvec(irow) = rx
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 3) then
        write (noutpt,1200)
 1200   format(/10x,'--- Right-Hand-Side Vector ---',/)
        do irow = 1,ibt
          nb = jjndex(irow)
          ns = nbasp(nb)
          write (noutpt,1210) irow,uspec(ns),rhsvec(irow)
 1210     format(1x,i5,2x,a24,2x,g12.5)
        enddo
        write (noutpt,1070)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Solve the matrix equation.
c
c     Calling sequence substitutions:
c       ibt for kdim
c
      call msolvr(aamatr,delvec,gmmatr,ier,ipivot,ibt,kmax,
     $ noutpt,nttyo,qpr,rhsvec)
c
      if (ier .gt. 0) then
        write (noutpt,1220)
        write (nttyo,1220)
 1220   format(/' * Error - (EQ3NR/arrsim) The matrix equation',
     $  /7x,'required by the speciation model appears to be',
     $  /7x,'singular. Check the solubility and other constraints',
     $  /7x,'on the input file for a violation of the mineralogic',
     $  /7x,'phase rule.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Apply change limits and range limits. Then load the results.
c
      irow = 0
      do krow = 1,kbt
        nb = iindx1(krow)
        if (kkndex(nb) .ge. 1) then
          irow = irow + 1
          zx = delvec(irow)
          zxo = zvclg1(krow)
          if (zxo .gt. -99999.) then
            dzx = zx - zxo
            if (ns1 .eq. narn1) then
              if (dzx .gt. 0.2) dzx = 0.2
              if (dzx .lt. -0.2) dzx = -0.2
            else
              if (dzx .gt. 2.0) dzx = 2.0
              if (dzx .lt. -2.0) dzx = -2.0
            endif
            zx = zxo + dzx
          endif
          if (ns1.eq.narn1 .and. zx.gt.0.) zx = 0.
          if (zx .gt. 2.0) zx = 2.0
          if (zx .lt. -100.) zx = -100.
          zvclg1(krow) = zx
          ns = nbasp(nb)
          conclg(ns) = zx
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 3) then
        write (noutpt,1050)
        do krow = 1,kbt
          nb = iindx1(krow)
          ns = nbasp(nb)
          if (kkndex(nb) .ge. 1)
     $    write (noutpt,1060) krow,uspec(ns),zvclg1(krow),acflg(ns)
        enddo
        write (noutpt,1070)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check possible effects of large concentrations resulting from
c     equilibrium constraints.
c
      i = ibt
      if (ielect.gt.0 .or. io2gaq.gt.0) i = i - 1
      if (i .le. 0) go to 999
      nerr = 0
c
      do krow = 1,kbt
        nb = iindx1(krow)
        ns = nbasp(nb)
        if (kkndex(nb).ge.1 .and. ns.ne.nelect .and. ns.ne.no2gaq) then
          cx = conclg(ns)
          cxx = texp(cx)
          if (cxx .gt. 20.) then
            ker = 1
            j2 = ilnobl(uspec(ns)(1:24))
            write (noutpt,1230) uspec(ns)(1:j2),cxx
            write (nttyo,1230) uspec(ns)(1:j2),cxx
 1230       format(/' * Note - (EQ3NR/arrsim) The species ',a,
     $      /7x,'appears to have a required molality near ',g12.5)
            if (npass.gt.1 .and. cxx.gt.100.) nerr = nerr + 1
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nerr .gt. 0) then
        write (noutpt,1240)
        write (nttyo,1240)
 1240   format(/' * Error - (EQ3NR/arrsim) Reconsider your choice ',
     $  /7x,'of input constraints.')
        ker = 2
        go to 999
      endif
c
      if (ker .eq. 1) then
        write (noutpt,1250)
        write (nttyo,1250)
 1250   format(/' * Warning - (EQ3NR/arrsim) This run may crash',
     $  /7x,'because of poor choice of input constraints.')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
