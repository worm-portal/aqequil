      subroutine rc3ocf(amu,jpfcmx,ifcphi1,ifcphi2,ifnnn,ifn2n,
     $ ifpsi1,ifpsi2,ifzeta,ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,
     $ ilpsi2,ilzeta,iodb,nmux,nmut,nmutmx,nodbmx,noutpt,nstmax,
     $ nttyo,uspec,zchar)
c
c     This subroutine recalculates the original coefficients for
c     the third-order Cphi, psi, and zeta parameters. Coefficients
c     for parameters originally in mu form are are unchanged.
c
c     After all the third-order data are calculated and otherwise
c     processed, the Cphi data are transformed into the equivalent
c     C data.
c
c     This routine supports the default option to us C, psi,
c     and zeta directly in the calculation of activity
c     coefficients. This can be changed to the former treatment
c     using the USEOLDPITZERMU option string in the input file
c     title (this changes qhawep from .true. to .false.).
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       amu    = array of coefficients for computing mu interaction
c                parameters
c       nmux   = array containing the indices of the species
c                  corresponding to data in the amu array
c       uspec  = array of species names
c       zchar  = array of charge numbers
c
c     Principal output:
c
c       amu    = modified, array of coefficients for computing
c                Cphi, psi, zeta, and some remaining mu parameters
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt,nttyo
c
      integer jpfcmx,nmutmx,nodbmx,nstmax
c
      integer ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,
     $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta
c
      integer nmut
c
      integer nmux(3,nmutmx),iodb(nodbmx)
c
      character(len=48) uspec(nstmax)
c
      real(8) amu(jpfcmx,nmutmx),zchar(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ilnobl
c
      integer iz1,iz3,inote,na,j,j1,j2,j3,n,nmu,nmu1,nmu2,nnote,
     $ ns1,ns2,ns3
c
      real(8) az1,az2,az3,cx,cx1,cx2,z1,z2,z3
c
c-----------------------------------------------------------------------
c
c     Recalculate the original coefficients for the Cphi, psi, and
c     zeta interaction parameters from the derived set of coefficients
c     for conventional mu interaction parameters.
c
c     Take advantage of the expected structure of the data for the
c     mu parameters. These occur in the following order:
c
c       cca data (obtained from Cphi data)
c       aac data (obtained from the same set of Cphi data)
c       cc'a data (obtained from the psi(cc'a) data)
c       aa'c data (obtained from the psi(aa'c) data)
c       nnn data (mu data, stays the same)
c       nnn' data (mu data, stays the same)
c       nca data (obtained from zeta data)
c
c     The species within a triplet follow the above pattern.
c     Thus, so see if one species appears twice, it is only
c     necessary to see if the first and second species are
c     identical. Thus, cca data will always be stored as
c     "cca", never as "cac" or "acc".
c
c     The strategy here is as follows: store all converted data
c     in the original amu array, as originally indexed. The nmxi
c     "pointer" array can retain its structure as well.
c
c     The first step is to convert the mu(cca) back to Cphi(ca)
c     data. The mu(aac) are likewise converted back to Cphi(ac) data.
c     This is the same set of Cphi data, just available with different
c     indexing (Cphi(ac) = Cphi(ca)). This duplication allows the amu
c     array to be used without changing its structure.
c
c     The mu(cc'a) data are converted back to psi(cc'a) data,
c     and the mu(aa'c) data are likewise converted back to
c     psi(aa'c) data. The two psi data sets here are separate. Note
c     that the conversion back to psi data uses the recalculated
c     Cphi data. Each psi datum has two "cognate" Cphi data.
c     These cognate Cphi data are for the constituent binary
c     pairs (e.g., the cognates for psi(cc'a) are Cphi(ca) and
c     Cphi(c'a). If the data for a cognate Cphi are zero, the
c     corresponding mu data may not be loaded for use during the
c     run. In that case, the data will not be found among the
c     recalculated Cphi data. Zero values will then be assumed.
c
c     The mu(nnn) and mu(nnn') data were originally in the same form
c     and remain that way.
c
c     The mu(nca) data are converted back to zeta(nca) data. This is
c     a simple conversion, as the potential cognate parameters
c     (mu(ncc) and mu(naa)) are not used in the standard Pitzer model
c     incorporated into EQ3/6 and other geochemical modeling codes.
c
      inote = 0
      nnote = 0
c
c     Make the mu(cca) to Cphi(ca) conversions.
c
      ifcphi1 = 1
      ilcphi1 = 0
c
      do nmu = ifcphi1,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        z1 = zchar(ns1)
        z3 = zchar(ns3)
        az1 = abs(z1)
        az3 = abs(z3)
c
        if ((ns1.eq.ns2 .and. ns2.ne.ns3) .and.
     $    (z1.gt.0 .and. z3.lt.0.)) then
c
c         Have an instance of mu(cca) data. Convert the coefficients
c         for mu(cca) to those for Cphi(ca).
c
          cx = 6.*sqrt(az3/az1)
          do j = 1,jpfcmx
            amu(j,nmu) = cx*amu(j,nmu)
          enddo
c
          ilcphi1 = nmu
c
        else
c
c         Have finished finding and processing all mu(cca) data.
c
          go to 110
c
        endif
      enddo
  110 continue
c
c     Make the mu(aac) to Cphi(ac) conversions.
c
      ifcphi2 = ilcphi1 + 1
      ilcphi2 = ilcphi1
c
      do nmu = ifcphi2,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        z1 = zchar(ns1)
        z3 = zchar(ns3)
        az1 = abs(z1)
        az3 = abs(z3)
c
        if ((ns1.eq.ns2 .and. ns2.ne.ns3) .and.
     $    (z1.lt.0 .and. z3.gt.0.)) then
c
c         Have an instance of mu(aac) data. Convert the coefficients
c         for mu(aac) to those for Cphi(ac).
c
          cx = 6.*sqrt(az3/az1)
          do j = 1,jpfcmx
            amu(j,nmu) = cx*amu(j,nmu)
          enddo
c
          ilcphi2 = nmu
c
        else
c
c         Have finished finding and processing all mu(aac) data.
c
          go to 120
c
        endif
      enddo
  120 continue
c
c     Make the mu(cc'a) to psi(cc'a) conversions.
c
      ifpsi1 = ilcphi2 + 1
      ilpsi1 = ilcphi2
c
      do nmu = ifpsi1,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        z1 = zchar(ns1)
        z2 = zchar(ns2)
        z3 = zchar(ns3)
        az1 = abs(z1)
        az2 = abs(z2)
        az3 = abs(z3)
c
        if ((ns1.ne.ns2 .and. ns2.ne.ns3) .and.
     $    (z1.gt.0. .and. z2.gt.0. .and. z3.lt.0.)) then
c
c         Have an instance of mu(cc'a) data. Convert the coefficients
c         for mu(cc'a) to those for psi(cc'a).
c
c         First find the indices of the cognate mu(cca) and
c         mu(c'c'a) data.
c
          inote = 0
c
          nmu1 = 0
          do n = ifcphi1,ilcphi1
            if (nmux(1,n).eq.ns1 .and. nmux(3,n).eq.ns3) then
              nmu1 = n
              go to 130
            endif
          enddo
c
          if (iodb(1) .gt. 0) then
            j1 = ilnobl(uspec(ns1)(1:24))
            j2 = ilnobl(uspec(ns2)(1:24))
            j3 = ilnobl(uspec(ns3)(1:24))
            write (noutpt,1000) uspec(ns1)(1:j1),uspec(ns3)(1:j3),
     $      uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
            write (nttyo,1000) uspec(ns1)(1:j1),uspec(ns3)(1:j3),
     $      uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
 1000       format(/' * Note - (EQLIBG/rc3ocf) Could not find the',
     $      ' Cphi data for',/7x,a,', ',a,'. They are needed to',
     $      ' recalculate',/7x,'the psi data for ',a,', ',a,', ',a,'.',
     $      /7x,'Zero values will be assumed.')
            inote = inote + 1
            nnote = nnote + 1
          endif
c
  130     continue
c
          nmu2 = 0
          do n = ifcphi1,ilcphi1
            if (nmux(1,n).eq.ns2 .and. nmux(3,n).eq.ns3) then
              nmu2 = n
              go to 140
            endif
          enddo
c
          if (iodb(1) .gt. 0) then
            j1 = ilnobl(uspec(ns1)(1:24))
            j2 = ilnobl(uspec(ns2)(1:24))
            j3 = ilnobl(uspec(ns3)(1:24))
            write (noutpt,1000) uspec(ns2)(1:j1),uspec(ns3)(1:j3),
     $      uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
            write (nttyo,1000) uspec(ns2)(1:j1),uspec(ns3)(1:j3),
     $      uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
            inote = inote + 1
            nnote = nnote + 1
          endif
c
  140     continue
c
          do j = 1,jpfcmx
            amu(j,nmu) = 6.*amu(j,nmu)
          enddo
          if (nmu1 .gt. 0) then
            cx1 = 0.5*az2/sqrt(az1*az3)
            do j = 1,jpfcmx
              amu(j,nmu) = amu(j,nmu) - cx1*amu(j,nmu1)
            enddo
          endif
          if (nmu2 .gt. 0) then
            cx2 = 0.5*az1/sqrt(az2*az3)
            do j = 1,jpfcmx
              amu(j,nmu) = amu(j,nmu) - cx2*amu(j,nmu2)
            enddo
          endif
c
          ilpsi1 = nmu
c
        else
c
c         Have finished finding and processing all mu(aac) data.
c
          go to 150
c
        endif
      enddo
  150 continue
c
c     Make the mu(aa'c) to psi(aa'c) conversions.
c
      ifpsi2 = ilpsi1 + 1
      ilpsi2 = ilpsi1
c
      do nmu = ifpsi2,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        z1 = zchar(ns1)
        z2 = zchar(ns2)
        z3 = zchar(ns3)
        az1 = abs(z1)
        az2 = abs(z2)
        az3 = abs(z3)
c
        if ((ns1.ne.ns2 .and. ns2.ne.ns3) .and.
     $    (z1.lt.0. .and. z2.lt.0. .and. z3.gt.0.)) then
c
c         Have an instance of mu(aa'c) data. Convert the coefficients
c         for mu(aa'c) to those for psi(aa'c).
c
c         First find the indices of the cognate mu(aac) and
c         mu(a'a'c) data.
c
          inote = 0
c
          nmu1 = 0
          do n = ifcphi2,ilcphi2
            if (nmux(1,n).eq.ns1 .and. nmux(3,n).eq.ns3) then
              nmu1 = n
              go to 160
            endif
          enddo
c
          if (iodb(1) .gt. 0) then
            j1 = ilnobl(uspec(ns1)(1:24))
            j2 = ilnobl(uspec(ns2)(1:24))
            j3 = ilnobl(uspec(ns3)(1:24))
            write (noutpt,1000) uspec(ns1)(1:j1),uspec(ns3)(1:j3),
     $      uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
            write (nttyo,1000) uspec(ns1)(1:j1),uspec(ns3)(1:j3),
     $      uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
            inote = inote + 1
            nnote = nnote + 1
          endif
c
  160     continue
c
          nmu2 = 0
          do n = ifcphi2,ilcphi2
            if (nmux(1,n).eq.ns2 .and. nmux(3,n).eq.ns3) then
              nmu2 = n
              go to 170
            endif
          enddo
c
          if (iodb(1) .gt. 0) then
            j1 = ilnobl(uspec(ns1)(1:24))
            j2 = ilnobl(uspec(ns2)(1:24))
            j3 = ilnobl(uspec(ns3)(1:24))
            write (noutpt,1000) uspec(ns2)(1:j1),uspec(ns3)(1:j3),
     $      uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
            write (nttyo,1000) uspec(ns2)(1:j1),uspec(ns3)(1:j3),
     $      uspec(ns1)(1:j1),uspec(ns2)(1:j2),uspec(ns3)(1:j3)
            inote = inote + 1
            nnote = nnote + 1
          endif
c
  170     continue
c
          do j = 1,jpfcmx
            amu(j,nmu) = 6.*amu(j,nmu)
          enddo
          if (nmu1 .gt. 0) then
            cx1 = 0.5*az2/sqrt(az3*az1)
            do j = 1,jpfcmx
              amu(j,nmu) = amu(j,nmu) - cx1*amu(j,nmu1)
            enddo
          endif
          if (nmu2 .gt. 0) then
            cx2 = 0.5*az1/sqrt(az3*az2)
            do j = 1,jpfcmx
              amu(j,nmu) = amu(j,nmu) - cx2*amu(j,nmu2)
            enddo
          endif
c
          ilpsi2 = nmu
c
        else
c
c         Have finished finding and processing all mu(aac) data.
c
          go to 180
        endif
      enddo
  180 continue
c
c     Skip past the mu(nnn) data.
c
      ifnnn = ilpsi2 + 1
      ilnnn = ilpsi2
c
      do nmu = ifnnn,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        z1 = zchar(ns1)
        iz1 = nint(z1)
c
        if ((ns1.eq.ns2 .and. ns2.eq.ns3) .and. iz1.eq.0) then
c
c         Have an instance of mu(nnn) data.
c
          ilnnn = nmu
c
        else
c
c         Have finished finding all mu(nnn) data.
c
          go to 190
c
        endif
      enddo
  190 continue
c
c     Skip past the mu(nnn') data.
c
      ifn2n = ilnnn + 1
      iln2n = ilnnn
c
      do nmu = ifn2n,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        z1 = zchar(ns1)
        z3 = zchar(ns3)
        iz1 = nint(z1)
        iz3 = nint(z3)
c
        if ((ns1.eq.ns2 .and. ns1.ne.ns3) .and.
     $    (iz1.eq.0 .and. iz3.eq.0)) then
c
c         Have an instance of mu(nnn') data.
c
          iln2n = nmu
c
        else
c
c         Have finished finding all mu(nnn') data.
c
          go to 200
c
        endif
      enddo
  200 continue
c
      ifzeta = iln2n + 1
      ilzeta = iln2n
c
c     Make the mu(nca) to zeta(nca) conversions.
c
      do nmu = ifzeta,nmut
        ns1 = nmux(1,nmu)
        ns2 = nmux(2,nmu)
        ns3 = nmux(3,nmu)
        z1 = zchar(ns1)
        z3 = zchar(ns3)
        az1 = abs(z1)
        az3 = abs(z3)
c
        if ((ns1.ne.ns2 .and. ns2.ne.ns3) .and.
     $    (iz1.eq.0 .and. z2.gt.0. .and. z3.lt.0.)) then
c
c         Have an instance of mu(nca) data. Convert the coefficients
c         for mu(nca) to those for zeta(nca).
c
          do j = 1,jpfcmx
            amu(j,nmu) = 6.*amu(j,nmu)
          enddo
c
          ilzeta = nmu
c
        else
c
c         Have finished finding and processing all mu(nca) data.
c
          go to 210
c
        endif
      enddo
  210 continue
c
      if (iodb(1) .gt. 0) then
c
c       Check results.
c
        write (noutpt,1100)
 1100   format(/1x,"Recalculated third-order Pitzer coefficients")
c
        write (noutpt,1110)
 1110   format(/1x,"Cphi(ca) data")
        do n = ifcphi1,ilcphi1
          ns1 = nmux(1,n)
          ns3 = nmux(3,n)
          j3 = ilnobl(uspec(ns3)(1:24))
          write (noutpt,1120) uspec(ns1),uspec(ns3)(1:j3)
 1120     format(/1x,a24,2x,a)
          write (noutpt,1130)
 1130     format(3x,"Cphi:")
          do j = 1,jpfcmx
            write (noutpt,1150) j,amu(j,n)
 1150       format (5x,"a",i1," = ",f10.6)
          enddo
        enddo
c
        write (noutpt,1200)
 1200   format(/1x,"Cphi(ac) data")
        do n = ifcphi2,ilcphi2
          ns1 = nmux(1,n)
          ns3 = nmux(3,n)
          j3 = ilnobl(uspec(ns3)(1:24))
          write (noutpt,1120) uspec(ns1),uspec(ns3)(1:j3)
          write (noutpt,1130)
          do j = 1,jpfcmx
            write (noutpt,1150) j,amu(j,n)
          enddo
        enddo
c
        write (noutpt,1310)
 1310   format(/1x,"psi(cc'a) data")
        do n = ifpsi1,ilpsi1
          ns1 = nmux(1,n)
          ns2 = nmux(2,n)
          ns3 = nmux(3,n)
          j3 = ilnobl(uspec(ns3)(1:24))
          write (noutpt,1320) uspec(ns1),uspec(ns2),uspec(ns3)(1:j3)
 1320     format(/1x,a24,2x,a24,2x,a)
          write (noutpt,1330)
 1330     format(3x,"psi:")
          do j = 1,jpfcmx
            write (noutpt,1150) j,amu(j,n)
          enddo
        enddo
c
        write (noutpt,1410)
 1410   format(/1x,"psi(aa'c) data")
        do n = ifpsi2,ilpsi2
          ns1 = nmux(1,n)
          ns2 = nmux(2,n)
          ns3 = nmux(3,n)
          j3 = ilnobl(uspec(ns3)(1:24))
          write (noutpt,1320) uspec(ns1),uspec(ns2),uspec(ns3)(1:j3)
          write (noutpt,1330)
          do j = 1,jpfcmx
            write (noutpt,1150) j,amu(j,n)
          enddo
        enddo
c
        write (noutpt,1510)
 1510   format(/1x,"mu(nnn) data")
        do n = ifnnn,ilnnn
          ns1 = nmux(1,n)
          ns2 = nmux(2,n)
          ns3 = nmux(3,n)
          j3 = ilnobl(uspec(ns3)(1:24))
          write (noutpt,1320) uspec(ns1),uspec(ns2),uspec(ns3)(1:j3)
          write (noutpt,1520)
 1520     format(3x,"mu:")
          do j = 1,jpfcmx
            write (noutpt,1150) j,amu(j,n)
          enddo
        enddo
c
        write (noutpt,1610)
 1610   format(/1x,"mu(nnn') data")
        do n = ifn2n,iln2n
          ns1 = nmux(1,n)
          ns2 = nmux(2,n)
          ns3 = nmux(3,n)
          j3 = ilnobl(uspec(ns3)(1:24))
          write (noutpt,1320) uspec(ns1),uspec(ns1),uspec(ns3)(1:j3)
          write (noutpt,1520)
          do j = 1,jpfcmx
            write (noutpt,1150) j,amu(j,n)
          enddo
        enddo
c
        write (noutpt,1710)
 1710   format(/1x,"zeta(nca) data")
        do n = ifzeta,ilzeta
          ns1 = nmux(1,n)
          ns2 = nmux(2,n)
          ns3 = nmux(3,n)
          j3 = ilnobl(uspec(ns3)(1:24))
          write (noutpt,1320) uspec(ns1),uspec(ns2),uspec(ns3)(1:j3)
          write (noutpt,1720)
 1720     format(3x,"zeta:")
          do j = 1,jpfcmx
            write (noutpt,1150) j,amu(j,n)
          enddo
        enddo
c
      endif
c
c     Transform the Cphi data into the equivalent C data.
c
      do n = ifcphi1,ilcphi1
        ns1 = nmux(1,n)
        ns3 = nmux(3,n)
        z1 = zchar(ns1)
        z3 = zchar(ns3)
        az1 = abs(z1)
        az3 = abs(z3)
        cx = 1./(2.*sqrt(az1*az3))
        do j = 1,jpfcmx
          amu(j,n) = cx*amu(j,n)
        enddo
      enddo
c
      do n = ifcphi2,ilcphi2
        ns1 = nmux(1,n)
        ns3 = nmux(3,n)
        z1 = zchar(ns1)
        z3 = zchar(ns3)
        az1 = abs(z1)
        az3 = abs(z3)
        cx = 1./(2.*sqrt(az1*az3))
        do j = 1,jpfcmx
          amu(j,n) = cx*amu(j,n)
        enddo
      enddo
c
  999 continue
      end
