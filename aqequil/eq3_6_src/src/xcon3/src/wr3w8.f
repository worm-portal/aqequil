      subroutine wr3w8(cgexj,cgexpi,covali,ehi,egexsi,fo2lgi,iebal3,
     $ ietmax,iodb,iopg,iopr,iopt,irdxc3,itdsf3,itermx,jgext,jetmax,
     $ jflgi,jgexti,jpres3,kxmod,mwtges,nbti,nbtmax,ncmpri,net,neti,
     $ netmax,newin,ngexti,ngexrt,nobswt,nodbmx,nopgmx,noprmx,noptmx,
     $ nsbswt,ntitl,ntitmx,nttyo,nxicmx,nxmdmx,nxmod,nxti,nxtimx,pei,
     $ press,qgexsh,rho,scamas,tdspkg,tdspl,tempc,tgexp,tolbt,toldl,
     $ tolspf,ucospi,uebal,ugexj,ugexji,ugexmo,ugexp,ugexpi,ugexr,
     $ ugexsi,umemi,uobsw,uredox,usbsw,usoli,uspeci,utitl,uhfgex,
     $ uvfgex,uxkgex,uxmod,xbari,xgexsi,xhfgex,xlkgex,xvfgex,
     $ xlkmod,zgexj)
c
c     This subroutine writes the EQ3NR input file in compact ("W")
c     format for version 8.0.
c
c     The calling sequence of this subroutine is identical to that of
c     XCON3/wr3d8.f.
c
c     The calling sequence of this subroutine is identical to that of
c     EQ3NR/rd3inw.f, EQ3NR/rd3ind.f, XCON3/rd3w8.f, and XCON3/rd3d8.f,
c     except that newin is added and ninpts, nprob, noutpt, qend, and
c     qrderr are deleted.
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
      integer ietmax,jetmax,nbtmax,netmax,nodbmx,nopgmx,noprmx,noptmx,
     $ ntitmx,nxmdmx,nxicmx,nxtimx
c
      integer newin,nttyo
c
      integer iodb(nodbmx),iopg(nopgmx),iopr(noprmx),iopt(noptmx),
     $ jgext(netmax),jgexti(netmax),jflgi(nbtmax),kxmod(nxmdmx),
     $ ncmpri(2,nxtimx),ngexrt(jetmax,netmax),ngexti(jetmax,netmax)
c
      integer iebal3,irdxc3,itdsf3,itermx,jpres3,nbti,net,neti,nobswt,
     $ nsbswt,ntitl,nxmod,nxti
c
      logical qgexsh
c
      character*80 utitl(ntitmx)
      character*56 ugexr(ietmax,jetmax,netmax)
      character*48 ucospi(nbtmax),uobsw(2,nbtmax),usbsw(2,nbtmax),
     $ uspeci(nbtmax),uxmod(nxmdmx)
      character*24 ugexmo(netmax),ugexp(netmax),ugexpi(netmax),
     $ ugexsi(ietmax,jetmax,netmax),umemi(nxicmx),usoli(nxtimx)
      character*24 uebal,uredox
      character*8 ugexj(jetmax,netmax),ugexji(jetmax,netmax),
     $ uhfgex(ietmax,jetmax,netmax),uvfgex(ietmax,jetmax,netmax),
     $ uxkgex(ietmax,jetmax,netmax)
c
      real*8 cgexj(jetmax,netmax),cgexpi(netmax),covali(nbtmax),
     $ egexsi(ietmax,jetmax,netmax),mwtges(netmax),tgexp(netmax),
     $ xbari(nxicmx),xhfgex(ietmax,jetmax,netmax),
     $ xlkgex(ietmax,jetmax,netmax),xvfgex(ietmax,jetmax,netmax),
     $ xlkmod(nxmdmx),zgexj(jetmax,netmax),xgexsi(ietmax,jetmax,netmax)
c
      real*8 ehi,fo2lgi,pei,press,rho,scamas,tdspkg,tdspl,tempc,tolbt,
     $ toldl,tolspf
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,iei,je,jei,j2,j3,j4,n,nbi,ne,nei,nr1,nr2,nssoti,nxi,
     $ nxic
c
      integer ilnobl
c
      logical qgexef
c
      character*8 uendit
c
c-----------------------------------------------------------------------
c
      data uendit /'endit.  '/
c
c-----------------------------------------------------------------------
c
c     The following is a bit of nonsense so compiler warnings will
c     not be generated saying that nttyo is not used.
c
      n = nttyo
      nttyo = n
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Title.
c
      do n = 1,ntitl
        j2 = ilnobl(utitl(n))
        write (newin,1020) utitl(n)(1:j2)
 1020   format(a)
      enddo
      j3 = ilnobl(uendit)
      if (ntitl .lt. ntitmx) write (newin,1020) uendit(1:j3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Special basis switches.
c
      write (newin,1100)
 1100 format('* Special basis switches')
c
      write (newin,1120) nsbswt
 1120 format(4x,'nsbswt= ',i3)
c
      do n = 1,nsbswt
        j2 = ilnobl(usbsw(1,n))
        write (newin,1140) usbsw(1,n)(1:j2)
 1140   format('species= ',a)
        j2 = ilnobl(usbsw(2,n))
        write (newin,1160) usbsw(2,n)(1:j2)
 1160   format(2x,'switch with= ',a)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Temperature.
c
      write (newin,1165)
 1165 format('* General')
c
      write (newin,1180) tempc
 1180 format(5x,'tempc= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Pressure.
c
      write (newin,1200) jpres3,press
 1200 format(4x,'jpres3= ',i3,/5x,'press= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Density.
c
      write (newin,1222) rho
 1222 format(7x,'rho= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Total dissolved solutes.
c
      write (newin,1227) itdsf3,tdspkg,tdspl
 1227 format(4x,'itdsf3= ',i3,/4x,'tdspkg= ',1pe12.5,5x,'tdspl= ',e12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Species to adjust for electrical balancing.
c
      j2 = ilnobl(uebal)
      write (newin,1510) iebal3,uebal(1:j2)
 1510 format(4x,'iebal3= 'i3,/5x,'uebal= ',a)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Default redox parameters.
c
      write (newin,1234) irdxc3
 1234 format(4x,'irdxc3= ',i3)
      write (newin,1242) fo2lgi,ehi
 1242 format(4x,'fo2lgi= ',1pe12.5,7x,'ehi= ',e12.5)
      j2 = ilnobl(uredox)
      write (newin,1252) pei,uredox(1:j2)
 1252 format(7x,'pei= ',1pe12.5,4x,'uredox= ',a)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Aqueous basis species and associated constraints.
c
      write (newin,2490)
 2490 format('* Aqueous basis species')
c
      do nbi = 1,nbti
        j2 = ilnobl(uspeci(nbi))
        write (newin,2530) uspeci(nbi)(1:j2)
 2530   format('species= ',a)
        write (newin,2550) jflgi(nbi),covali(nbi)
 2550   format(3x,'jflgi= ',i2,4x,'covali= ',1pg12.5)
        if (jflgi(nbi).eq.17 .or. jflgi(nbi).eq.18 .or.
     $  jflgi(nbi).eq.25) then
          j2 = ilnobl(ucospi(nbi))
          write (newin,2580) ucospi(nbi)(1:j2)
 2580     format(2x,'ucospi= ',a)
        endif
      enddo
      j3 = ilnobl(uendit)
      write (newin,2510) uendit(1:j3)
 2510 format(a)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ion exchanger creation.
c
      write (newin,1590)
 1590 format('* Ion exchangers')
c
      write (newin,1600) qgexsh
 1600 format(4x,'qgexsh= ',l8)
c
      write (newin,1610) net
 1610 format(7x,'net= ',i3)
c
      do ne = 1,net
        j2 = ilnobl(ugexp(ne))
        write (newin,1630) ugexp(ne)(1:j2)
 1630   format(5x,'ugexp= ',a)
c
        write (newin,1660) mwtges(ne)
 1660   format(4x,'mwtges= ',1pe12.5)
c
        j3 = ilnobl(ugexmo(ne))
        write (newin,1670) ugexmo(ne)(1:j3)
 1670   format(4x,'ugexmo= ',a)
c
        write (newin,1680) tgexp(ne)
 1680   format(5x,'tgexp= ',1pe12.5)
c
        write (newin,1700) jgext(ne)
 1700   format(5x,'jgext= ',i3)
c
        do je = 1,jgext(ne)
          j3 = ilnobl(ugexj(je,ne))
          write (newin,1740) ugexj(je,ne)(1:j3)
 1740     format(5x,'ugexj= ',a)
c
          write (newin,1760) cgexj(je,ne),zgexj(je,ne)
 1760     format(5x,'cgexj= '1pe12.5,5x,'zgexj= ',e12.5)
c
          write (newin,1810) ngexrt(je,ne)
 1810     format(4x,'ngexrt= ',i3)
c
          do n = 1,ngexrt(je,ne)
            j2 = ilnobl(ugexr(n,je,ne))
            write (newin,1840) ugexr(n,je,ne)(1:j2)
 1840       format(5x,'ugexr= ',a)
            j4 = ilnobl(uxkgex(n,je,ne))
            write (newin,1860) xlkgex(n,je,ne),uxkgex(n,je,ne)(1:j4)
 1860       format(4x,'xlkgex= ',1pe12.5,5x,'units= ',a)
            j4 = ilnobl(uhfgex(n,je,ne))
            write (newin,1870) xhfgex(n,je,ne),uhfgex(n,je,ne)(1:j4)
 1870       format(4x,'xhfgex= ',1pe12.5,5x,'units= ',a)
            j4 = ilnobl(uvfgex(n,je,ne))
            write (newin,1880) xvfgex(n,je,ne),uvfgex(n,je,ne)(1:j4)
 1880       format(4x,'xvfgex= ',1pe12.5,5x,'units= ',a)
          enddo
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Specified generic ion exchanger compositions.
c
      write (newin,1900)
 1900 format('* Ion exchanger compositions')
c
      write (newin,1920) neti
 1920 format(6x,'neti= ',i3)
      if (neti .le. 0) go to 250
c
      do nei = 1,neti
c
c       Find the corresponding ne index.
c
        qgexef = .true.
        do ne = 1,net
          j2 = ilnobl(ugexp(ne))
          j3 = ilnobl(ugexpi(nei))
          if (ugexp(nei)(1:j3) .eq. ugexp(ne)(1:j2)) then
            j2 = ilnobl(ugexmo(ne))
            if (ugexmo(ne)(1:j2) .eq. 'Gapon' .or.
     $        ugexmo(ne)(1:6) .eq. 'Gapon-' .or.
     $        ugexmo(ne)(1:j2) .eq. 'Vanselow' .or.
     $        ugexmo(ne)(1:9) .eq. 'Vanselow-') then
c
c             Input composition is described in terms of equivalent
c             fractions on the sites.
c
              qgexef = .true.
            else
c
c             Input composition is described in terms of mole fractions
c             on the sites.
c
              qgexef = .false.
            endif
            go to 240
          endif
        enddo
  240   continue
c
        j2 = ilnobl(ugexpi(nei))
        write (newin,1950) ugexpi(nei)(1:j2)
 1950   format(4x,'ugexpi= ',a)
c
        write (newin,1970) cgexpi(nei)
 1970   format(4x,'cgexpi= ',1pe12.5)
c
        write (newin,1990) jgexti(nei)
 1990   format(4x,'jgexti= ',i3)
c
        do jei = 1,jgexti(nei)
          j3 = ilnobl(ugexji(jei,nei))
          write (newin,2020) ugexji(jei,nei)(1:j3)
 2020     format(4x,'ugexji= ',a)
c
          write (newin,2030) ngexti(jei,nei)
 2030     format(4x,'ngexti= ',i3)
c
          do iei = 1,ngexti(jei,nei)
            if (qgexef) then
              write (newin,2060) ugexsi(iei,jei,nei),
     $        egexsi(iei,jei,nei)
 2060         format(4x,'ugexsi= ',a24,2x,'egexsi= ',1pe12.5)
            else
              write (newin,2070) ugexsi(iei,jei,nei),
     $        xgexsi(iei,jei,nei)
 2070         format(4x,'ugexsi= ',a24,2x,'xgexsi= ',1pe12.5)
            endif
          enddo
        enddo
      enddo
c
  250 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (newin,2100)
 2100 format('* Solid solution compositions')
c
      write (newin,2104) nxti
 2104 format(6x,'nxti= ',i3)
c
      do nxi = 1,nxti
        j2 = ilnobl(usoli(nxi))
        write (newin,2130) usoli(nxi)(1:j2)
 2130   format(5x,'usoli= ',a)
c
        nr1 = ncmpri(1,nxi)
        nr2 = ncmpri(2,nxi)
        nssoti = nr2 - nr1 + 1
c
        write (newin,2142) nssoti
 2142   format(4x,'nssoti= ',i3)
c
        do nxic = nr1,nr2
          write (newin,2170) umemi(nxic),xbari(nxic)
 2170     format(5x,'umemi= ',a24,5x,'xbari= ',1pe12.5)
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Number of nxmod options.
c
      write (newin,2330)
 2330 format('* Alter/suppress options')
c
      write (newin,2340) nxmod
 2340 format(5x,'nxmod= ',i3)
c
c     Nxmod options.
c
      do n = 1,nxmod
        j2 = ilnobl(uxmod(n))
        write (newin,2380) uxmod(n)(1:j2),kxmod(n),xlkmod(n)
 2380   format(3x,'species= ',a,
     $  /4x,'option= ',i2,14x,'xlkmod= ',1pe12.5)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopt model option switches.
c     Note: iopt(1) = iopt1, etc.
c
      write (newin,1365)
 1365 format('* Iopt, iopg, iopr, and iodb options')
c
      write (newin,1370)
 1370 format('*',15x,'1    2    3    4    5    6    7    8    9   10')
c
      write (newin,1390) (iopt(i), i = 1,20)
 1390 format(2x,'iopt1-10= ',10i5,/1x,'iopt11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopg activity coefficient option switches.
c     Note: iopg(1) = iopg1, etc.
c
      write (newin,1420) (iopg(i), i = 1,20)
 1420 format(2x,'iopg1-10= ',10i5,/1x,'iopg11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopr print option switches.
c     Note: iopr(1) = iopr1, etc.
c
      write (newin,1440) (iopr(i), i = 1,20)
 1440 format(2x,'iopr1-10= ',10i5,/1x,'iopr11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iodb debugging print option switches.
c     Note: iodb(1) = iodb1, etc.
c
      write (newin,1480) (iodb(i), i = 1,20)
 1480 format(2x,'iodb1-10= ',10i5,/1x,'iodb11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Numerical parameters. Convergence tolerances (tolbt and toldl)
c     and maximum number of Newton-Raphson iterations (itermx).
c
      write (newin,1270)
 1270 format('* Numerical parameters')
c
      write (newin,1290) tolbt,toldl
 1290 format (5x,'tolbt= ',1pe12.5,5x,'toldl= ',e12.5)
c
      write (newin,1310) itermx
 1310 format (4x,'itermx= ',i3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ordinary basis switches.
c
      write (newin,2700)
 2700 format('* Ordinary basis switches')
c
      write (newin,2710) nobswt
 2710 format(4x,'nobswt= ',i3)
c
      do n = 1,nobswt
        j2 = ilnobl(uobsw(1,n))
        write (newin,1140) uobsw(1,n)(1:j2)
        j2 = ilnobl(uobsw(2,n))
        write (newin,1160) uobsw(2,n)(1:j2)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Saturation flag tolerance (tolspf).
c
      write (newin,1292)
 1292 format('* Saturation flag tolerance')
c
      write (newin,1296) tolspf
 1296 format (4x,'tolspf= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Scale factor for the mass of aqueous solution to write
c     on the PICKUP file.
c
      write (newin,1315)
 1315 format('* Aqueous phase scale factor')
c
      write (newin,1330) scamas
 1330 format (4x,'scamas= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
