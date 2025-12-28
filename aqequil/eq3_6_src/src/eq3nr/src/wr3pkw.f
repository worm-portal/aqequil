      subroutine wr3pkw(electr,cgexj,ietmax,iopg,jetmax,jflgi,
     $ jgext,kbt,kct,kdim,kmax,kmt,kxmod,kxt,mtbi,mtbaqi,mwtges,
     $ nbti,nbtmax,net,netmax,newin,ngexrt,nobswt,nopgmx,nsbswt,
     $ ntitl2,ntitmx,nxmdmx,nxmod,qgexsh,pressi,tempci,tgexp,
     $ ubmtbi,ugexj,ugexmo,ugexp,ugexr,uhfgex,uobsw,usbsw,utitl2,
     $ uvfgex,uxkgex,uxmod,uzveci,xhfgex,xlkgex,xlkmod,xvfgex,
     $ zgexj,zvclgi)
c
c     This subroutine writes the pickup file in compact ("W") format.
c     This file is used to communicate data from EQ3NR to EQ6 (it
c     comprises the bottom half of an EQ6 input file).
c
c     This subroutine is a near clone of parts of XCON6/wr6w8.f and
c     EQ6/wr6pkw.f.
c
c     The calling sequence of this subroutine is identical to that of
c     EQ3NR/wr3pkd.f.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
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
      integer ietmax,jetmax,kmax,nbtmax,netmax,nopgmx,ntitmx,nxmdmx
c
      integer newin
c
      integer iopg(nopgmx),jgext(netmax),jflgi(nbtmax),kxmod(nxmdmx),
     $ ngexrt(jetmax,netmax)
c
      integer kbt,kct,kdim,kmt,kxt,nbti,net,nobswt,nsbswt,ntitl2,nxmod
c
      logical qgexsh
c
      character(len=80) utitl2(ntitmx)
      character(len=56) ugexr(ietmax,jetmax,netmax)
      character(len=48) ubmtbi(nbtmax),uobsw(2,nbtmax),usbsw(2,nbtmax),
     $ uxmod(nxmdmx),uzveci(kmax)
      character(len=24) ugexmo(netmax),ugexp(netmax)
      character(len=8) ugexj(jetmax,netmax),
     $ uhfgex(ietmax,jetmax,netmax),uvfgex(ietmax,jetmax,netmax),
     $ uxkgex(ietmax,jetmax,netmax)
c
      real(8) cgexj(jetmax,netmax),mtbaqi(nbtmax),mtbi(nbtmax),
     $ mwtges(netmax),tgexp(netmax),xhfgex(ietmax,jetmax,netmax),
     $ xlkgex(ietmax,jetmax,netmax),xlkmod(nxmdmx),
     $ xvfgex(ietmax,jetmax,netmax),zgexj(jetmax,netmax),zvclgi(kmax)
c
      real(8) electr,pressi,tempci
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,je,j2,j3,j4,kcol,kprs,krow,n,nbi,ne
c
      integer ilnobl
c
      character(len=8) uendit
c
c-----------------------------------------------------------------------
c
      data uendit /'endit.  '/
c
c-----------------------------------------------------------------------
c
c     The items in this block correspond to things defined in the
c     first part of XCON6/wr6w8.f or EQ6/wr6pkw.f.
c
      j3 = ilnobl(uendit)
c
c     Comment header for iopt, iopr, iodb, and iopg option switches.
c
 1200 format('*',15x,'1    2    3    4    5    6    7    8    9   10')
c
 1020 format(a)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Old title.
c
      write (newin,2590)
 2590 format('* Title of previous run')
c
      do n = 1,ntitl2
        j2 = ilnobl(utitl2(n))
        write (newin,1020) utitl2(n)(1:j2)
      enddo
      if (ntitl2 .lt. ntitmx) write (newin,1020) uendit(1:j3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Special basis switches.
c
      write (newin,2610)
 2610 format('* Special basis switches')
c
      write (newin,2630) nsbswt
 2630 format(4x,'nsbswt= ',i3)
c
      do n = 1,nsbswt
        j2 = ilnobl(usbsw(1,n))
        write (newin,2650) usbsw(1,n)(1:j2)
 2650   format('species= ',a)
        j2 = ilnobl(usbsw(2,n))
        write (newin,2670) usbsw(2,n)(1:j2)
 2670   format(2x,'switch with= ',a)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original temperature.
c
      write (newin,2680)
 2680 format('* Original temperature and pressure')
      write (newin,2690) tempci
 2690 format(4x,'tempci= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original pressure.
c
      write (newin,2710) pressi
 2710 format(4x,'pressi= ',1pe12.5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ion exchanger creation.
c
      write (newin,2790)
 2790 format('* Ion exchangers')
c
      write (newin,2810) net
 2810 format(7x,'net= ',i3)
c
      do ne = 1,net
        j2 = ilnobl(ugexp(ne))
        write (newin,2840) ugexp(ne)(1:j2)
 2840   format(5x,'ugexp= ',a)
c
        write (newin,2860) mwtges(ne)
 2860   format(4x,'mwtges= ',1pe12.5)
c
        j3 = ilnobl(ugexmo(ne))
        write (newin,2870) ugexmo(ne)(1:j3)
 2870   format(4x,'ugexmo= ',a)
c
        write (newin,2880) tgexp(ne)
 2880   format(5x,'tgexp= ',1pe12.5)
c
        write (newin,2900) jgext(ne)
 2900   format(5x,'jgext= ',i3)
c
        do je = 1,jgext(ne)
          j3 = ilnobl(ugexj(je,ne))
          write (newin,2940) ugexj(je,ne)(1:j3)
 2940     format(5x,'ugexj= ',a)
c
          write (newin,2960) cgexj(je,ne),zgexj(je,ne)
 2960     format(5x,'cgexj= ',1pe12.5,5x,'zgexj= ',e12.5)
c
          write (newin,3010) ngexrt(je,ne)
 3010     format(4x,'ngexrt= ',i3)
c
          do n = 1,ngexrt(je,ne)
            j2 = ilnobl(ugexr(n,je,ne))
            write (newin,3040) ugexr(n,je,ne)(1:j2)
 3040       format(5x,'ugexr= ',a)
            j4 = ilnobl(uxkgex(n,je,ne))
            write (newin,3060) xlkgex(n,je,ne),uxkgex(n,je,ne)(1:j4)
 3060       format(4x,'xlkgex= ',1pe12.5,5x,'units= ',a)
            j4 = ilnobl(uhfgex(n,je,ne))
            write (newin,3070) xhfgex(n,je,ne),uhfgex(n,je,ne)(1:j4)
 3070       format(4x,'xhfgex= ',1pe12.5,5x,'units= ',a)
            j4 = ilnobl(uvfgex(n,je,ne))
            write (newin,3080) xvfgex(n,je,ne),uvfgex(n,je,ne)(1:j4)
 3080       format(4x,'xvfgex= ',1pe12.5,5x,'units= ',a)
          enddo
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Number of nxmod options.
c
      write (newin,3110)
 3110 format('* Alter/suppress options')
c
      write (newin,3130) nxmod
 3130 format(5x,'nxmod= ',i2)
c
c     Nxmod options.
c
      do n = 1,nxmod
        j2 = ilnobl(uxmod(n))
        write (newin,3160) uxmod(n)(1:j2),kxmod(n),xlkmod(n)
 3160   format(3x,'species= ',a,
     $  /4x,'option= ',i2,14x,'xlkmod= ',1pe12.5)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Iopg options.
c     Note: iopg(1) = iopg1, etc.
c
      write (newin,3172)
 3172 format('* Iopg options')
      write (newin,1200)
c
      write (newin,3180) (iopg(i), i = 1,20)
 3180 format(2x,'iopg1-10= ',10i5,/1x,'iopg11-20= ',10i5)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Index limits.
c
      write (newin,3210)
 3210 format('* Matrix index limits')
c
      kprs = 0
      write (newin,3230) kct,kbt,kmt,kxt,kdim,kprs
 3230 format(7x,'kct= ',i2,17x,'kbt= ',i2,17x,'kmt= ',i2,/
     $ 7x,'kxt= ',i2,16x,'kdim= ',i2,16x,'kprs= ',i2)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Species for which mass balances are defined.
c
c       ubmtbi = names of the data file basis species for which mass
c                  balances are defined
c       jflgi  = jflag input for basis species
c          0 = Retain as an active basis species
c         30 = Convert to a dependent species; fold the mass balance
c                total for this species into the mass balance totals
c                of basis species which remain active
c
      write (newin,3330)
 3330 format('* Data file basis species and jflgi values')
c
      do nbi = 1,nbti
        write (newin,3350) ubmtbi(nbi),jflgi(nbi)
 3350   format(3x,a48,3x,i2)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Mass balance totals.
c
      write (newin,3400)
 3400 format('* Mass balance totals (ES and aqueous solution)')
c
      do nbi = 1,nbti
        write (newin,3420) mtbi(nbi),mtbaqi(nbi)
 3420   format(3x,1pe25.15,3x,e25.15)
      enddo
c
      write (newin,3440) electr
 3440 format(9x,'Electrical imbalance= ',1pe25.15)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ordinary basis switches.
c
      write (newin,3450)
 3450 format('* Ordinary basis switches')
c
      write (newin,3460) nobswt
 3460 format(4x,'nobswt= ',i3)
c
      do n = 1,nobswt
        j2 = ilnobl(uobsw(1,n))
        write (newin,2650) uobsw(1,n)(1:j2)
        j2 = ilnobl(uobsw(2,n))
        write (newin,2670) uobsw(2,n)(1:j2)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Matrix column variables.
c
      write (newin,3470)
 3470 format('* Matrix species or entities')
c
      do krow = 1,kdim
        j2 = ilnobl(uzveci(krow))
        write (newin,3490) uzveci(krow)(1:j2)
 3490   format(3x,a)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Matrix variable values.
c
      write (newin,3510)
 3510 format('* Values of matrix variables')
c
      do kcol = 1,kdim
        write (newin,3530) zvclgi(kcol)
 3530   format(3x,1pe22.15)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The is no PRS in EQ3NR, so the block of data for phases and
c     species in the PRS is not written here.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
