      subroutine tpadv(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,
     $ abdoth,abdot,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,
     $ adh,adhfe,adhh,adhv,adhfs,adhfsd,advfe,advfs,advfsd,afcnst,
     $ al10,amu,aslm,aphi,aprehw,apresg,apresh,apx,avcnst,axhfe,axhfs,
     $ axhfsd,axlke,axlks,axlksd,axvfe,axvfs,axvfsd,bdh,bdhh,bdhv,
     $ bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,dhfe,dhfs,
     $ dhfsd,dvfe,dvfs,dvfsd,eact,ehfac,farad,hact,iact,iapxmx,iktmax,
     $ imchmx,imech,iopg,iopt,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,ixrn1,
     $ ixrn2,jpfcmx,jpress,jptffl,jsol,jtemp,narxmx,narxt,narxth,nbasp,
     $ nbaspd,nbt,nbtd,nbtmax,ncmpr,ndrsr,ndrsrd,nmut,nmutmx,nopgmx,
     $ noptmx,noutpt,nptkmx,nptmax,nrct,nrctmx,nrk,nslt,nsltmx,nst,
     $ nstmax,ntpr,ntprmx,ntprt,nttkmx,nttyo,nweope,nwndpc,nxt,nxtmax,
     $ pmu,presg,presh,press,pressb,pressd,pslamn,ptk,rcnstv,rconst,rk,
     $ rkb,rtcnst,tempc,tempcb,tempcd,tempcu,tempk,time1,trkb,ttk,
     $ uphase,uspec,wfac,xhfe,xhfs,xhfsd,xi1,xlke,xlks,xlksd,xvfe,
     $ xvfs,xvfsd)
c
c     This subroutine changes the temperature and pressure and
c     recomputes all temperature-pressure dependent data. Here tempcd
c     and pressd are the last temperature and pressure, respectively,
c     at which the data were evaluated.
c
c     This subroutine is called by:
c
c       EQ6/path.f
c       EQ6/eqshel.f
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
      integer iapxmx,iktmax,imchmx,ipbtmx,ipchmx,ipcvmx,jpfcmx,narxmx,
     $ nbtmax,nmutmx,nopgmx,noptmx,nptkmx,nptmax,nrctmx,nsltmx,nstmax,
     $ ntprmx,nttkmx,nxtmax
c
      integer noutpt,nttyo
c
      integer iact(imchmx,2,nrctmx),imech(2,nrctmx),iopg(nopgmx),
     $ iopt(noptmx),jsol(nxtmax),narxt(ntprmx),narxth(2),nbasp(nbtmax),
     $ nbaspd(nbtmax),ncmpr(2,nptmax),ndrsr(2,nstmax),ndrsrd(2,nstmax),
     $ nrk(2,nrctmx)
c
      integer ipch,ipcv,ixrn1,ixrn2,jpress,jptffl,jtemp,nbt,nbtd,nmut,
     $ nrct,nslt,nst,ntpr,ntprt,nweope,nwndpc,nxt
c
      character*48 uspec(nstmax)
      character*24 uphase(nptmax)
c
      real*8 tempcu(ntprmx)
c
      real*8 aadh(narxmx,ntprmx),aadhh(narxmx,ntprmx),
     $ aadhv(narxmx,ntprmx),aaphi(narxmx,ntprmx),
     $ abdh(narxmx,ntprmx),abdhh(narxmx,ntprmx),abdhv(narxmx,ntprmx),
     $ abdot(narxmx,ntprmx),abdoth(narxmx,ntprmx),abdotv(narxmx,ntprmx),
     $ adadhh(narxmx,ntprmx,ipchmx),adadhv(narxmx,ntprmx,ipcvmx),
     $ adbdhh(narxmx,ntprmx,ipchmx),adbdhv(narxmx,ntprmx,ipcvmx),
     $ adbdth(narxmx,ntprmx,ipchmx),adbdtv(narxmx,ntprmx,ipcvmx),
     $ adhfe(narxmx,ntprmx,ipchmx),adhfs(narxmx,ntprmx,ipchmx,nstmax),
     $ adhfsd(narxmx,ntprmx,ipchmx,nstmax),advfe(narxmx,ntprmx,ipcvmx),
     $ advfs(narxmx,ntprmx,ipcvmx,nstmax),
     $ advfsd(narxmx,ntprmx,ipcvmx,nstmax),
     $ aprehw(narxmx,ntprmx),apresg(narxmx,ntprmx),apresh(5,2),
     $ axhfe(narxmx,ntprmx),axhfs(narxmx,ntprmx,nstmax),
     $ axhfsd(narxmx,ntprmx,nstmax),axlke(narxmx,ntprmx),
     $ axlks(narxmx,ntprmx,nstmax),axlksd(narxmx,ntprmx,nstmax),
     $ axvfe(narxmx,ntprmx),axvfs(narxmx,ntprmx,nstmax),
     $ axvfsd(narxmx,ntprmx,nstmax)
c
      real*8 dadhh(ipchmx),dadhv(ipcvmx),dbdhh(ipchmx),dbdhv(ipcvmx),
     $ dbdth(ipchmx),dbdtv(ipcvmx),dhfe(ipchmx),dhfs(ipchmx,nstmax),
     $ dhfsd(ipchmx,nstmax),dvfe(ipcvmx),dvfs(ipcvmx,nstmax),
     $ dvfsd(ipcvmx,nstmax),xhfs(nstmax),xhfsd(nstmax),xlks(nstmax),
     $ xlksd(nstmax),xvfs(nstmax),xvfsd(nstmax)
c
      real*8 adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,bdotv,
     $ prehw,presg,presh,press,xhfe,xlke,xvfe
c
      real*8 amu(jpfcmx,nmutmx),aslm(jpfcmx,0:ipbtmx,nsltmx),
     $ pmu(nmutmx),pslamn(0:ipbtmx,nsltmx)
c
      real*8 apx(iapxmx,nxtmax),wfac(iktmax,nxtmax)
c
      real*8 eact(imchmx,2,nrctmx),hact(imchmx,2,nrctmx),
     $ ptk(nptkmx),rk(imchmx,2,nrctmx),rkb(imchmx,2,nrctmx),
     $ trkb(imchmx,2,nrctmx),ttk(nttkmx)
c
      real*8 afcnst,al10,avcnst,ehfac,farad,pressb,pressd,rcnstv,
     $ rconst,rtcnst,tempc,tempcb,tempcd,tempk,time1,xi1
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nwpclm
c
      logical qnewp,qnewt
c
      real*8 dp,dt,pxl,pxu,toldt,toldp
c
c-----------------------------------------------------------------------
c
c     nwpclm  = maximum number of repeated warnings regarding problems
c                associated with pressure corrections
c     nwndpc = number of warnings of no data to support pressure
c                corrections
c     nweope = number of warnings of excursions outside the recommended
c                pressure envelope
c
      data nwpclm /5/
c
c     toldt = tolerance for recalculating temperature-dependent data (C)
c     toldp = tolerance for recalculating pressure-dependent data (bars)
c
      data toldt / 1.e-4 /, toldp /1.e-4 /
c
c-----------------------------------------------------------------------
c
c     Compute the new temperature.
c
      call gtemp(afcnst,al10,iopt,jtemp,noptmx,noutpt,nttkmx,
     $ nttyo,rconst,rtcnst,tempc,tempcb,tempk,time1,ttk,xi1)
c
      dt = tempc - tempcd
      qnewt = .false.
      if (abs(dt) .gt. toldt) then
c
c       Determine the corresponding temperature range flag.
c
        call gntpr(ntpr,ntprmx,ntprt,tempc,tempcu)
c
        qnewt = .true.
c
c       Compute thermodynamic data for the current temperature.
c
        call evdata(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,
     $  abdot,abdoth,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,
     $  adh,adhfe,adhh,adhv,adhfs,adhfsd,advfe,advfs,advfsd,afcnst,
     $  al10,amu,aslm,aphi,aprehw,apresg,apresh,apx,avcnst,axhfe,
     $  axhfs,axhfsd,axlke,axlks,axlksd,axvfe,axvfs,axvfsd,bdh,bdhh,
     $  bdhv,bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,
     $  dhfe,dhfs,dhfsd,dvfe,dvfs,dvfsd,ehfac,farad,iapxmx,iktmax,
     $  iopg,iopt,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,ixrn1,ixrn2,jpfcmx,
     $  jptffl,jsol,narxmx,narxt,narxth,ncmpr,nmut,nmutmx,nopgmx,
     $  noptmx,noutpt,nptmax,nslt,nsltmx,nst,nstmax,ntpr,ntprmx,nttyo,
     $  nxt,nxtmax,pmu,prehw,presg,presh,press,pslamn,rconst,rcnstv,
     $  rtcnst,tempc,tempk,uphase,uspec,wfac,xhfe,xhfs,xhfsd,xlke,xlks,
     $  xlksd,xvfe,xvfs,xvfsd)
c
c       Compute kinetic data for the current temperature.
c
        call evratc(eact,hact,iact,imchmx,imech,nrct,nrctmx,nrk,
     $  rk,rkb,rtcnst,tempk,trkb)
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the new pressure.
c
      call gpress(iopt,jpress,noptmx,noutpt,nptkmx,nttyo,presg,
     $ presh,press,pressb,time1,ptk,xi1)
c
c     Skip pressure corrections if the pressure is tracking on the
c     data file reference pressure curve.
c
      if (jpress .le. 0) go to 990
c
c     Skip pressure corrections if the temperature and pressure
c     haven't changed significantly.
c
      dp = press - pressd
      qnewp = abs(dp) .gt. toldp
      if (.not.qnewt .and. .not.qnewp) go to 990
c
c     Skip pressure corrections if the pressure happens to match the
c     data file reference pressure curve value.
c
      dp = press - presg
      if (abs(dp) .le. toldp) go to 990
c
c     Make pressure corrections, if the requisite data are available.
c
      if (ipcv .lt. 0) then
c
c       There are no data to support needed pressure corrections.
c
        if (nwndpc .le. nwpclm) then
          write (noutpt,1000) press,presg,dp
          write (nttyo,1000) press,presg,dp
 1000     format(/' * Warning - (EQ6/tpadv) The supporting data file',
     $    /7x,'contains no data to support making thermodynamic',
     $    /7x,'pressure corrections. No such corrections will be made.',
     $    /7x,'The current pressure is ',1pg12.5,' bars, the standard',
     $    /7x,'grid pressure is ',g12.5,' bars, and the pressure',
     $    /7x,'difference is ',g12.5,' bars.')
          if (nwndpc .eq. nwpclm) then
            write (noutpt,1010)
            write (nttyo,1010)
 1010       format(/' This warning will not be repeated during the',
     $      /7x,'rest of this run.')
          endif
          nwndpc = nwndpc + 1
        endif
      else
        if (abs(dp) .gt. prehw) then
c
c         The pressure is outside the recommended envelope.
c
          pxu = presg + prehw
          pxl = presg - prehw
          if (pxl .le. 0.) pxl = 0.
          write (noutpt,1020) press,pxl,pxu,tempc
          write (nttyo,1020) press,pxl,pxu,tempc
 1020     format(/' * Warning - (EQ6/tpadv) The current pressure',
     $    /7x,'of ',1pg12.5,' bars is outside the recommended',
     $    /7x,'pressure envelope of ',g12.5,' to ',g12.5,' bars',
     $    /7x,'at ',0pf6.2,' C.')
          if (nweope .eq. nwpclm) then
            write (noutpt,1010)
            write (nttyo,1010)
          endif
          nweope = nweope + 1
        endif
c
c       Make pressure corrections to the thermodynamic data.
c
        call pcorrx(avcnst,dhfs,dvfs,ipch,ipchmx,ipcv,ipcvmx,
     $  nbasp,nbt,nbtmax,ndrsr,nst,nstmax,presg,press,xhfs,
     $  xlks,xvfs)
c
c       Calling sequence substitutions:
c         dhfsd for dhfs
c         dvfsd for dvfs
c         nbaspd for nbasp
c         nbtd for nbt
c         ndrsrd for ndrsr
c         xhfsd for xhfs
c         xlksd for xlks
c         xvfsd for xvfs
c
        call pcorrx(avcnst,dhfsd,dvfsd,ipch,ipchmx,ipcv,ipcvmx,
     $  nbaspd,nbtd,nbtmax,ndrsrd,nst,nstmax,presg,press,xhfsd,
     $  xlksd,xvfsd)
c
c       Make pressure corrections to the kinetic data.
c
c         Presently there are no such corrections.
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 tempcd = tempc
      pressd = press
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
c
      end
