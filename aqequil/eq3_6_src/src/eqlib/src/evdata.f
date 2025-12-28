      subroutine evdata(aadh,aadhh,aadhv,aaphi,abdh,abdhh,abdhv,
     $ abdot,abdoth,abdotv,adadhh,adadhv,adbdhh,adbdhv,adbdth,adbdtv,
     $ adh,adhfe,adhh,adhv,adhfs,adhfsd,advfe,advfs,advfsd,afcnst,
     $ al10,amu,aslm,aphi,aprehw,apresg,apresh,apx,avcnst,axhfe,
     $ axhfs,axhfsd,axlke,axlks,axlksd,axvfe,axvfs,axvfsd,bdh,bdhh,
     $ bdhv,bdot,bdoth,bdotv,dadhh,dadhv,dbdhh,dbdhv,dbdth,dbdtv,
     $ dhfe,dhfs,dhfsd,dvfe,dvfs,dvfsd,ehfac,farad,iapxmx,iktmax,
     $ iopg,iopt,ipbtmx,ipch,ipchmx,ipcv,ipcvmx,ixrn1,ixrn2,jpfcmx,
     $ jptffl,jsol,narxmx,narxt,narxth,ncmpr,nmut,nmutmx,nopgmx,
     $ noptmx,noutpt,nptmax,nslt,nsltmx,nst,nstmax,ntpr,ntprmx,nttyo,
     $ nxt,nxtmax,pmu,prehw,presg,presh,press,pslamn,rconst,rcnstv,
     $ rtcnst,tempc,tempk,uphase,uspec,wfac,xhfe,xhfs,xhfsd,xlke,xlks,
     $ xlksd,xvfe,xvfs,xvfsd)
c
c     This subroutine evaluates thermodynamic properties as a function
c     of temperature and the standard grid pressure. Pressure
c     corrections are made by EQLIB/pcorrx.f.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c       EQ6/tpadv.f
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
      integer iapxmx,iktmax,ipbtmx,ipchmx,ipcvmx,jpfcmx,narxmx,nmutmx,
     $ nopgmx,noptmx,nptmax,nsltmx,nstmax,ntprmx,nxtmax
c
      integer noutpt,nttyo
c
      integer iopg(nopgmx),iopt(noptmx),jsol(nxtmax),narxt(ntprmx),
     $ narxth(2),ncmpr(2,nptmax)
c
      integer ipch,ipcv,ixrn1,ixrn2,jptffl,nmut,nslt,nst,ntpr,nxt
c
      character*48 uspec(nstmax)
      character*24 uphase(nptmax)
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
      real*8 afcnst,al10,avcnst,ehfac,farad,rconst,rcnstv,rtcnst,
     $ tempc,tempk
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ipc,ntprh
c
      real*8 prop
c
c-----------------------------------------------------------------------
c
c     Compute temperature dependent constants.
c
      rtcnst = 0.001*rconst*tempk
      afcnst = al10*rtcnst
      avcnst = al10*rcnstv*tempk
      ehfac = (al10*rconst*tempk)/farad
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the standard grid pressure.
c
c     Calling sequence substitutions:
c       apresg for arr
c       presg for prop
c
      call evdat2(apresg,narxmx,narxt,ntpr,ntprmx,presg,tempc)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the 1.013-bar/steam-saturation curve pressure at the
c     initial temperature.
c
      if (tempc .le. 100.) then
        ntprh = 1
      else
        ntprh = 2
      endif
c
c     Calling sequence substitutions:
c       apresh for arr
c       5 for narxmx
c       narxth for narxt
c       ntprh for ntpr
c       2 for ntprmx
c       presh for prop
c
      call evdat2(apresh,5,narxth,ntprh,2,presh,tempc)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ipcv .ge. 0) then
c
c       Compute the half-width of the recommended pressure envelope.
c
c       Calling sequence substitutions:
c         aprehw for arr
c         prehw for prop
c
        call evdat2(aprehw,narxmx,narxt,ntpr,ntprmx,prehw,tempc)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopg(1) .le. 0) then
c
c       Davies equation or B-dot equation.
c
c       Compute the Debye-Huckel A(gamma,10) parameter.
c
c       Calling sequence substitutions:
c         aadh for arr
c         adh for prop
c
        call evdat2(aadh,narxmx,narxt,ntpr,ntprmx,adh,tempc)
c
        aphi = adh*al10/3.
      endif
c
      if (iopg(1) .eq. 1) then
c
c       Pitzer's equations.
c
c       Compute the Debye-Huckel A(phi) parameter.
c
c       Calling sequence substitutions:
c         aaphi for arr
c         aphi for prop
c
        call evdat2(aaphi,narxmx,narxt,ntpr,ntprmx,aphi,tempc)
        adh= 3.*aphi/al10
      endif
c
      if (iopg(1) .eq. 2) then
c
c       HC + DH equations.
c
c       Compute the Debye-Huckel A(gamma,10) parameter.
c
c       Calling sequence substitutions:
c         aadh for arr
c         adh for prop
c
        call evdat2(aadh,narxmx,narxt,ntpr,ntprmx,adh,tempc)
c
        aphi = adh*al10/3.
      endif
c
      if (ipch .ge. 0) then
c
c       Compute the Debye-Huckel A(H) parameter.
c
c       Calling sequence substitutions:
c         aadhh for arr
c         adhh for prop
c
        call evdat2(aadhh,narxmx,narxt,ntpr,ntprmx,adhh,tempc)
c
        do ipc = 1,ipch
c
c         Compute the pressure derivatives of the Debye-Huckel A(H)
c         parameter.
c
c         Calling sequence substitutions:
c           adadhh for arr
c           ipc for k
c           ipchmx for nmax
c
          call evdat3(adadhh,ipc,ipchmx,narxmx,narxt,ntpr,ntprmx,
     $    prop,tempc)
          dadhh(ipc) = prop
        enddo
      endif
c
      if (ipcv .ge. 0) then
c
c       Compute the Debye-Huckel A(V) parameter.
c
c       Calling sequence substitutions:
c         aadhv for arr
c         adhv for prop
c
        call evdat2(aadhv,narxmx,narxt,ntpr,ntprmx,adhv,tempc)
c
        do ipc = 1,ipcv
c
c         Compute the pressure derivatives of the Debye-Huckel A(V)
c         parameter.
c
c         Calling sequence substitutions:
c           adadhv for arr
c           ipc for k
c           ipcvmx for nmax
c
          call evdat3(adadhv,ipc,ipcvmx,narxmx,narxt,ntpr,ntprmx,
     $    prop,tempc)
          dadhv(ipc) = prop
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopg(1).eq.0 .or. iopg(1).eq.2) then
c
c       B-dot equation or HC + DH equations.
c
c       Compute the Debye-Huckel B(gamma) parameter.
c
c       Calling sequence substitutions:
c         abdh for arr
c         bdh for prop
c
        call evdat2(abdh,narxmx,narxt,ntpr,ntprmx,bdh,tempc)
c
        if (ipch .ge. 0) then
c
c         Compute the Debye-Huckel B(H) parameter.
c
c         Calling sequence substitutions:
c           abdhh for arr
c           bdhh for prop
c
          call evdat2(abdhh,narxmx,narxt,ntpr,ntprmx,bdhh,tempc)
c
          do ipc = 1,ipch
c
c           Compute the pressure derivatives of the Debye-Huckel B(H)
c           parameter.
c
c           Calling sequence substitutions:
c             adbdhh for arr
c             ipc for k
c             ipchmx for nmax
c
            call evdat3(adbdhh,ipc,ipchmx,narxmx,narxt,ntpr,ntprmx,
     $      prop,tempc)
            dbdhh(ipc) = prop
          enddo
        endif
c
        if (ipcv .ge. 0) then
c
c         Compute the Debye-Huckel B(V) parameter.
c
c         Calling sequence substitutions:
c           abdhv for arr
c           bdhv for prop
c
          call evdat2(abdhv,narxmx,narxt,ntpr,ntprmx,bdhv,tempc)
c
          do ipc = 1,ipcv
c
c           Compute the pressure derivatives of the Debye-Huckel B(V)
c           parameter.
c
c           Calling sequence substitutions:
c             adbdhv for arr
c             ipc for k
c             ipcvmx for nmax
c
            call evdat3(adbdhv,ipc,ipcvmx,narxmx,narxt,ntpr,ntprmx,
     $      prop,tempc)
            dbdhv(ipc) = prop
          enddo
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopg(1) .eq. 0) then
c
c       B-dot equation.
c
c       Compute the Helgeson (1969) B-dot parameter.
c
c       Calling sequence substitutions:
c         abdot for arr
c         bdot for prop
c
        call evdat2(abdot,narxmx,narxt,ntpr,ntprmx,bdot,tempc)
c
        if (ipch .ge. 0) then
c
c         Compute the Debye-Huckel B-dot(H) parameter.
c
c         Calling sequence substitutions:
c           abdoth for arr
c           bdoth for prop
c
          call evdat2(abdoth,narxmx,narxt,ntpr,ntprmx,bdoth,tempc)
c
          do ipc = 1,ipch
c
c           Compute the pressure derivatives of the Debye-Huckel
c           B-dot(H) parameter.
c
c           Calling sequence substitutions:
c             adbdth for arr
c             ipc for k
c             ipchmx for nmax
c
            call evdat3(adbdth,ipc,ipchmx,narxmx,narxt,ntpr,ntprmx,
     $      prop,tempc)
            dbdth(ipc) = prop
          enddo
        endif
c
        if (ipcv .ge. 0) then
c
c         Compute the Debye-Huckel B-dot(V) parameter.
c
c         Calling sequence substitutions:
c           abdotv for arr
c           bdotv for prop
c
          call evdat2(abdotv,narxmx,narxt,ntpr,ntprmx,bdotv,tempc)
c
          do ipc = 1,ipcv
c
c           Compute the pressure derivatives of the Debye-Huckel
c           B-dot(V) parameter.
c
c           Calling sequence substitutions:
c             adbdtv for arr
c             ipc for k
c             ipcvmx for nmax
c
            call evdat3(adbdtv,ipc,ipcvmx,narxmx,narxt,ntpr,ntprmx,
     $      prop,tempc)
            dbdtv(ipc) = prop
          enddo
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the log K for the 'Eh' reaction.
c
      if (axlke(1,ntpr) .lt. 9999999.) then
c
c       Calling sequence substitutions:
c         axlke for arr
c         xlke for prop
c
        call evdat2(axlke,narxmx,narxt,ntpr,ntprmx,xlke,tempc)
      else
        xlke = 9999999.
      endif
c
      if (ipch .ge. 0) then
c
c       Compute the enthalpy of reaction for the 'Eh' reaction.
c
        if (axhfe(1,ntpr) .lt. 9999999.) then
c
c         Calling sequence substitutions:
c           axhfe for arr
c           xhfe for prop
c
          call evdat2(axhfe,narxmx,narxt,ntpr,ntprmx,xhfe,tempc)
        else
          xhfe = 9999999.
        endif
c
        do ipc = 1,ipch
c
c         Compute the pressure derivatives of the enthalpy of
c         reaction for the 'Eh' reaction.
c
c         Calling sequence substitutions:
c           adhfe for arr
c           ipc for k
c           ipchmx for nmax
c
          call evdat3(adhfe,ipc,ipchmx,narxmx,narxt,ntpr,ntprmx,
     $    prop,tempc)
          dhfe(ipc) = prop
        enddo
      endif
c
      if (ipcv .ge. 0) then
c
c       Compute the volume of reaction for the 'Eh' reaction.
c
        if (axvfe(1,ntpr) .lt. 9999999.) then
c
c         Calling sequence substitutions:
c           axvfe for arr
c           xvfe for prop
c
          call evdat2(axvfe,narxmx,narxt,ntpr,ntprmx,xvfe,tempc)
        else
          xvfe = 9999999.
        endif
c
        do ipc = 1,ipcv
c
c         Compute the pressure derivatives of the volume of
c         reaction for the 'Eh' reaction.
c
c         Calling sequence substitutions:
c           advfe for arr
c           ipc for k
c           ipcvmx for nmax
c
          call evdat3(advfe,ipc,ipcvmx,narxmx,narxt,ntpr,ntprmx,
     $    prop,tempc)
          dvfe(ipc) = prop
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the log K values for all reactions as they are
c     currently written.
c
      call evdatr(adhfs,advfs,axhfs,axlks,axvfs,dhfs,dvfs,
     $ ipch,ipchmx,ipcv,ipcvmx,narxmx,narxt,nst,nstmax,ntpr,ntprmx,
     $ tempc,xhfs,xlks,xvfs)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the log K values for all reactions as they were
c     written on the data file.
c
c     Calling sequence substitutions:
c       adhfsd for adhfs
c       advfsd for advfs
c       axhfsd for axhfs
c       axlksd for axlks
c       axvfsd for avhfs
c       dhfsd for dhfs
c       dvfsd for dvfs
c       xhfsd for xhfs
c       xlksd for xlks
c       xvfsd for xvfs
c
      call evdatr(adhfsd,advfsd,axhfsd,axlksd,axvfsd,dhfsd,dvfsd,
     $ ipch,ipchmx,ipcv,ipcvmx,narxmx,narxt,nst,nstmax,ntpr,ntprmx,
     $ tempc,xhfsd,xlksd,xvfsd)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute remaining temperature-dependent aqueous species activity
c     coefficient parameters.
c
      if (iopg(1) .eq. 1) then
c
c       Computer the S-lambda(n) and mu coefficients for Pitzer's
c       equations.
c
        call evptzc(amu,aslm,ipbtmx,jpfcmx,jptffl,nmut,nmutmx,
     $  noutpt,nttyo,nslt,nsltmx,pmu,pslamn,tempc)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute temperature dependent solid solution activity
c     coefficient parameters.
c
      if (iopt(4) .ge. 1) then
        call wterm(apx,iapxmx,iktmax,ixrn1,ixrn2,jsol,ncmpr,
     $  noutpt,nptmax,nstmax,nttyo,nxt,nxtmax,press,tempk,uphase,
     $  uspec,wfac)
      endif
c
      end
