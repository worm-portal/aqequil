      subroutine satchk(acflg,act,actlg,afcnst,affp,affs,apx,bpx,
     $ cdrs,eps100,iindx1,iodb,iopt,iapxmx,ibpxmx,iktmax,ixrn1,jflag,
     $ jpflag,jsflag,jsol,kmax,km1,kpsat,kpsst,kxt,nbasp,nbt,nbtmax,
     $ ncmpr,ndrs,ndrsmx,ndrsr,nodbmx,noptmx,noutpt,npchk,npt,nptmax,
     $ nstmax,nttyo,nxrn1,nxrn2,nxtmax,qxknph,sidrph,sidrsp,tolsat,
     $ uphase,uspec,wfac,xbar,xbarlg,xlks)
c
c     This subroutine checks for newly saturated phases and
c     supersaturated phases.
c
c     This subroutine is called by:
c
c       EQ6/eqphas.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c       affp   = array of affinities associated with phases
c       affs   = array of affinities associated with species
c       jpflag = status flag array for phases
c       kpsat  = number of newly saturated phases
c       kpsst  = number of supersaturated phases
c       sidrph = array of saturation indices (log Q/K) associated
c                  with phases
c       sidrsp = array of saturation indices (log Q/K) associated
c                  with species
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer iapxmx,ibpxmx,iktmax,kmax,nbtmax,ndrsmx,nodbmx,noptmx,
     $ nptmax,nstmax,nxtmax
c
      integer noutpt,nttyo
c
      integer iindx1(kmax),iodb(nodbmx),iopt(noptmx),jflag(nstmax),
     $ jpflag(nptmax),jsflag(nstmax),jsol(nxtmax),nbasp(nbtmax),
     $ ncmpr(2,nptmax),ndrs(ndrsmx),ndrsr(2,nstmax),npchk(nptmax)
c
      integer ixrn1,km1,kpsat,kpsst,kxt,nbt,npt,nxrn1,nxrn2
c
      logical qxknph(nptmax)
c
      character*48 uspec(nstmax)
      character*24 uphase(nptmax)
c
      real*8 acflg(nstmax),act(nstmax),actlg(nstmax),affp(nptmax),
     $ affs(nstmax),apx(iapxmx,nxtmax),bpx(ibpxmx,nxtmax),cdrs(ndrsmx),
     $ sidrph(nptmax),sidrsp(nstmax),wfac(iktmax,nxtmax),xbar(nstmax),
     $ xbarlg(nstmax),xlks(nstmax)
c
      real*8 afcnst,eps100,tolsat
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ier,j2,kcol,nb,np,nr1,nr2,ns,nt
c
      integer ilnobl
c
      logical qtest
c
      character*24 uaqsln,ugas
c
      real*8 af,si
c
c-----------------------------------------------------------------------
c
      data uaqsln /'Aqueous solution        '/
     $     ugas   /'Gas                     '/
c
c-----------------------------------------------------------------------
c
cXX   This double-checking of the jsflag and jpflag arrays shouldn't be
cXX   necessary.
c     Double-check flags on status of non-aqueous phases which are
c     in the matrix.
c
      do np = 1,npt
        if (uphase(np)(1:24) .ne. uaqsln(1:24)) then
          if (jpflag(np) .le. 0) then
            jpflag(np) = 0
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)
            do ns = nr1,nr2
              if (jsflag(ns) .le. 0) then
                jsflag(ns) = 0
                do kcol = km1,kxt
                  if (ns .eq. iindx1(kcol)) then
                    jsflag(ns) = -1
                    jpflag(np) = -1
                    go to 30
                  endif
                enddo
              endif
   30       continue
            enddo
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      kpsat = 0
      kpsst = 0
c
c     Compute affinities.
c
      do np = 1,npt
        nr1 = ncmpr(1,np)
        nr2 = ncmpr(2,np)
        nt = nr2 - nr1 + 1
c
        if (uphase(np)(1:24) .eq. ugas(1:24)) then
c
c         Case of a gas phase.
c
          affp(np) = -9999999.
          sidrph(np) = -9999999.
          do ns = nr1,nr2
            affs(ns) = -9999999.
            sidrsp(ns) = -9999999.
          enddo
c
        elseif (nt .eq. 1) then
c
c         Case of a pure phase.
c
          ns = nr1
          if (jpflag(np) .le. 1) then
            call afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,ndrsmx,
     $      ndrsr,ns,nstmax,si,xlks)
            affs(ns) = af
            affp(np) = af
          else
            af = -9999999.
            si = -9999999.
          endif
          affs(ns) = af
          affp(np) = af
          sidrsp(ns) = si
          sidrph(np) = si
c
        elseif (uphase(np)(1:24) .eq. uaqsln(1:24)) then
c
c         Case of an aqueous solution.
c
          if (jpflag(np) .le. 1) then
            affp(np) = 0.
            sidrph(np) = 0.
            do nb = 1,nbt
              ns = nbasp(nb)
              if (jsflag(ns) .le. 1) then
                call afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,
     $          ndrsmx,ndrsr,ns,nstmax,si,xlks)
                affs(ns) = af
                sidrsp(ns) = si
                affp(np) = affp(np) + xbar(ns)*affs(ns)
                sidrph(np) = sidrph(np) + xbar(ns)*sidrsp(ns)
              else
                affs(ns) = -9999999.
                sidrsp(ns) = -9999999.
              endif
            enddo
          else
            affp(np) = -9999999.
            sidrph(np) = -9999999.
          endif
c
        elseif (iopt(4) .ge. 1) then
c
c         If not doing non-aqueous solutions, skip.
c
          if (jpflag(np) .eq. -1) then
c
c           Case of non-aqueous solutions in the ES.
c
            qxknph(np) = .true.
            affp(np) = 0.
            sidrph(np) = 0.
            qtest = .false.
            do ns = nr1,nr2
              if (jsflag(ns) .le. 1) then
                call afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,
     $          ndrsmx,ndrsr,ns,nstmax,si,xlks)
                affs(ns) = af
                sidrsp(ns) = si
                affp(np) = affp(np) + xbar(ns)*affs(ns)
                sidrph(np) = sidrph(np) + xbar(ns)*sidrsp(ns)
                qtest = .true.
              else
                affs(ns) = -9999999.
                sidrsp(ns) = -9999999.
              endif
            enddo
            if (.not.qtest) then
              affp(np) = -9999999.
              sidrph(np) = -9999999.
            endif
          elseif (jpflag(np) .le. 1) then
c
c           Case of non-aqueous solutions not in the ES.
c
            call hpsat(acflg,act,actlg,afcnst,affp,affs,apx,bpx,
     $      cdrs,eps100,iapxmx,ibpxmx,ier,iktmax,ixrn1,jflag,jpflag,
     $      jsflag,jsol,ncmpr,ndrs,ndrsmx,ndrsr,noutpt,np,nptmax,
     $      nstmax,nttyo,nxrn1,nxrn2,nxtmax,sidrsp,sidrph,uphase,
     $      uspec,wfac,xbar,xbarlg,xlks)
c
            if (ier .le. 0) then
              qxknph(np) = .true.
            else
              qxknph(np) = .false.
              j2 = ilnobl(uphase(np))
              write (noutpt,2005) uphase(np)(1:j2)
              write (nttyo,2005) uphase(np)(1:j2)
 2005         format(/' * Warning - (EQ6/satchk) A hypothetical',
     $        ' affinity calculation',/7x,'failed for solid',
     $        ' solution ',a,'.')
              affp(np) = -9999999.
              sidrph(np) = -9999999.
              do ns = nr1,nr2
                affs(ns) = -9999999.
                sidrsp(ns) = -9999999.
              enddo
            endif
          endif
        endif
c
c       Set the jpflag array to mark supersaturated and newly saturated
c       phases.
c
        if (jpflag(np).le.0 .and. jpflag(np).ne.-1 .and.
     $    npchk(np) .ne. 1) then
          if (abs(affp(np)) .le. tolsat) then
            jpflag(np) = -10
            kpsat = kpsat + 1
          elseif (affp(np) .gt. tolsat) then
            jpflag(np) = -2
            kpsst = kpsst + 1
          endif
        endif
c
c       Set the jsflag array to mark the corresponding species.
c
        if (uphase(np)(1:24) .eq. uaqsln(1:24)) then
          do nb = 1,nbt
            ns = nbasp(nb)
            if (jsflag(ns).le.0 .and. jsflag(ns).ne.-1) then
              if (abs(affs(ns)) .le. tolsat) then
                jsflag(ns) = -10
              elseif (affs(ns) .gt. tolsat) then
                jsflag(ns) = -2
              endif
            endif
          enddo
        else
          do ns = nr1,nr2
            if (jsflag(ns).le.0 .and. jsflag(ns).ne.-1) then
              if (abs(affs(ns)) .le. tolsat) then
                jsflag(ns) = -10
              elseif (affs(ns) .gt. tolsat) then
                jsflag(ns) = -2
              endif
            endif
          enddo
        endif
c
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
