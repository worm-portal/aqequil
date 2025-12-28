      subroutine raff(acflg,actlg,afcnst,affp,afrc1,bpx,cdrs,cgexj,
     $ ibpxmx,ibpxt,iern1,ietmax,iktmax,ixrn1,ixrn2,jcode,jern1,jern2,
     $ jetmax,jflag,jgext,jpflag,jsflag,jsol,ncmpr,ndrs,ndrsmx,ndrsr,
     $ nertmx,net,netmax,ngext,noutpt,nptmax,nrct,nrctmx,nrndex,
     $ nstmax,nttyo,nxridx,nxrtmx,nxtmax,rxbar,uphase,uspec,wfac,
     $ xbar,xbarlg,xgers,xlks)
c
c     This subroutine calculates the affinities of irreversible
c     reactions (afrc1).
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
      include 'eqlib/eqlpar.h'
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ibpxmx,ietmax,iktmax,jetmax,ndrsmx,nertmx,netmax,nptmax,
     $ nrctmx,nstmax,nxrtmx,nxtmax
c
      integer noutpt,nttyo
c
      integer ibpxt(nxtmax),jcode(nrctmx),jern1(jetmax,netmax),
     $ jern2(jetmax,netmax),jflag(nstmax),jgext(netmax),jpflag(nptmax),
     $ jsflag(nstmax),jsol(nxtmax),ncmpr(2,nptmax),ndrs(ndrsmx),
     $ ndrsr(2,nstmax),ngext(jetmax,netmax),nrndex(nrctmx),
     $ nxridx(nrctmx)
c
      integer iern1,ixrn1,ixrn2,net,nrct
c
      character*48 uspec(nstmax)
      character*24 uphase(nptmax)
c
      real*8 acflg(nstmax),actlg(nstmax),affp(nptmax),afrc1(nrctmx),
     $ bpx(ibpxmx,nxtmax),cdrs(ndrsmx),cgexj(jetmax,netmax),
     $ rxbar(iktmax,nxrtmx),wfac(iktmax,nxtmax),xbar(nstmax),
     $ xbarlg(nstmax),xgers(ietmax,jetmax,nertmx),xlks(nstmax)
c
      real*8 afcnst
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with fixed global dimensioning.
c
      integer jflgx(iet_par,net_par),jsflgx(iet_par,net_par)
c
      real(8) acflge(iet_par,jet_par),actlge(iet_par,jet_par),
     $ xbarle(iet_par,jet_par),xbare(iet_par,jet_par)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with variable global dimensioning.
c
      integer isv_iktmax
c
      SAVE isv_iktmax
c
      real(8), dimension(:), allocatable :: acflgs,actlgs,xbarls,xbars
c
      SAVE acflgs,actlgs,xbarls,xbars
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ie,ik,je,jpflgx,k,ne,ner,np,nrc,nr1,nr2,ns,nxr
c
      real*8 af,affpr,si,xx
c
      real*8 tlg
c
c-----------------------------------------------------------------------
c
c     Allocate or reallocate local work arrays as needed.
c
      if (.not.ALLOCATED(acflgs)) then
c
c       Local work arrays are not allocated. Zero the saved
c       array size variables. Note that only one array is tested
c       to see if it is allocated. It is assumed that all local
c       work arrays are either allocated or not.
c
        isv_iktmax = 0
      else
c
c       Local work arrays are allocated. Check to see if any of the
c       array size variables have changed. If so, deallocate
c       the corresponding local work arrays and zero the corresponding
c       saved size variables.
c
        if (iktmax .ne. isv_iktmax) then
          DEALLOCATE(acflgs,actlgs,xbarls,xbars)
          isv_iktmax = 0
        endif
      endif
c
c     At this point, the saved array size values are zero if the
c     corresponding arrays need to be allocated.
c
      if (isv_iktmax .eq. 0) then
        ALLOCATE(acflgs(iktmax),actlgs(iktmax))
        ALLOCATE(xbarls(iktmax),xbars(iktmax))
        isv_iktmax = iktmax
      endif
c
c     Zero the contents of the local work arrays.
c
      do k = 1,iktmax
        acflgs(k) = 0.
        actlgs(k) = -99999.
        xbars(k) = 0.
        xbarls(k) = -99999.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Recall the following jreac reactant status flag conventions:
c
c       jreac =  0: set to react
c       jreac =  1: exhausted
c       jreac = -1: saturated, but the remaining reactant mass
c                     continues to react irreversibly
c       jreac =  2: saturated; the status of any remaining reactant
c                     mass is changed to that of a product phase
c
c     When jreac = -1 or 2, the affinity is zero by definition.
c     However, the corresponding computed affinities returned from the
c     present subroutine may be non-zero, owing to finite convergence
c     tolerances. Avoid calculating finite differences from such
c     computed affinity values in EQ6/stepfd.f. However, do not set
c     the corresponding afrc1 values to zero in the present subroutine.
c
      do nrc = 1,nrct
        if (jcode(nrc) .eq. 0) then
c
c          Pure minerals.
c
           np = nrndex(nrc)
           afrc1(nrc) = -affp(np)
c
        elseif (jcode(nrc) .eq. 1) then
c
c         Solid solutions.
c
          np = nrndex(nrc)
          nxr = nxridx(nrc)
c
          nr1 = ncmpr(1,np)
          nr2 = ncmpr(2,np)
c
          ik = 0
          do ns = nr1,nr2
            ik = ik + 1
            acflgs(ik) = acflg(ns)
            actlgs(ik) = actlg(ns)
            xbarls(ik) = xbarlg(ns)
            xbars(ik) = xbar(ns)
            xx = rxbar(ik,nxr)
            xbar(ns) = xx
            xbarlg(ns) = tlg(xx)
          enddo
c
c         Calling sequence substitutions:
c           acflg for acflgc
c
          call lambda(acflg,afcnst,bpx,ibpxmx,ibpxt,iktmax,ixrn1,
     $    ixrn2,jsol,ncmpr,noutpt,np,nptmax,nstmax,nttyo,nxtmax,wfac,
     $    xbar,xbarlg,uphase,uspec)
c
          do ns = nr1,nr2
            actlg(ns) = xbarlg(ns) + acflg(ns)
          enddo
c
          affpr = 0.
          ik = 0
          do ns = nr1,nr2
            ik = ik + 1
            xx = xbar(ns)
            if (xx .gt. 0.) then
              call afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,ndrsmx,
     $        ndrsr,ns,nstmax,si,xlks)
              if (af .gt. -9999999.) then
                affpr = affpr + xx*af
              else
                affpr = -9999999.
                go to 100
              endif
            endif
          enddo
  100     continue
          afrc1(nrc) = -affpr
c
          ik = 0
          do ns = nr1,nr2
            ik = ik + 1
            acflg(ns) = acflgs(ik)
            actlg(ns) = actlgs(ik)
            xbarlg(ns) = xbarls(ik)
            xbar(ns) = xbars(ik)
          enddo
c
        elseif (jcode(nrc) .eq. 2) then
c
c         Special reactants.
c
          afrc1(nrc) = 9999999.
c
        elseif (jcode(nrc) .eq. 3) then
c
c         Aqueous species.
c
          afrc1(nrc) = 9999999.
c
        elseif (jcode(nrc) .eq. 4) then
c
c         Gases.
c
          afrc1(nrc) = 9999999.
c
        elseif (jcode(nrc) .eq. 5) then
c
c         Generic ion exchangers.
c
          np = nrndex(nrc)
          ner = nxridx(nrc)
          ne = np - iern1 + 1
c
          nr1 = ncmpr(1,np)
          nr2 = ncmpr(2,np)
c
          jpflgx = jpflag(np)
          jpflag(np) = 0
c
          do je = 1,jgext(ne)
            ns = jern1(je,ne) - 1
            do ie = 1,ngext(je,ne)
              ns = ns + 1
              jsflgx(ie,je) = jsflag(ns)
              jsflag(ns) = 0
              jflgx(ie,je) = jflag(ns)
              jflag(ns) = 0
              acflge(ie,je) = acflg(ns)
              actlge(ie,je) = actlg(ns)
              xbarle(ie,je) = xbarlg(ns)
              xbare(ie,je) = xbar(ns)
              xx = xgers(ie,je,ner)
              xbar(ns) = xx
              xbarlg(ns) = tlg(xx)
            enddo
          enddo
c
c         Calculate activity coefficients. Recall that these may
c         non-zero even in the case of ideal exchangers, due to
c         site-mixing.
c
c         Calling sequence substitutions:
c           acflg for acflgc
c
          call lamgex(acflg,cgexj,jern1,jern2,jetmax,jgext,net,
     $    netmax,nstmax,xbarlg)
c
          do ns = nr1,nr2
            actlg(ns) = xbarlg(ns) + acflg(ns)
          enddo
c
          affpr = 0.
          do je = 1,jgext(ne)
            ns = jern1(je,ne) - 1
            do ie = 1,ngext(je,ne)
              ns = ns + 1
              xx = xbar(ns)
              if (xx .gt. 0.) then
                call afcalc(actlg,af,afcnst,cdrs,jflag,jsflag,ndrs,
     $          ndrsmx,ndrsr,ns,nstmax,si,xlks)
                if (af .gt. -9999999.) then
                  affpr = affpr + xx*af
                else
                  affpr = -9999999.
                  go to 110
                endif
              endif
            enddo
          enddo
  110     continue
          afrc1(nrc) = -affpr
c
          jpflag(np) = jpflgx
c
          do je = 1,jgext(ne)
            ns = jern1(je,ne) - 1
            do ie = 1,ngext(je,ne)
              ns = ns + 1
              jsflag(ns) = jsflgx(ie,je)
              jflag(ns) = jflgx(ie,je)
              acflg(ns) = acflge(ie,je)
              actlg(ns) = actlge(ie,je)
              xbarlg(ns) = xbarle(ie,je)
              xbar(ns) = xbare(ie,je)
            enddo
          enddo
        endif
c
      enddo
c
      end
