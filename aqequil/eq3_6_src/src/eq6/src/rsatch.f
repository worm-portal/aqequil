      subroutine rsatch(csts,egers,egexs,iern1,ietmax,iindx1,iktmax,
     $ iopt,ipndx1,jcode,jern1,jern2,jetmax,jgext,jpflag,jreac,kmax,
     $ km1,kmt,kx1,kxt,loph,losp,moph,morr,mosp,mrgers,mtb,mtb0,
     $ nbaspd,nbtmax,ncmpr,nern1,nern2,nert,nertmx,netmax,ngext,
     $ noptmx,noutpt,nptmax,nrct,nrctmx,nrk,nrndex,nstmax,nsts,nstsmx,
     $ nstsr,nttyo,nxridx,nxrt,nxrtmx,qreq,rxbar,tolxsf,uphase,ureac,
     $ uspec,xbar,xbarlg,zvclg1,zvec1)
c
c     This subroutine tests reactants for saturation.
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
      integer ietmax,iktmax,jetmax,kmax,nbtmax,nertmx,netmax,noptmx,
     $ nptmax,nrctmx,nstmax,nstsmx,nxrtmx
c
      integer noutpt,nttyo
c
      integer iindx1(kmax),iopt(noptmx),ipndx1(kmax),jcode(nrctmx),
     $ jern1(jetmax,netmax),jern2(jetmax,netmax),jgext(netmax),
     $ jpflag(nptmax),jreac(nrctmx),nbaspd(nbtmax),ncmpr(2,nptmax),
     $ ngext(jetmax,netmax),nrk(2,nrctmx),nrndex(nrctmx),nsts(nstsmx),
     $ nstsr(2,nstmax),nxridx(nrctmx)
c
      integer iern1,km1,kmt,kx1,kxt,nern1,nern2,nert,nrct,nxrt
c
      logical qreq
c
      character*48 uspec(nstmax)
      character*24 uphase(nptmax),ureac(nrctmx)
c
      real*8 csts(nstsmx),egexs(ietmax,jetmax,netmax),loph(nptmax),
     $ egers(ietmax,jetmax,netmax),losp(nstmax),moph(nptmax),
     $ morr(nrctmx),mosp(nstmax),mrgers(ietmax,jetmax,nertmx),
     $ mtb(nbtmax),mtb0(nbtmax),rxbar(iktmax,nxrtmx),xbar(nstmax),
     $ xbarlg(nstmax),zvclg1(kmax),zvec1(kmax)
c
      real*8 tolxsf
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ie,ik,je,jlen,j2,kcol,kount,n,nb,ne,ner,np,npp,nrc,nrn1,
     $ nrn2,nr1,nr2,ns,nss,nt,nxr
c
      integer ilnobl
c
      character*56 uspn56
c
      real*8 dce,dm,dmx,dx,ecomp,ed,lx,lxx,mx,mxx,xcomp,xd
c
      real*8 texp,tlg
c
c-----------------------------------------------------------------------
c
      qreq = .false.
c
c     Pure minerals.
c
      do nrc = 1,nrct
        if (jcode(nrc).eq.0 .and. nrk(2,nrc).eq.0) then
          if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
            np = nrndex(nrc)
            if (jpflag(np) .gt. 0) go to 230
            jreac(nrc) = 0
            if (jpflag(np) .ne. -1) go to 230
            jreac(nrc) = -1
            if (iopt(1) .eq. 0) then
              jreac(nrc) = 2
              qreq = .true.
              j2 = ilnobl(ureac(nrc))
              write (noutpt,1000) ureac(nrc)(1:j2)
 1000         format(/' Reactant ',a,' has become saturated and any',
     $        ' remaining',/' unreacted mass has been transferred to',
     $        ' the ES.',/)
              dm = morr(nrc)
              morr(nrc) = 0.
              mx = moph(np) + dm
              lx = tlg(mx)
              moph(np) = mx
              loph(np) = lx
              ns = ncmpr(1,np)
              mosp(ns) = mx
              losp(ns) = lx
              do kcol = km1,kmt
                if (ns .eq. iindx1(kcol)) then
                  zvclg1(kcol) = losp(ns)
                  zvec1(kcol) = mx
                  nr1 = nstsr(1,ns)
                  nr2 = nstsr(2,ns)
                  do n = nr1,nr2
                    nb = nsts(n)
                    dce = csts(n)*dm
                    mtb(nb) = mtb(nb) + dce
                    mtb0(nb) = mtb0(nb) + dce
                    enddo
                  go to 230
                endif
              enddo
              j2 = ilnobl(ureac(nrc))
              write (noutpt,1010) ureac(nrc)(1:j2)
              write (nttyo,1010) ureac(nrc)(1:j2)
 1010         format(/' * Error - (EQ6/rsatch) The reactant ',a,
     $        'is saturated,',/7x,"but it isn't in the ES.")
              stop
            endif
          endif
        endif
  230   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Solid solutions.
c
      if (nxrt .le. 0) go to 340
c
      do nrc = 1,nrct
        if (jcode(nrc).eq.1 .and. nrk(2,nrc).eq.0) then
            if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
            np = nrndex(nrc)
            if (jpflag(np) .ge. 2) go to 330
            jreac(nrc) = 0
            if (jpflag(np) .ne. -1) go to 330
            jreac(nrc) = -1
            if (iopt(1) .ne. 0) go to 330
c
            xcomp = 0.
            nrn1 = ncmpr(1,np)
            nrn2 = ncmpr(2,np)
            nt = nrn2 - nrn1 + 1
            nxr = nxridx(nrc)
            ik = 0
              do ns = nrn1,nrn2
              ik = ik + 1
              xd = xbar(ns) - rxbar(ik,nxr)
              xcomp = xd*xd + xcomp
            enddo
            xcomp = sqrt(xcomp/nt)
            if (xcomp .gt. tolxsf) then
              j2 = ilnobl(ureac(nrc))
              write (noutpt,1020) uphase(np)(1:j2)
              write (nttyo,1020) uphase(np)(1:j2)
 1020         format(/' * Note - (EQ6/rsatch) The solid solution',
     $        ' reactant ',a,/7x,'has become saturated, but it',
     $        ' differs in composition from the',/7x,'corresponding',
     $        " product phase. It can't be moved into the ES.")
              go to 330
            endif
c
            jreac(nrc) = 2
            qreq = .true.
            write (noutpt,1000) ureac(nrc)
            dm = morr(nrc)
            morr(nrc) = 0.
            mx = moph(np) + dm
            moph(np) = mx
            loph(np) = tlg(mx)
c
            kount = 0
            do kcol = kx1,kxt
              ns = iindx1(kcol)
              npp = ipndx1(kcol)
              if (npp .eq. np) then
                kount = kount + 1
                lxx = loph(np) + xbarlg(ns)
                losp(ns) = lxx
                mxx = texp(lxx)
                mosp(ns) = mxx
                zvclg1(kcol) = lxx
                zvec1(kcol) = texp(lxx)
                dmx = dm*xbar(ns)
                nr1 = nstsr(1,ns)
                nr2 = nstsr(2,ns)
                do n = nr1,nr2
                  nb = nsts(n)
                  dce = csts(n)*dmx
                  mtb(nb) = mtb(nb) + dce
                  mtb0(nb) = mtb0(nb) + dce
                enddo
              endif
            enddo
c
            if (kount .le. 0) then
              j2 = ilnobl(ureac(nrc))
              write (noutpt,1010) ureac(nrc)(1:j2)
              write (nttyo,1010) ureac(nrc)(1:j2)
              stop
            endif
c
          endif
        endif
  330   continue
      enddo
c
  340 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Generic ion exchangers.
c
      if (nert .le. 0) go to 440
c
      do nrc = 1,nrct
        if (jcode(nrc).eq.5 .and. nrk(2,nrc).eq.0) then
            if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
            np = nrndex(nrc)
            if (jpflag(np) .ge. 2) go to 430
            jreac(nrc) = 0
            if (jpflag(np) .ne. -1) go to 430
            jreac(nrc) = -1
c
            ne = np - iern1 + 1
            ner = nxridx(nrc)
c
            ecomp = 0.
            do je = 1,jgext(ne)
              do ie = 1,ngext(je,ne)
                ed = egexs(ie,je,ne) - egers(ie,je,ner)
                ecomp = ed*ed + ecomp
              enddo
            enddo
            ecomp = sqrt(ecomp/nt)
            if (ecomp .gt. tolxsf) then
              j2 = ilnobl(ureac(nrc))
              write (noutpt,1040) uphase(np)(1:j2)
              write (nttyo,1040) uphase(np)(1:j2)
 1040         format(/' * Note - (EQ6/rsatch) The generic ion',
     $        ' exchanger reactant'/7x,a,' has become saturated, but',
     $        ' it differs in composition',/7x,'from the corresponding',
     $        " product phase. It can't be moved into the ES.")
              go to 440
            endif
c
            jreac(nrc) = 2
            qreq = .true.
            write (noutpt,1000) ureac(nrc)
            dm = morr(nrc)
            morr(nrc) = 0.
            mx = moph(np) + dm
            moph(np) = mx
            loph(np) = tlg(mx)
c
            kount = 0
c
c           Create matrix index range markers for exchange species.
c
            do kcol = 2,km1 - 1
              npp = ipndx1(kcol)
              if (npp .eq. np) then
                kount = kount + 1
                ns = iindx1(kcol)
                do je = 1,jgext(ne)
                  if (ns.ge.jern1(je,ne) .and. ns.le.jern2(je,ne))
     $            go to 410
                enddo
c
c               Calling sequence substitutions:
c                 uspec(ns) for unam48
c
                call fmspnm(jlen,uspec(ns),uspn56)
                write (noutpt,1050) uspn56(1:jlen)
                write (nttyo,1050) uspn56(1:jlen)
 1050           format(/' * Error - (EQ6/rsatch) Programming error',
     $          " trap: Can't determine",/7x,'the exchange site index',
     $          ' je for the species',/7x,a,'.')
                stop
c
  410           do ie = 1,ngext(je,ne)
                  nss = ie + jern1(je,ne) -1
                  if (ns .eq. nss) go to 420
                enddo
c
c               Calling sequence substitutions:
c                 uspec(ns) for unam48
c
                call fmspnm(jlen,uspec(ns),uspn56)
                write (noutpt,1060) uspn56(1:jlen)
                write (nttyo,1060) uspn56(1:jlen)
 1060           format(/' * Error - (EQ6/rsatch) Programming error',
     $          " trap: Can't determine",/7x,'the exchange species',
     $          ' index ie of the species',/7x,a,'.')
                stop
c
  420           dx = dm*mrgers(ie,je,ner)
                mxx = mosp(ns) + dx
                mosp(ns) = mxx
                lxx = tlg(mxx)
                losp(ns) = lxx
                zvclg1(kcol) = lxx
                zvec1(kcol) = mxx
                nr1 = nstsr(1,ns)
                nr2 = nstsr(2,ns)
                if (je .le. 1) then
c
c                 Have the first site. Increment the mass balance
c                 totals in a straightforward manner.
c
                  do n = nr1,nr2
                    nb = nsts(n)
                    dce = dx*csts(n)
                    mtb(nb) = mtb(nb) + dce
                    mtb0(nb) = mtb0(nb) + dce
                  enddo
                else
c
c                 Have a site beyond the first. Increment the mass
c                 balance totals in the usual manner, except do not
c                 increment here the total for the exchanger
c                 substrate. The complete increment for the substrate
c                 is obtained by considering only one site.
c
                  do n = nr1,nr2
                    nb = nsts(n)
                    ns = nbaspd(nb)
                    if (ns.lt.nern1 .or. ns.gt.nern2) then
                      dce = dx*csts(n)
                      mtb(nb) = mtb(nb) + dce
                      mtb0(nb) = mtb0(nb) + dce
                    endif
                  enddo
                endif
              endif
            enddo
c
            if (kount .le. 0) then
              j2 = ilnobl(ureac(nrc))
              write (noutpt,1010) ureac(nrc)(1:j2)
              write (nttyo,1010) ureac(nrc)(1:j2)
              stop
            endif
c
          endif
        endif
  430   continue
      enddo
c
  440 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
