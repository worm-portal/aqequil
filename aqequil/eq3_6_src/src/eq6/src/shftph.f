      subroutine shftph(emop,emop0,emos,emos0,fdpe0,fdpem1,
     $ fdse0,fdsem1,iemop,iemos,iern1,iern2,ietmax,iindx1,
     $ imrn1,imrn2,ipndx1,ixrn1,ixrn2,jern1,jetmax,jgext,
     $ jpflag,jsflag,kbt,kmax,km1,kmt,kx1,kxt,loph,losp,moph,
     $ mosp,mprph,mprsp,mrgexs,nbtmax,ncmpe,ncmpr,netmax,
     $ ngext,nordmx,noutpt,np,npet,npetmx,nptmax,nsetmx,
     $ nstmax,nttyo,qshftd,qtotsh,uphase,xbar,xbarlg,zklgmn,
     $ zklogl,zvclg0,zvclg1,zvec0,zvec1)
c
c     This subroutine shifts mass of a single phase from the equilibrium
c     system (ES) to the physically removed subsystem (PRS), which is
c     a system conceptually out of contact with the aqueous solution.
c     Here np is the index of the desired phase.
c
c     Presently, the following kinds of phases may be shifted:
c
c       Generic ion exchange phases
c       Pure minerals
c       Solid solutions
c
c     A shift may be total or partial. If qtotsh = .true., the shift
c     is total and the entire mineral phase is transferred. Otherwise,
c     the shift is partial and the shift transfers only some of the
c     mass. At least 10**zklgmn moles of mass will remain in the ES.
c     If the mass is less than that to start with, no shift takes
c     place. Otherwise, the ES mass is reduced by multiplying the
c     existing mass by a factor of 10**-zklogl.
c
c     This subroutine decrements the ES mass of the phase and
c     correspondingly increments its PRS mass. However, it does
c     not recalculate the component total masses for the ES. After
c     this subroutine has been called, once or in a sequence, a call
c     must be made to EQ6/escalc.f, which recalculates the component
c     mass balance totals.
c
c     If a total shift is made, it will subsequently be necessary to
c     update the iemop indexing array and associated arrays used in
c     tracking the masses of phases in the ES.
c
c     If a shift is made, qshftd is returned as .true.; otherwise,
c     it is returned as .false.
c
c     The way the shift is made presumes that the calculation is
c     currently at the "zero" point, the value of Xi from which the
c     calculation will then step out.
c
c     This subroutine is called by:
c
c       EQ6/dumpdp.f
c       EQ6/pshfta.f
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
      integer ietmax,jetmax,kmax,nbtmax,netmax,nordmx,npetmx,nptmax,
     $ nsetmx,nstmax
c
      integer noutpt,nttyo
c
      integer iemop(npetmx),iemos(nsetmx),iindx1(kmax),ipndx1(kmax),
     $ jern1(jetmax,netmax),jgext(netmax),jpflag(nptmax),jsflag(nstmax),
     $ ncmpe(2,npetmx),ncmpr(2,nptmax),
     $ ngext(jetmax,netmax)
c
      integer iern1,iern2,imrn1,imrn2,ixrn1,ixrn2,kbt,km1,kmt,
     $ kx1,kxt,np,npet
c
      logical qshftd,qtotsh
c
      character*24 uphase(nptmax)
c
      real*8 emop(npetmx),emop0(npetmx),emos(nsetmx),emos0(nsetmx),
     $ fdpe0(nordmx,npetmx),fdpem1(nordmx,npetmx),
     $ fdse0(nordmx,nsetmx),fdsem1(nordmx,nsetmx),
     $ loph(nptmax),losp(nstmax),mrgexs(ietmax,jetmax,netmax),
     $ moph(nptmax),mosp(nstmax),mprph(nptmax),mprsp(nstmax),
     $ xbar(nstmax),xbarlg(nstmax),zvclg0(kmax),zvclg1(kmax),
     $ zvec0(kmax),zvec1(kmax)
c
      real*8 zklgmn,zklogl
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ie,iptype,je,j2,kcol,kstart,kend,n,ne,npe,nr1,nr2,ns,nse
c
      integer ilnobl
c
      real*8 delmp,delms,lgexpl,lnew,lnew1,lold,lrx,mgexpl,mnew,mnew1,
     $ mold,mrx,zx
c
      real*8 texp,tlg
c
c-----------------------------------------------------------------------
c
      qshftd = .false.
c
c     Set the minimum mass of an exchanger phase that must remain in
c     the ES.
c
      lgexpl = -14.0
      mgexpl = texp(lgexpl)
c
c     Get the type of phase.
c
      iptype = 0
      if (np.ge.imrn1 .and. np.le.imrn2) then
c
c       Have a pure mineral.
c
        iptype = 1
      elseif (np.ge.ixrn1 .and. np.le.ixrn2) then
c
c       Have a solid solution.
c
        iptype = 2
      elseif (np.ge.iern1 .and. np.le.iern2) then
c
c       Have a generic ion exchange phase.
c
        iptype = 3
      endif
c
      if (iptype .le. 0) then
        j2 = ilnobl(uphase(np))
        write (noutpt,1000) np,uphase(np)(1:j2)
        write (nttyo,1000) np,uphase(np)(1:j2)
 1000   format(/' * Error - (EQ6/shftph) Programming error trap: The',
     $  ' phase',/7x,a," is not of a recognized type and can't have",
     $  ' mass shifted',/7x,'from the ES to the PRS.')
        stop
      endif
c
c     Save the current numbers of moles values.
c
      lold = loph(np)
      mold = moph(np)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the shift increment.
c
      if (qtotsh) then
c
c       Do a total shift.
c
        if (iptype .eq. 3) then
c
c         Have a generic ion exchanger phase.
c
cXX       Presently a generic ion exchange phase can not be
cXX       totally shifted.
c
          if (mold .le. mgexpl) then
            write (noutpt,1010) uphase(np)(1:j2),mold
 1010       format(3x,'Not shifted: ',a,/15x,'Mass= ',1pe12.5,' moles')
            go to 999
          endif
c
          lnew = lgexpl
          mnew = mgexpl
          delmp = mold - mgexpl
          jpflag(np) = 0
          j2 = ilnobl(uphase(np))
          write (noutpt,1020) uphase(np)(1:j2),mold,mnew
 1020     format(7x,'Shifted: ',a,/
     $    15x,'Old mass= ',1pe12.5,' moles, new mass= ',1pe12.5,
     $    ' moles')
        else
c
c         Have a pure mineral or a solid solution.
c
          lnew = -99999.
          mnew = 0.
          delmp = mold
          jpflag(np) = 0
          j2 = ilnobl(uphase(np))
          write (noutpt,1020) uphase(np)(1:j2),mold,mnew
        endif
      else
c
c       Do a partial shift.
c
        if (lold .le. zklgmn) then
c
c         Do not shift if the mass is already small.
c
          j2 = ilnobl(uphase(np))
          write (noutpt,1010) uphase(np)(1:j2),mold
          go to 999
        endif
c
        zx = lold - zklogl
        lnew = max(zx,zklgmn)
        mnew = texp(lnew)
        delmp = mold - mnew
        j2 = ilnobl(uphase(np))
        write (noutpt,1020) uphase(np)(1:j2),mold,mnew
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Complete the shift calculations.
c
      qshftd = .true.
c
c     Reset the primary mass variables.
c
      loph(np) = lnew
      moph(np) = mnew
c
c     Increment the PRS mass.
c
      mprph(np) = mprph(np) + delmp
c
      if (iptype.eq.1 .or. iptype.eq.2) then
c
c       Have a pure mineral or a solid solution.
c
        nr1 = ncmpr(1,np)
        nr2 = ncmpr(2,np)
        do ns = nr1,nr2
          delms = delmp*xbar(ns)
          mprsp(ns) = mprsp(ns) + delms
          if (qtotsh) then
c
c           Do a total shift.
c
            lnew1 = -99999.
            jsflag(ns) = 0
          else
c
c           Do a partial shift.
c
            lnew1 = lnew + xbarlg(ns)
          endif
          mnew1 = texp(lnew1)
          losp(ns) = lnew1
          mosp(ns) = mnew1
        enddo
      elseif (iptype .eq. 3) then
c
c       Have a generic ion exchange phase.
c
        ne = np - iern1 + 1
        do je = 1,jgext(ne)
          ns = jern1(je,ne) - 1
          do ie = 1,ngext(je,ne)
            ns = ns + 1
            mrx = mrgexs(ie,je,ne)
            lrx = tlg(mrx)
            delms = delmp*mrx
            mprsp(ns) = mprsp(ns) + delms
            if (qtotsh) then
c
c             Do a total shift of the component species.
c
cXX           Presently a generic ion exchange phase can not be
cXX           totally shifted; hence, neither can the component
cXX           species.
c
              lnew1 = lgexpl
            else
c
c             Do a partial shift.
c
              lnew1 = lnew + lrx
            endif
            mnew1 = texp(lnew1)
            losp(ns) = lnew1
            mosp(ns) = mnew1
          enddo
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Reset the variables used for tracking mole numbers. If a phase
c     is totally purged from the ES, it will subsequently be necessary
c     to update the associated index arrays and reload the
c     corresponding data.
c
      do npe = 1,npet
        if (np .eq. iemop(npe)) go to 110
      enddo
      npe = 0
  110 continue
c
      if (npe .gt. 0) then
        mnew = moph(np)
        emop(npe) = mnew
        emop0(npe) = mnew
c
c       Zero the associated finite differences as they are
c       no longer valid. The derivatives (demop0) previously
c       derived from them remain valid and are not zeroed here.
c       Zeroing the finite differences here will effectively
c       drop the order for the mole number tracking of the
c       affected phase.
c
        do n = 1,nordmx
          fdpe0(n,npe) = 0.
          fdpem1(n,npe) = 0.
        enddo
c
        nr1 = ncmpe(1,npe)
        nr2 = ncmpe(2,npe)
        do nse = nr1,nr2
          ns = iemos(nse)
          mnew1 = mosp(ns)
          emos(nse) = mnew1
          emos0(nse) = mnew1
c
c         Zero the associated finite differences as they are
c         no longer valid. The derivatives (demos0) previously
c         derived from them remain valid and are not zeroed here.
c         Zeroing the finite differences here will effectively
c         drop the order for the mole number tracking of the
c         affected species.
c
          do n = 1,nordmx
            fdse0(n,nse) = 0.
            fdsem1(n,nse) = 0.
          enddo
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find the corresponding matrix positions. As coded here, there
c     need not be any such.
c
      kstart = 0
      kend  = -1
c
      if (iptype .eq. 1) then
c
c       Have a pure mineral.
c
        do kcol = km1,kmt
          if (ipndx1(kcol) .eq. np) then
            kstart = kcol
            go to 130
          endif
        enddo
  130   if (kstart .gt. 0) kend = kstart
      elseif (iptype .eq. 2) then
c
c       Have a solid solution.
c
        do kcol = kx1,kxt
          if (ipndx1(kcol) .eq. np) then
            kstart = kcol
            go to 150
          endif
        enddo
  150   if (kstart .gt. 0) then
          do kcol = kxt,kstart,-1
            if (ipndx1(kcol) .eq. np) then
              kend = kcol
              go to 160
            endif
          enddo
  160     continue
        endif
      elseif (iptype .eq. 3) then
c
c       Have a generic ion exchange phase.
c
        do kcol = 1,kbt
          if (ipndx1(kcol) .eq. np) then
            kstart = kcol
            go to 170
          endif
        enddo
  170   if (kstart .gt. 0) then
          do kcol = kbt,kstart,-1
            if (ipndx1(kcol) .eq. np) then
              kend = kcol
              go to 180
            endif
          enddo
  180     continue
        endif
      endif
c
      do kcol = kstart,kend
        ns = iindx1(kcol)
        lnew1 = losp(ns)
        mnew1 = mosp(ns)
        zvclg0(kcol) = lnew1
        zvclg1(kcol) = lnew1
        zvec0(kcol) = mnew1
        zvec1(kcol) = mnew1
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
c
      end
