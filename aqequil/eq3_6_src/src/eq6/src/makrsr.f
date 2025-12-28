      subroutine makrsr(cbsri,cesri,cess,eps100,ibsrti,iesrti,jcode,
     $ nbt1mx,nct,nctmax,ness,nessmx,nessr,noutpt,nrct,nrctmx,nsrtmx,
     $ nstmax,nttyo,ubsri,uelem,uesri,ureac,uspec,zchar)
c
c     This subroutine makes a reaction for a special reactant. The
c     reaction is written in terms of the special reactant and strict
c     basis species only, and is charge balanced.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       cesri  = array of coefficients for chemical elements composing
c                  special reactants
c       ibsrti = array of numbers of species in reactions for special
c                  reactants
c       iesrti = array of numbers of chemical elements composing special
c                  reactants
c       jcode  = flag array denoting the types of reactants
c       nbt1mx = the maximum number of basis species plus 1
c       nct    = the number of chemical elements
c       nctmax = the maximum number of chemical elements
c       ness   = array of indices of chemical elements composing
c                  species
c       nessmx = the maximum number of entries in the cess or ness array
c       nessr  = pointer array giving the ranges in the cess and ness
c                  arrays giving the compositions of species
c       nrctmx = the maximum number of reactants
c       nsrtmx = the maximum number of special reactants
c       uelem  = array of names of chemical elements
c       uesri  = array of names of chemical elements composing special
c                  reactants
c       ureac  = array of names of reactants
c       uspec  = array of names of species
c       zchar  = array of electrical charges of species
c
c     Principal output:
c
c       cbsri  = array of coefficients for reactions for special
c                  reactants
c       ibsrti = array of numbers of species in reactions for special
c                  reactants
c       ubsri  = array of names of species in reactions for special
c                  reactants
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbt1mx,nctmax,nessmx,nrctmx,nsrtmx,nstmax
c
      integer noutpt,nttyo
c
      integer ibsrti(nsrtmx),iesrti(nsrtmx),jcode(nrctmx),
     $ ness(nessmx),nessr(2,nstmax)
c
      integer nct,nrct
c
      character*48 uspec(nstmax)
      character*24 ubsri(nbt1mx,nsrtmx),ureac(nrctmx)
      character*8 uelem(nctmax),uesri(nctmax,nsrtmx)
c
      real*8 cbsri(nbt1mx,nsrtmx),cess(nessmx),cesri(nctmax,nsrtmx),
     $ zchar(nstmax)
c
      real*8 eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3,n,nc,nch,nci,ncih,ncio,nco,nerr,nn,nrc,nr1,nr2,
     $ ns,nsr
c
      integer ilnobl
c
      real*8 ctxo,ctxh,cx,cxe,cxh,cxo,ztx
c
c-----------------------------------------------------------------------
c
      nerr = 0
c
      nco = 0
      do nc = 1,nct
        if (uelem(nc)(1:2) .eq. 'O ') then
          nco = nc
          go to 100
        endif
      enddo
  100 continue
c
      nch = 0
      do nc = 1,nct
        if (uelem(nc)(1:2) .eq. 'H ') then
          nch = nc
          go to 110
        endif
      enddo
  110 continue
c
      nsr = 0
      do nrc = 1,nrct
        if (jcode(nrc) .eq. 2) then
          nsr = nsr + 1
          if (ibsrti(nsr) .eq. 0) then
            cbsri(1,nsr) = -1.0
            ubsri(1,nsr) = ureac(nrc)
            n = 1
c
            ncio = 0
            do nci = 1,iesrti(nsr)
              if (uesri(nci,nsr)(1:2) .eq. 'O ') then
                ncio = nci
                go to 120
              endif
            enddo
  120       continue
c
            ncih = 0
            do nci = 1,iesrti(nsr)
              if (uesri(nci,nsr)(1:2) .eq. 'H ') then
                ncih = nci
                go to 130
              endif
            enddo
  130       continue
c
            ctxo = 0.
            ctxh = 0.
            if (ncio .gt. 0) ctxo = -cesri(ncio,nsr)
            if (ncih .gt. 0) ctxh = -cesri(ncih,nsr)
            ztx = 0.
c
c           Map all chemical elements except H and O.
c
            do nci = 1,iesrti(nsr)
              if (nci.ne.ncio .and. nci.ne.ncih) then
                do nc = 1,nct
                  if (uesri(nci,nsr)(1:8) .eq. uelem(nc)(1:8)) then
                    n = n + 1
                    ns = nc
                    nr1 = nessr(1,ns)
                    nr2 = nessr(2,ns)
c
                    cxe = 0.
                    cxo = 0.
                    cxh = 0.
                    do nn = nr1,nr2
                      if (ness(nn) .eq. nco) then
                        cxo = cess(nn)
                      elseif (ness(nn) .eq. nch) then
                        cxh = cess(nn)
                      else
                        cxe = cess(nn)
                      endif
                    enddo
                    if (cxe .eq. 0.) then
                      if (cxh.ne.0. .and. cxo.eq.0.) then
                        cxe = cxh
                        cxh = 0.
                      elseif (cxo.ne.0. .and. cxh.eq.0.) then
                        cxe = cxo
                        cxo = 0.
                      elseif (cxo.ne.0. .and. cxh.ne.0.) then
                        cxe = cxo
                        cxo = 0.
                      endif
                    endif
c
                    cx = cesri(nci,nsr)/cxe
                    cbsri(n,nsr) = cx
                    ubsri(n,nsr) = uspec(ns)(1:24)
                    ctxo = ctxo + cxo*cx
                    ctxh = ctxh + cxh*cx
                    ztx = ztx + zchar(ns)*cx
                    go to 140
                  endif
                enddo
                j2 = ilnobl(uesri(nci,nsr))
                j3 = ilnobl(ureac(nrc))
                write (noutpt,1000) uesri(nci,nsr)(1:j2),
     $          ureac(nrc)(1:j3)
                write (nttyo,1000) uesri(nci,nsr)(1:j2),
     $          ureac(nrc)(1:j3)
 1000           format(/" * Error - (EQ6/makrsr) Can't map the",
     $          /7x,'chemical element ',a,', which appears in the'
     $          /7x,'composition for the special reactant ',a,',',
     $          /7x,'into a corresponding basis species to use in',
     $          /7x,'composing reaction for that reactant.')
                nerr = nerr + 1
  140           continue
              endif
            enddo
c
c           Get H+ from charge balance.
c
            if (abs(ztx) .gt. eps100) then
              n = n + 1
              cx = -ztx
              cbsri(n,nsr) = cx
              ubsri(n,nsr) = 'H+ '
              ztx = ztx + cx
              ctxh = ctxh + cx
            endif
c
c           Get H2O from H balance.
c
            if (abs(ctxh) .gt. eps100) then
              n = n + 1
              cx = -0.5*ctxh
              cbsri(n,nsr) = cx
              ubsri(n,nsr) = 'H2O '
              ctxo = ctxo + cx
              ctxh = ctxh + 2.*cx
            endif
c
c           Get O2(g) from O balance.
c
            if (abs(ctxo) .gt. eps100) then
              n = n + 1
              cx = -0.5*ctxo
              cbsri(n,nsr) = cx
              ubsri(n,nsr) = 'O2(g)'
              ctxo = ctxo + 2.*cx
            endif
c
c           Check the balances on H, O, and charge.
c
            if (abs(ctxh).gt.eps100 .or. abs(ctxo).gt.eps100 .or.
     $        abs(ztx).gt.eps100) then
              j2 = ilnobl(ureac(nrc))
              write (noutpt,1110) ureac(nrc)(1:j2)
              write (nttyo,1110) ureac(nrc)(1:j2)
 1110         format(/' * Error - (EQ6/makrsr) The reaction which',
     $        /7x,'was composed for the special reactant "',a,'"',
     $        /7x,'fails to satisfy all balance conditions.')
              nerr = nerr + 1
            endif
c
            ibsrti(nsr) = n
          endif
        endif
      enddo
c
      if (nerr .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
