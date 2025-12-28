      subroutine gpheh(acflg,actlg,actwlg,adh,ah,ahmes,ahnbs,conc,
     $ eh,ehfac,ehmes,ehnbs,farad,fo2lg,fxi,iopg,mrmlra,nchlor,nhydr,
     $ nopgmx,noutpt,nstmax,nttyo,pch,pe,pemes,penbs,ph,phcl,phmes,
     $ phnbs,qphcl,qredox,qrho,xlke)
c
c     This subroutine computes the pH, redox potential, and pe
c     (electron activity function) on the operational pH scale used
c     in the calculations (see iopg(2)), the NBS pH scale,
c     and the scale on which pH is numerically equal to -log m(H+).
c     It also computes the pHCl, which is equivalent to pH + pCl.
c
c     This subroutine is called by:
c
c       EQ3NR/scripx.f
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       acflg  = array of logarithms of activity coefficients
c       actlg  = array of logarithms of aqueous species activities;
c                  actlg(nchlor) is the log activity of Cl-
c       actwlg = the log activity of water
c       adh    = the Debye-Huckel A(gamma) parameter
c       conc   = array of species concentrations
c       ehfac  = the Eh factor
c       farad  = the Faraday constant
c       iopg   = array of activity coefficient option switches;
c                  iopg(2) defines the operational pH scale used in
c                  the computation:
c                    -1 = Internal scale corresponding to the activity
c                           coefficient model defined by iopg(1)
c                     0 = The NBS pH scale (log gamma Cl- is defined
c                           by the Bates-Guggenheim expression)
c                     1 = The Mesmer pH scale (log gamma H+ is defined
c                           as zero)
c       fo2lg  = log oxygen fugacity
c       nchlor = species index of the chloride ion
c       nhydr  = species index of the hydrogen ion
c       qredox = logical flag, = .true. if there is oxidation-reduction
c                  in the modeled system
c       qrho   = .true. if a solution density value is available
c       fxi    = the ionic strength
c       xlke   = log K for special reaction for calculating Eh
c                 from log oxygen fugacity
c
c     Principal output:
c
c       ah     = Ah corresponding to the pH used in the computation
c       ahmes  = Ah corresponding to the Mesmer pH scale
c       ahnbs  = Ah corresponding to the NBS pH scale
c       eh     = Eh corresponding to the pH used in the computation
c       ehmes  = Eh corresponding to the Mesmer pH scale
c       ehnbs  = Eh corresponding to the NBS pH scale
c       ph     = pH on scale used in computation
c       phcl   = pHCl
c       phmes  = pH on the Mesmer scale
c       phnbs  = pH on the NBS scale
c       pe     = pe function corresonding to Eh
c       pemes  = pe function corresonding to the Mesmer pH scale
c       penbs  = pe function corresonding to the NBS pH scale
c       qhcl   = logical flag, = .true. if pHCl is defined
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nopgmx,nstmax
c
      integer noutpt,nttyo
c
      integer iopg(nopgmx)
      integer nchlor,nhydr
c
      logical qphcl,qredox,qrho
c
      real(8) tlg
c
      real(8) acflg(nstmax),actlg(nstmax),conc(nstmax)
      real(8) actwlg,adh,ah,ahmes,ahnbs,eh,ehfac,ehmes,ehnbs,farad,
     $ fo2lg,fxi,mrmlra,pch,pe,pemes,penbs,ph,phcl,phmes,phnbs,xlke
c
c----------------------------------------------------------------------
c
c     Local variable declarations.
c
      real*8 acfnbs,ahfx,delacf,ehfx,eterm,eterm0,eterm1
c
c----------------------------------------------------------------------
c
c     Compute pH, Eh, and pe for the operational pH scale.
c
      ph = -actlg(nhydr)
      eh = -99999.
      pe = -99999.
      ah = -99999.
      ehfx = 0.25*ehfac
      ahfx = 0.001*farad
      if (qredox) then
        eterm = fo2lg - 4.*ph - 2.*actwlg - xlke
        eh = ehfx*eterm
        pe = eh/ehfac
        ah = ahfx*eh
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute pH, Eh, and pe consistent with the NBS pH scale.
c     Note that this requires the activity coefficient of the chloride
c     ion on the operational scale.
c
      if (iopg(2) .eq. 0)  then
        phnbs = ph
        ehnbs = eh
        penbs = pe
        ahnbs = ah
      else
        phnbs = -99999.
        ehnbs = -99999.
        penbs = -99999.
        ahnbs = -99999.
        if (nchlor. gt. 0) then
          call nbsgam(acfnbs,adh,fxi,nchlor,noutpt,nttyo)
          delacf = acfnbs - acflg(nchlor)
          phnbs = ph + delacf
          if (qredox) then
            eterm0 = eterm + 4.*(ph - phnbs)
            ehnbs = ehfx*eterm0
            penbs = ehnbs/ehfac
            ahnbs = ahfx*ehnbs
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute pH, Eh, and pe consistent with the Mesmer pH scale.
c     Mesmer pH is the same as pmH.
c
      ehmes = -99999.
      pemes = -99999.
      if (iopg(2) .eq. 1) then
        phmes = ph
        ehmes = eh
        pemes = pe
        ahmes = ah
      else
        delacf = - acflg(nhydr)
        phmes = ph - delacf
        if (qredox) then
          eterm1 = eterm + 4.*(ph - phmes)
          ehmes = ehfx*eterm1
          pemes = ehmes/ehfac
          ahmes = ahfx*ehmes
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute pcH.
c
      pch = 99999.
      if (qrho) then
        pch = phmes - tlg(mrmlra)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute pHCl = pH + pCl
c
      qphcl = .false.
      phcl = -99999.
      if (nchlor .gt. 0) then
        if (conc(nchlor) .gt. 0.) then
          qphcl = .true.
          phcl = ph - actlg(nchlor)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
