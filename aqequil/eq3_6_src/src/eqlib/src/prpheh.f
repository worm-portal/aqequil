      subroutine prpheh(ah,ahmes,ahnbs,eh,ehmes,ehnbs,iopg,
     $ nopgmx,noutpt,pch,pe,pemes,penbs,ph,phcl,phmes,phnbs,
     $ qphcl,qredox,qrho)
c
c     This subroutine prints the pH, Eh (redox potential), and pe
c     (log electron activity function) on the following pH scales:
c
c       1. The scale corresponding to the activity coefficient
c          model equations, unless these are rescaled to be
c          consistent with the NBS pH scale or the Mesmer pH scale
c          (see iopg(2), discussed below).
c
c       2. The NBS pH scale.
c
c       3. The Mesmer pH scale (pH is numerically equal to -log m(H+)).
c
c     Note the following options for the unscaled model pH scale,
c     which are chosen by the option switch iopg(2):
c
c       = -1   The internal scale is defined for the activity
c              coefficient model defined by iopg(1)
c
c       =  0   The modified NBS pH scale (log gamma Cl-
c              is defined by the Bates-Guggenheim expression)
c
c       =  2   The Mesmer pH (pmH) scale
c
c     This subroutine also prints the pHCl, which is equivalent to
c     pH + pCl and independent of the choice of pH scale.
c
c     The quantities printed by this subroutine are computed by
c     EQLIB/gpheh.f.
c
c     This subroutine is called by:
c
c       EQ3NR/scripx.f
c       EQ6/scripz.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ah     = Ah corresponding to the pH used in the computation
c       ahmes  = Ah corresponding to the phmes
c       ahnbs  = Ah corresponding to the phnbs
c       pch    = pcH, pH based on the molarity of the hydrogen ion
c       ph     = pH on scale used in computation
c       phcl   = pHCl
c       phmes  = pH on the Mesmer scale (same as the pmH)
c       phnbs  = pH on the modified NBS scale
c       eh     = Eh corresponding to the pH used in the computation
c       ehmes  = Eh corresponding to the phmes
c       ehnbs  = Eh corresponding to the phnbs
c       pe     = pe function corresonding to Eh
c       pemes  = pe function corresonding to ehmes
c       penbs  = pe function corresonding to ehnbs
c       iopg   = array of activity coefficient option switches
c       noutpt = the unit number of the output file
c       qrho   = .true. if a solution density value is available
c       qpchl  = logical flag, = .true. if pHCl is defined
c       qredox = logical flag, = .true. if there is oxidation-reduction
c                  in the modeled system
c
c     Principal output:
c
c       None
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nopgmx
c
      integer noutpt
c
      integer iopg(nopgmx)
c
      logical qphcl,qredox,qrho
c
      real(8) ah,ahmes,ahnbs,eh,ehmes,ehnbs,pch,pe,pemes,penbs,ph,
     $ phcl,phmes,phnbs
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,isc
c
      integer ilnobl
c
      character*24 uphscs(6)
c
c-----------------------------------------------------------------------
c
      data (uphscs(isc), isc = 1,6)  /'Internal pH scale       ',
     $     'NBS pH scale            ','Mesmer pH (pmH) scale   ',
     $     'Pitzer pH scale         ','B-dot pH scale          ',
     $     'Davies pH scale         '/
c
c-----------------------------------------------------------------------
c
      isc = 1
      if (iopg(2) .ge. 0) then
        isc = iopg(2) + 2
      else
        if (iopg(1) .eq. 1) then
          isc = 4
        elseif (iopg(1) .eq. 0) then
          isc = 5
        elseif (iopg(1) .eq. -1) then
          isc = 6
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qredox) then
        write (noutpt,1000)
 1000   format(/11x,'--- The pH, Eh, pe-, and Ah on various pH',
     $  ' scales ---')
        write (noutpt,1010)
 1010   format(/31x,'pH',5x,'Eh, volts',8x,'pe-',7x,'Ah, kcal',/)
c
        if (isc.ne.2 .and. isc.ne.3) then
          write (noutpt,1020) uphscs(isc),ph,eh,pe,ah
 1020     format(1x,a24,2x,f8.4,3x,f8.4,4x,1pe11.4,2x,0pf10.4)
        endif
c
        write (noutpt,1020) uphscs(2),phnbs,ehnbs,penbs,ahnbs
        write (noutpt,1020) uphscs(3),phmes,ehmes,pemes,ahmes
      else
        write (noutpt,1030)
 1030   format(/11x,'--- The pH on various pH scales ---')
        write (noutpt,1040)
 1040   format(/40x,'pH',/)
        if (isc.ne.2 .and. isc.ne.3) then
          write (noutpt,1050) uphscs(isc),ph
 1050     format(6x,a24,3x,f11.4)
        endif
c
        write (noutpt,1050) uphscs(2),phnbs
        write (noutpt,1050) uphscs(3),phmes
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1100)
 1100 format(/1x)
c
      if (qrho) then
        write (noutpt,1110) pch
 1110   format(6x,'pcH= ',f11.4)
      endif
c
      if (qphcl) then
        write (noutpt,1140) phcl
 1140   format(5x,'pHCl= ',f11.4)
      else
        write (noutpt,1150)
 1150   format(5x,'The pHCl is undefined because no Cl- is present.')
      endif
c
      write (noutpt,1160)
 1160 format(1x)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      j2 = ilnobl(uphscs(isc))
      write (noutpt,1190) uphscs(isc)(1:j2)
 1190 format(/3x,'The single ion activities and activity coefficients',
     $ ' listed below',/3x,'are consistent with the ',a,'.',/)
c
      end
