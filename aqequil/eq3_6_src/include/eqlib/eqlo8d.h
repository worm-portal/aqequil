c eqlo8d.h
c
c     Strings, arrays, and such for the menu-style ("D") input format
c     for INPUT file options common to EQ3NR and EQ6 in Version 8.
c     This INCLUDE file is referenced by EQ6, EQ3NR, XCON3, and XCON6.
c
c     The relevant arrays are declared in the EQLIB INCLUDE file
c     eqlo8.h.
c
c-----------------------------------------------------------------------
c
c     Nxmod alter/supress option strings:
c
        data ukxm(-1) / 'Suppress        '/
        data ukxm(0)  / 'Replace         '/
        data ukxm(1)  / 'AugmentLogK     '/
        data ukxm(2)  / 'AugmentG        '/
c
c-----------------------------------------------------------------------
c
c     Exchange model option strings for generic ion exchangers:
c
        data ugexmv(1)  / 'Vanselow                '/
        data ugexmv(2)  / 'Gapon                   '/
        data ugexmv(3)  / 'ERROR                   '/
        data ugexmv(4)  / 'ERROR                   '/
        data ugexmv(5)  / 'ERROR                   '/
        data ugexmv(6)  / 'ERROR                   '/
        data ugexmv(7)  / 'ERROR                   '/
        data ugexmv(8)  / 'ERROR                   '/
        data ugexmv(9)  / 'ERROR                   '/
        data ugexmv(10) / 'ERROR                   '/
c
c-----------------------------------------------------------------------
c
c     Units option strings for thermodynamic parameters of generic ion
c         exchangers:
c
        data uxfuni(1) / 'LogK/eq '/
        data uxfuni(2) / 'kcal/eq '/
        data uxfuni(3) / 'kJ/eq   '/
        data uxfuni(4) / 'cm3/eq  '/
c
c-----------------------------------------------------------------------
c
c     Iopt model option strings.
c
c       iopt(1):
c
          data uopttx(1) /'Physical System Model Selection'/
          data uoptcx(1) /'EQ6 only'/
c
          data uoptox(1,1) /'Closed system'/
          data ioptox(1,1) /0/
          data uoptox(2,1) /'Titration system'/
          data ioptox(2,1) /1/
          data uoptox(3,1) /'Fluid-centered flow-through open system'/
          data ioptox(3,1) /2/
          data uoptox(4,1) /' '/
          data ioptox(4,1) /0/
          data uoptox(5,1) /' '/
          data ioptox(5,1) /0/
c
c       iopt(2):
c
          data uopttx(2) /'Kinetic Mode Selection'/
          data uoptcx(2) /'EQ6 only'/
c
          data uoptox(1,2) /'Reaction progress mode (arbitrary kinetics)
     $'/
          data ioptox(1,2) /0/
          data uoptox(2,2) /'Reaction progress/time mode (true kinetics)
     $'/
          data ioptox(2,2) /1/
          data uoptox(3,2) /' '/
          data ioptox(3,2) /0/
          data uoptox(4,2) /' '/
          data ioptox(4,2) /0/
          data uoptox(5,2) /' '/
          data ioptox(5,2) /0/
c
c       iopt(3):
c
          data uopttx(3) /'Phase Boundary Searches'/
          data uoptcx(3) /'EQ6 only'/
c
          data uoptox(1,3) /'Search for phase boundaries and constrain t
     $he step size to match'/
          data ioptox(1,3) /0/
          data uoptox(2,3) /"Search for phase boundaries and print their
     $ locations"/
          data ioptox(2,3) /1/
          data uoptox(3,3) /"Don't search for phase boundaries"/
          data ioptox(3,3) /2/
          data uoptox(4,3) /' '/
          data ioptox(4,3) /0/
          data uoptox(5,3) /' '/
          data ioptox(5,3) /0/
c
c       iopt(4):
c
          data uopttx(4) /'Solid Solutions'/
          data uoptcx(4) /'EQ3NR, EQ6'/
c
          data uoptox(1,4) /'Ignore'/
          data ioptox(1,4) /0/
          data uoptox(2,4) /'Permit'/
          data ioptox(2,4) /1/
          data uoptox(3,4) /' '/
          data ioptox(3,4) /0/
          data uoptox(4,4) /' '/
          data ioptox(4,4) /0/
          data uoptox(5,4) /' '/
          data ioptox(5,4) /0/
C
c       iopt(5):
c
          data uopttx(5) /'Clear the ES Solids Read from the INPUT File'
     $/
          data uoptcx(5) /'EQ6 only'/
c
          data uoptox(1,5) /"Don't do it"/
          data ioptox(1,5) /0/
          data uoptox(2,5) /'Do it'/
          data ioptox(2,5) /1/
          data uoptox(3,5) /' '/
          data ioptox(3,5) /0/
          data uoptox(4,5) /' '/
          data ioptox(4,5) /0/
          data uoptox(5,5) /' '/
          data ioptox(5,5) /0/
c
c       iopt(6):
c
          data uopttx(6) /'Clear the ES Solids at the Initial Value of R
     $eaction Progress'/
          data uoptcx(6) /'EQ6 only'/
c
          data uoptox(1,6) /"Don't do it"/
          data ioptox(1,6) /0/
          data uoptox(2,6) /'Do it'/
          data ioptox(2,6) /1/
          data uoptox(3,6) /' '/
          data ioptox(3,6) /0/
          data uoptox(4,6) /' '/
          data ioptox(4,6) /0/
          data uoptox(5,6) /' '/
          data ioptox(5,6) /0/
c
c       iopt(7):
c
          data uopttx(7) /'Clear the ES Solids at the End of the Run'/
          data uoptcx(7) /'EQ6 only'/
c
          data uoptox(1,7) /"Don't do it"/
          data ioptox(1,7) /0/
          data uoptox(2,7) /'Do it'/
          data ioptox(2,7) /1/
          data uoptox(3,7) /' '/
          data ioptox(3,7) /0/
          data uoptox(4,7) /' '/
          data ioptox(4,7) /0/
          data uoptox(5,7) /' '/
          data ioptox(5,7) /0/
c
c       iopt(8):
c
          data uopttx(8) /'Not used'/
          data uoptcx(8) /'None'/
c
          data uoptox(1,8) /' '/
          data ioptox(1,8) /0/
          data uoptox(2,8) /' '/
          data ioptox(2,8) /0/
          data uoptox(3,8) /' '/
          data ioptox(3,8) /0/
          data uoptox(4,8) /' '/
          data ioptox(4,8) /0/
          data uoptox(5,8) /' '/
          data ioptox(5,8) /0/
c
c       iopt(9):
c
          data uopttx(9) /'Clear the PRS Solids Read from the INPUT file
     $'/
          data uoptcx(9) /'EQ6 only'/
c
          data uoptox(1,9) /"Don't do it"/
          data ioptox(1,9) /0/
          data uoptox(2,9) /'Do it'/
          data ioptox(2,9) /1/
          data uoptox(3,9) /' '/
          data ioptox(3,9) /0/
          data uoptox(4,9) /' '/
          data ioptox(4,9) /0/
          data uoptox(5,9) /' '/
          data ioptox(5,9) /0/
c
c       iopt(10):
c
          data uopttx(10) /'Clear the PRS Solids at the End of the Run'/
          data uoptcx(10) /'EQ6 only'/
c
          data uoptox(1,10) /"Don't do it"/
          data ioptox(1,10) /0/
          data uoptox(2,10) /'Do it, unless numerical problems cause ear
     $ly termination'/
          data ioptox(2,10) /1/
          data uoptox(3,10) /' '/
          data ioptox(3,10) /0/
          data uoptox(4,10) /' '/
          data ioptox(4,10) /0/
          data uoptox(5,10) /' '/
          data ioptox(5,10) /0/
c
c       iopt(11):
c
          data uopttx(11) /'Auto Basis Switching in pre-N-R Optimization
     $'/
          data uoptcx(11) /'EQ3NR, EQ6'/
c
          data uoptox(1,11) /'Turn off'/
          data ioptox(1,11) /0/
          data uoptox(2,11) /'Turn on'/
          data ioptox(2,11) /1/
          data uoptox(3,11) /' '/
          data ioptox(3,11) /0/
          data uoptox(4,11) /' '/
          data ioptox(4,11) /0/
          data uoptox(5,11) /' '/
          data ioptox(5,11) /0/
c
c       iopt(12):
c
          data uopttx(12) /'Auto Basis Switching after Newton-Raphson It
     $eration'/
          data uoptcx(12) /'EQ6 only'/
c
          data uoptox(1,12) /'Turn off'/
          data ioptox(1,12) /0/
          data uoptox(2,12) /'Turn on'/
          data ioptox(2,12) /1/
          data uoptox(3,12) /' '/
          data ioptox(3,12) /0/
          data uoptox(4,12) /' '/
          data ioptox(4,12) /0/
          data uoptox(5,12) /' '/
          data ioptox(5,12) /0/
c
c       iopt(13):
c
          data uopttx(13) /'Calculational Mode Selection'/
          data uoptcx(13) /'EQ6 only'/
c
          data uoptox(1,13) /'Normal path tracing'/
          data ioptox(1,13) /0/
          data uoptox(2,13) /'Economy mode (if permissible)'/
          data ioptox(2,13) /1/
          data uoptox(3,13) /'Super economy mode (if permissible)'/
          data ioptox(3,13) /2/
          data uoptox(4,13) /' '/
          data ioptox(4,13) /0/
          data uoptox(5,13) /' '/
          data ioptox(5,13) /0/
c
c       iopt(14):
c
          data uopttx(14) /'ODE Integrator Corrector Mode Selection'/
          data uoptcx(14) /'EQ6 only'/
c
          data uoptox(1,14) /'Allow Stiff and Simple Correctors'/
          data ioptox(1,14) /0/
          data uoptox(2,14) /'Allow Only the Simple Corrector'/
          data ioptox(2,14) /1/
          data uoptox(3,14) /'Allow Only the Stiff Corrector'/
          data ioptox(3,14) /2/
          data uoptox(4,14) /'Allow No Correctors'/
          data ioptox(4,14) /3/
          data uoptox(5,14) /' '/
          data ioptox(5,14) /0/
c
c       iopt(15):
c
          data uopttx(15) /'Force the Suppression of All Redox Reactions
     $'/
          data uoptcx(15) /'EQ6 only'/
c
          data uoptox(1,15) /"Don't do it"/
          data ioptox(1,15) /0/
          data uoptox(2,15) /'Do it'/
          data ioptox(2,15) /1/
          data uoptox(3,15) /' '/
          data ioptox(3,15) /0/
          data uoptox(4,15) /' '/
          data ioptox(4,15) /0/
          data uoptox(5,15) /' '/
          data ioptox(5,15) /0/
c
c       iopt(16):
c
          data uopttx(16) /'BACKUP File Options'/
          data uoptcx(16) /'EQ6 only'/
c
          data uoptox(1,16) /"Don't write a BACKUP file"/
          data ioptox(1,16) /-1/
          data uoptox(2,16) /'Write BACKUP files'/
          data ioptox(2,16) /0/
          data uoptox(3,16) /'Write a sequential BACKUP file'/
          data ioptox(3,16) /1/
          data uoptox(4,16) /' '/
          data ioptox(4,16) /0/
          data uoptox(5,16) /' '/
          data ioptox(5,16) /0/
c
c       iopt(17):
c
          data uopttx(17) /'PICKUP File Options'/
          data uoptcx(17) /'EQ3NR, EQ6'/
c
          data uoptox(1,17) /"Don't write a PICKUP file"/
          data ioptox(1,17) /-1/
          data uoptox(2,17) /'Write a PICKUP file'/
          data ioptox(2,17) /0/
          data uoptox(3,17) /' '/
          data ioptox(3,17) /0/
          data uoptox(4,17) /' '/
          data ioptox(4,17) /0/
          data uoptox(5,17) /' '/
          data ioptox(5,17) /0/
c
c       iopt(18):
c
          data uopttx(18) /'TAB File Options'/
          data uoptcx(18) /'EQ6 only'/
c
          data uoptox(1,18) /"Don't write a TAB file"/
          data ioptox(1,18) /-1/
          data uoptox(2,18) /'Write a TAB file'/
          data ioptox(2,18) /0/
          data uoptox(3,18) /'Write a TAB file, prepending TABX file dat
     $a from a previous run'/
          data ioptox(3,18) /1/
          data uoptox(4,18) /' '/
          data ioptox(4,18) /0/
          data uoptox(5,18) /' '/
          data ioptox(5,18) /0/
c
c       iopt(19):
c
          data uopttx(19) /'Advanced EQ3NR PICKUP File Options'/
          data uoptcx(19) /'EQ3NR only'/
c
          data uoptox(1,19) /'Write a normal EQ3NR PICKUP file'/
          data ioptox(1,19) /0/
          data uoptox(2,19) /'Write an EQ6 INPUT file with Quartz dissol
     $ving, relative rate law'/
          data ioptox(2,19) /1/
          data uoptox(3,19) /'Write an EQ6 INPUT file with Albite dissol
     $ving, TST rate law'/
          data ioptox(3,19) /2/
          data uoptox(4,19) /'Write an EQ6 INPUT file with Fluid 1 set u
     $p for fluid mixing'/
          data ioptox(4,19) /3/
          data uoptox(5,19) /' '/
          data ioptox(5,19) /0/
c
c       iopt(20):
c
          data uopttx(20) /'Advanced EQ6 PICKUP File Options'/
          data uoptcx(20) /'EQ6 only'/
c
          data uoptox(1,20) /'Write a normal EQ6 PICKUP file'/
          data ioptox(1,20) /0/
          data uoptox(2,20) /'Write an EQ6 INPUT file with Fluid 1 set u
     $p for fluid mixing'/
          data ioptox(2,20) /1/
          data uoptox(3,20) /' '/
          data ioptox(3,20) /0/
          data uoptox(4,20) /' '/
          data ioptox(4,20) /0/
          data uoptox(5,20) /' '/
          data ioptox(5,20) /0/
c
c-----------------------------------------------------------------------
c
c     Iopg activity coefficient option strings.
c
c       iopg(1):
c
          data uopgtx(1) /'Aqueous Species Activity Coefficient Model'/
          data uopgcx(1) /'EQ3NR, EQ6'/
c
          data uopgox(1,1) /'The Davies equation'/
          data iopgox(1,1) /-1/
          data uopgox(2,1) /'The B-dot equation'/
          data iopgox(2,1) /0/
          data uopgox(3,1) /"Pitzer's equations"/
          data iopgox(3,1) /1/
          data uopgox(4,1) /'HC + DH equations'/
          data iopgox(4,1) /2/
          data uopgox(5,1) /' '/
          data iopgox(5,1) /0/
c
c       iopg(2):
c
          data uopgtx(2) /'Choice of pH Scale (Rescales Activity Coeffic
     $ients)'/
          data uopgcx(2) /'EQ3NR, EQ6'/
c
          data uopgox(1,2) /'"Internal" pH scale (no rescaling)'/
          data iopgox(1,2) /-1/
          data uopgox(2,2) /'NBS pH scale (uses the Bates-Guggenheim equ
     $ation)'/
          data iopgox(2,2) /0/
          data uopgox(3,2) /'Mesmer pH scale (numerically, pH = -log m(H
     $+))'/
          data iopgox(3,2) /1/
          data uopgox(4,2) /' '/
          data iopgox(4,2) /0/
          data uopgox(5,2) /' '/
          data iopgox(5,2) /0/
c
c       iopg(3):
c
          data uopgtx(3) /'Not Used'/
          data uopgcx(3) /'None'/
c
          data uopgox(1,3) /' '/
          data iopgox(1,3) /0/
          data uopgox(2,3) /' '/
          data iopgox(2,3) /0/
          data uopgox(3,3) /' '/
          data iopgox(3,3) /0/
          data uopgox(4,3) /' '/
          data iopgox(4,3) /0/
          data uopgox(5,3) /' '/
          data iopgox(5,3) /0/
c
c       iopg(4):
c
          data uopgtx(4) /'Not Used'/
          data uopgcx(4) /'None'/
c
          data uopgox(1,4) /' '/
          data iopgox(1,4) /0/
          data uopgox(2,4) /' '/
          data iopgox(2,4) /0/
          data uopgox(3,4) /' '/
          data iopgox(3,4) /0/
          data uopgox(4,4) /' '/
          data iopgox(4,4) /0/
          data uopgox(5,4) /' '/
          data iopgox(5,4) /0/
c
c       iopg(5):
c
          data uopgtx(5) /'Not Used'/
          data uopgcx(5) /'None'/
c
          data uopgox(1,5) /' '/
          data iopgox(1,5) /0/
          data uopgox(2,5) /' '/
          data iopgox(2,5) /0/
          data uopgox(3,5) /' '/
          data iopgox(3,5) /0/
          data uopgox(4,5) /' '/
          data iopgox(4,5) /0/
          data uopgox(5,5) /' '/
          data iopgox(5,5) /0/
c
c       iopg(6):
c
          data uopgtx(6) /'Not Used'/
          data uopgcx(6) /'None'/
c
          data uopgox(1,6) /' '/
          data iopgox(1,6) /0/
          data uopgox(2,6) /' '/
          data iopgox(2,6) /0/
          data uopgox(3,6) /' '/
          data iopgox(3,6) /0/
          data uopgox(4,6) /' '/
          data iopgox(4,6) /0/
          data uopgox(5,6) /' '/
          data iopgox(5,6) /0/
c
c       iopg(7):
c
          data uopgtx(7) /'Not Used'/
          data uopgcx(7) /'None'/
c
          data uopgox(1,7) /' '/
          data iopgox(1,7) /0/
          data uopgox(2,7) /' '/
          data iopgox(2,7) /0/
          data uopgox(3,7) /' '/
          data iopgox(3,7) /0/
          data uopgox(4,7) /' '/
          data iopgox(4,7) /0/
          data uopgox(5,7) /' '/
          data iopgox(5,7) /0/
c
c       iopg(8):
c
          data uopgtx(8) /'Not Used'/
          data uopgcx(8) /'None'/
c
          data uopgox(1,8) /' '/
          data iopgox(1,8) /0/
          data uopgox(2,8) /' '/
          data iopgox(2,8) /0/
          data uopgox(3,8) /' '/
          data iopgox(3,8) /0/
          data uopgox(4,8) /' '/
          data iopgox(4,8) /0/
          data uopgox(5,8) /' '/
          data iopgox(5,8) /0/
c
c       iopg(9):
c
          data uopgtx(9) /'Not Used'/
          data uopgcx(9) /'None'/
c
          data uopgox(1,9) /' '/
          data iopgox(1,9) /0/
          data uopgox(2,9) /' '/
          data iopgox(2,9) /0/
          data uopgox(3,9) /' '/
          data iopgox(3,9) /0/
          data uopgox(4,9) /' '/
          data iopgox(4,9) /0/
          data uopgox(5,9) /' '/
          data iopgox(5,9) /0/
c
c       iopg(10):
c
          data uopgtx(10) /'Not Used'/
          data uopgcx(10) /'None'/
c
          data uopgox(1,10) /' '/
          data iopgox(1,10) /0/
          data uopgox(2,10) /' '/
          data iopgox(2,10) /0/
          data uopgox(3,10) /' '/
          data iopgox(3,10) /0/
          data uopgox(4,10) /' '/
          data iopgox(4,10) /0/
          data uopgox(5,10) /' '/
          data iopgox(5,10) /0/
c
c       iopg(11):
c
          data uopgtx(11) /'Not Used'/
          data uopgcx(11) /'None'/
c
          data uopgox(1,11) /' '/
          data iopgox(1,11) /0/
          data uopgox(2,11) /' '/
          data iopgox(2,11) /0/
          data uopgox(3,11) /' '/
          data iopgox(3,11) /0/
          data uopgox(4,11) /' '/
          data iopgox(4,11) /0/
          data uopgox(5,11) /' '/
          data iopgox(5,11) /0/
c
c       iopg(12):
c
          data uopgtx(12) /'Not Used'/
          data uopgcx(12) /'None'/
c
          data uopgox(1,12) /' '/
          data iopgox(1,12) /0/
          data uopgox(2,12) /' '/
          data iopgox(2,12) /0/
          data uopgox(3,12) /' '/
          data iopgox(3,12) /0/
          data uopgox(4,12) /' '/
          data iopgox(4,12) /0/
          data uopgox(5,12) /' '/
          data iopgox(5,12) /0/
c
c       iopg(13):
c
          data uopgtx(13) /'Not Used'/
          data uopgcx(13) /'None'/
c
          data uopgox(1,13) /' '/
          data iopgox(1,13) /0/
          data uopgox(2,13) /' '/
          data iopgox(2,13) /0/
          data uopgox(3,13) /' '/
          data iopgox(3,13) /0/
          data uopgox(4,13) /' '/
          data iopgox(4,13) /0/
          data uopgox(5,13) /' '/
          data iopgox(5,13) /0/
c
c       iopg(14):
c
          data uopgtx(14) /'Not Used'/
          data uopgcx(14) /'None'/
c
          data uopgox(1,14) /' '/
          data iopgox(1,14) /0/
          data uopgox(2,14) /' '/
          data iopgox(2,14) /0/
          data uopgox(3,14) /' '/
          data iopgox(3,14) /0/
          data uopgox(4,14) /' '/
          data iopgox(4,14) /0/
          data uopgox(5,14) /' '/
          data iopgox(5,14) /0/
c
c       iopg(15):
c
          data uopgtx(15) /'Not Used'/
          data uopgcx(15) /'None'/
c
          data uopgox(1,15) /' '/
          data iopgox(1,15) /0/
          data uopgox(2,15) /' '/
          data iopgox(2,15) /0/
          data uopgox(3,15) /' '/
          data iopgox(3,15) /0/
          data uopgox(4,15) /' '/
          data iopgox(4,15) /0/
          data uopgox(5,15) /' '/
          data iopgox(5,15) /0/
c
c       iopg(16):
c
          data uopgtx(16) /'Not Used'/
          data uopgcx(16) /'None'/
c
          data uopgox(1,16) /' '/
          data iopgox(1,16) /0/
          data uopgox(2,16) /' '/
          data iopgox(2,16) /0/
          data uopgox(3,16) /' '/
          data iopgox(3,16) /0/
          data uopgox(4,16) /' '/
          data iopgox(4,16) /0/
          data uopgox(5,16) /' '/
          data iopgox(5,16) /0/
c
c       iopg(17):
c
          data uopgtx(17) /'Not Used'/
          data uopgcx(17) /'None'/
c
          data uopgox(1,17) /' '/
          data iopgox(1,17) /0/
          data uopgox(2,17) /' '/
          data iopgox(2,17) /0/
          data uopgox(3,17) /' '/
          data iopgox(3,17) /0/
          data uopgox(4,17) /' '/
          data iopgox(4,17) /0/
          data uopgox(5,17) /' '/
          data iopgox(5,17) /0/
c
c       iopg(18):
c
          data uopgtx(18) /'Not Used'/
          data uopgcx(18) /'None'/
c
          data uopgox(1,18) /' '/
          data iopgox(1,18) /0/
          data uopgox(2,18) /' '/
          data iopgox(2,18) /0/
          data uopgox(3,18) /' '/
          data iopgox(3,18) /0/
          data uopgox(4,18) /' '/
          data iopgox(4,18) /0/
          data uopgox(5,18) /' '/
          data iopgox(5,18) /0/
c
c       iopg(19):
c
          data uopgtx(19) /'Not Used'/
          data uopgcx(19) /'None'/
c
          data uopgox(1,19) /' '/
          data iopgox(1,19) /0/
          data uopgox(2,19) /' '/
          data iopgox(2,19) /0/
          data uopgox(3,19) /' '/
          data iopgox(3,19) /0/
          data uopgox(4,19) /' '/
          data iopgox(4,19) /0/
          data uopgox(5,19) /' '/
          data iopgox(5,19) /0/
c
c       iopg(20):
c
          data uopgtx(20) /'Not Used'/
          data uopgcx(20) /'None'/
c
          data uopgox(1,20) /' '/
          data iopgox(1,20) /0/
          data uopgox(2,20) /' '/
          data iopgox(2,20) /0/
          data uopgox(3,20) /' '/
          data iopgox(3,20) /0/
          data uopgox(4,20) /' '/
          data iopgox(4,20) /0/
          data uopgox(5,20) /' '/
          data iopgox(5,20) /0/
c
c-----------------------------------------------------------------------
c
c     Iopr print opgion strings.
c
c       iopr(1):
c
          data uoprtx(1) /'Print All Species Read from the Data File'/
          data uoprcx(1) /'EQ3NR, EQ6'/
c
          data uoprox(1,1) /"Don't print"/
          data ioprox(1,1) /0/
          data uoprox(2,1) /'Print'/
          data ioprox(2,1) /1/
          data uoprox(3,1) /' '/
          data ioprox(3,1) /0/
          data uoprox(4,1) /' '/
          data ioprox(4,1) /0/
          data uoprox(5,1) /' '/
          data ioprox(5,1) /0/
c
c       iopr(2):
c
          data uoprtx(2) /'Print All Reactions'/
          data uoprcx(2) /'EQ3NR, EQ6'/
c
          data uoprox(1,2) /"Don't print"/
          data ioprox(1,2) /0/
          data uoprox(2,2) /'Print the reactions'/
          data ioprox(2,2) /1/
          data uoprox(3,2) /'Print the reactions and log K values'/
          data ioprox(3,2) /2/
          data uoprox(4,2) /'Print the reactions, log K values, and asso
     $ciated data'/
          data ioprox(4,2) /3/
          data uoprox(5,2) /' '/
          data ioprox(5,2) /0/
c
c       iopr(3):
c
          data uoprtx(3) /'Print the Aqueous Species Hard Core Diameters
     $'/
          data uoprcx(3) /'EQ3NR, EQ6'/
c
          data uoprox(1,3) /"Don't print"/
          data ioprox(1,3) /0/
          data uoprox(2,3) /'Print'/
          data ioprox(2,3) /1/
          data uoprox(3,3) /' '/
          data ioprox(3,3) /0/
          data uoprox(4,3) /' '/
          data ioprox(4,3) /0/
          data uoprox(5,3) /' '/
          data ioprox(5,3) /0/
c
c       iopr(4):
c
          data uoprtx(4) /'Print a Table of Aqueous Species Concentratio
     $ns, Activities, etc.'/
          data uoprcx(4) /'EQ3NR, EQ6'/
c
          data uoprox(1,4) /'Omit species with molalities < 1.e-8'/
          data ioprox(1,4) /-3/
          data uoprox(2,4) /'Omit species with molalities < 1.e-12'/
          data ioprox(2,4) /-2/
          data uoprox(3,4) /'Omit species with molalities < 1.e-20'/
          data ioprox(3,4) /-1/
          data uoprox(4,4) /'Omit species with molalities < 1.e-100'/
          data ioprox(4,4) /-0/
          data uoprox(5,4) /'Include all species'/
          data ioprox(5,4) /1/
c
c       iopr(5):
c
          data uoprtx(5) /'Print a Table of Aqueous Species/H+ Activity
     $Ratios'/
          data uoprcx(5) /'EQ3NR, EQ6'/
c
          data uoprox(1,5) /"Don't print"/
          data ioprox(1,5) /0/
          data uoprox(2,5) /'Print cation/H+ activity ratios only'/
          data ioprox(2,5) /1/
          data uoprox(3,5) /'Print cation/H+ and anion/H+ activity ratio
     $s'/
          data ioprox(3,5) /2/
          data uoprox(4,5) /'Print ion/H+ activity ratios and neutral sp
     $ecies activities'/
          data ioprox(4,5) /3/
          data uoprox(5,5) /' '/
          data ioprox(5,5) /0/
c
c       iopr(6):
c
          data uoprtx(6) /'Print a Table of Aqueous Mass Balance Percent
     $ages'/
          data uoprcx(6) /'EQ3NR, EQ6'/
c
          data uoprox(1,6) /"Don't print"/
          data ioprox(1,6) /-1/
          data uoprox(2,6) /'Print those species comprising at least 99%
     $ of each mass balance'/
          data ioprox(2,6) /0/
          data uoprox(3,6) /'Print all contributing species'/
          data ioprox(3,6) /1/
          data uoprox(4,6) /' '/
          data ioprox(4,6) /0/
          data uoprox(5,6) /' '/
          data ioprox(5,6) /0/
c
c       iopr(7):
c
          data uoprtx(7) /'Print Tables of Saturation Indices and Affini
     $ties'/
          data uoprcx(7) /'EQ3NR, EQ6'/
c
          data uoprox(1,7) /"Don't print"/
          data ioprox(1,7) /-1/
          data uoprox(2,7) /'Print, omitting those phases undersaturated
     $ by more than 10 kcal'/
          data ioprox(2,7) /0/
          data uoprox(3,7) /'Print for all phases'/
          data ioprox(3,7) /1/
          data uoprox(4,7) /' '/
          data ioprox(4,7) /0/
          data uoprox(5,7) /' '/
          data ioprox(5,7) /0/
c
c       iopr(8):
c
          data uoprtx(8) /'Print a Table of Fugacities'/
          data uoprcx(8) /'EQ3NR, EQ6'/
c
          data uoprox(1,8) /"Don't print"/
          data ioprox(1,8) /-1/
          data uoprox(2,8) /'Print'/
          data ioprox(2,8) /0/
          data uoprox(3,8) /' '/
          data ioprox(3,8) /0/
          data uoprox(4,8) /' '/
          data ioprox(4,8) /0/
          data uoprox(5,8) /' '/
          data ioprox(5,8) /0/
c
c       iopr(9):
c
          data uoprtx(9) /'Print a Table of Mean Molal Activity Coeffici
     $ents'/
          data uoprcx(9) /'EQ3NR, EQ6'/
c
          data uoprox(1,9) /"Don't print"/
          data ioprox(1,9) /0/
          data uoprox(2,9) /'Print'/
          data ioprox(2,9) /1/
          data uoprox(3,9) /' '/
          data ioprox(3,9) /0/
          data uoprox(4,9) /' '/
          data ioprox(4,9) /0/
          data uoprox(5,9) /' '/
          data ioprox(5,9) /0/
c
c       iopr(10):
c
          data uoprtx(10) /'Print a Tabulation of the Pitzer Interaction
     $ Coefficients'/
          data uoprcx(10) /'EQ3NR, EQ6'/
c
          data uoprox(1,10) /"Don't print"/
          data ioprox(1,10) /0/
          data uoprox(2,10) /'Print a summary tabulation'/
          data ioprox(2,10) /1/
          data uoprox(3,10) /'Print a more detailed tabulation'/
          data ioprox(3,10) /2/
          data uoprox(4,10) /' '/
          data ioprox(4,10) /0/
          data uoprox(5,10) /' '/
          data ioprox(5,10) /0/
c
c       iopr(11):
c
          data uoprtx(11) /'Not Used'/
          data uoprcx(11) /'None'/
c
          data uoprox(1,11) /' '/
          data ioprox(1,11) /0/
          data uoprox(2,11) /' '/
          data ioprox(2,11) /0/
          data uoprox(3,11) /' '/
          data ioprox(3,11) /0/
          data uoprox(4,11) /' '/
          data ioprox(4,11) /0/
          data uoprox(5,11) /' '/
          data ioprox(5,11) /0/
c
c       iopr(12):
c
          data uoprtx(12) /'Not Used'/
          data uoprcx(12) /'None'/
c
          data uoprox(1,12) /' '/
          data ioprox(1,12) /0/
          data uoprox(2,12) /' '/
          data ioprox(2,12) /0/
          data uoprox(3,12) /' '/
          data ioprox(3,12) /0/
          data uoprox(4,12) /' '/
          data ioprox(4,12) /0/
          data uoprox(5,12) /' '/
          data ioprox(5,12) /0/
c
c       iopr(13):
c
          data uoprtx(13) /'Not Used'/
          data uoprcx(13) /'None'/
c
          data uoprox(1,13) /' '/
          data ioprox(1,13) /0/
          data uoprox(2,13) /' '/
          data ioprox(2,13) /0/
          data uoprox(3,13) /' '/
          data ioprox(3,13) /0/
          data uoprox(4,13) /' '/
          data ioprox(4,13) /0/
          data uoprox(5,13) /' '/
          data ioprox(5,13) /0/
c
c       iopr(14):
c
          data uoprtx(14) /'Not Used'/
          data uoprcx(14) /'None'/
c
          data uoprox(1,14) /' '/
          data ioprox(1,14) /0/
          data uoprox(2,14) /' '/
          data ioprox(2,14) /0/
          data uoprox(3,14) /' '/
          data ioprox(3,14) /0/
          data uoprox(4,14) /' '/
          data ioprox(4,14) /0/
          data uoprox(5,14) /' '/
          data ioprox(5,14) /0/
c
c       iopr(15):
c
          data uoprtx(15) /'Not Used'/
          data uoprcx(15) /'None'/
c
          data uoprox(1,15) /' '/
          data ioprox(1,15) /0/
          data uoprox(2,15) /' '/
          data ioprox(2,15) /0/
          data uoprox(3,15) /' '/
          data ioprox(3,15) /0/
          data uoprox(4,15) /' '/
          data ioprox(4,15) /0/
          data uoprox(5,15) /' '/
          data ioprox(5,15) /0/
c
c       iopr(16):
c
          data uoprtx(16) /'Not Used'/
          data uoprcx(16) /'None'/
c
          data uoprox(1,16) /' '/
          data ioprox(1,16) /0/
          data uoprox(2,16) /' '/
          data ioprox(2,16) /0/
          data uoprox(3,16) /' '/
          data ioprox(3,16) /0/
          data uoprox(4,16) /' '/
          data ioprox(4,16) /0/
          data uoprox(5,16) /' '/
          data ioprox(5,16) /0/
c
c       iopr(17):
c
          data uoprtx(17) /'PICKUP file format ("W" or "D")'/
          data uoprcx(17) /'EQ3NR, EQ6'/
c
          data uoprox(1,17) /'Use the format of the INPUT file'/
          data ioprox(1,17) /0/
          data uoprox(2,17) /'Use "W" format'/
          data ioprox(2,17) /1/
          data uoprox(3,17) /'Use "D" format'/
          data ioprox(3,17) /2/
          data uoprox(4,17) /' '/
          data ioprox(4,17) /0/
          data uoprox(5,17) /' '/
          data ioprox(5,17) /0/
c
c       iopr(18):
c
          data uoprtx(18) /'Not Used'/
          data uoprcx(18) /'None'/
c
          data uoprox(1,18) /' '/
          data ioprox(1,18) /0/
          data uoprox(2,18) /' '/
          data ioprox(2,18) /0/
          data uoprox(3,18) /' '/
          data ioprox(3,18) /0/
          data uoprox(4,18) /' '/
          data ioprox(4,18) /0/
          data uoprox(5,18) /' '/
          data ioprox(5,18) /0/
c
c       iopr(19):
c
          data uoprtx(19) /'Not Used'/
          data uoprcx(19) /'None'/
c
          data uoprox(1,19) /' '/
          data ioprox(1,19) /0/
          data uoprox(2,19) /' '/
          data ioprox(2,19) /0/
          data uoprox(3,19) /' '/
          data ioprox(3,19) /0/
          data uoprox(4,19) /' '/
          data ioprox(4,19) /0/
          data uoprox(5,19) /' '/
          data ioprox(5,19) /0/
c
c       iopr(20):
c
          data uoprtx(20) /'Not Used'/
          data uoprcx(20) /'None'/
c
          data uoprox(1,20) /' '/
          data ioprox(1,20) /0/
          data uoprox(2,20) /' '/
          data ioprox(2,20) /0/
          data uoprox(3,20) /' '/
          data ioprox(3,20) /0/
          data uoprox(4,20) /' '/
          data ioprox(4,20) /0/
          data uoprox(5,20) /' '/
          data ioprox(5,20) /0/
c
c-----------------------------------------------------------------------
c
c     Iodb print option strings.
c
c       iodb(1):
c
          data uodbtx(1) /'Print General Diagnostic Messages'/
          data uodbcx(1) /'EQ3NR, EQ6'/
c
          data uodbox(1,1) /"Don't print"/
          data iodbox(1,1) /0/
          data uodbox(2,1) /'Print Level 1 diagnostic messages'/
          data iodbox(2,1) /1/
          data uodbox(3,1) /'Print Level 1 and Level 2 diagnostic messag
     $es'/
          data iodbox(3,1) /2/
          data uodbox(4,1) /' '/
          data iodbox(4,1) /0/
          data uodbox(5,1) /' '/
          data iodbox(5,1) /0/
c
c       iodb(2):
c
          data uodbtx(2) /'Kinetics Related Diagnostic Messages'/
          data uodbcx(2) /'EQ6 only'/
c
          data uodbox(1,2) /"Don't print"/
          data iodbox(1,2) /0/
          data uodbox(2,2) /'Print Level 1 kinetics diagnostic messages'
     $/
          data iodbox(2,2) /1/
          data uodbox(3,2) /'Print Level 1 and Level 2 kinetics diagnost
     $ic messages'/
          data iodbox(3,2) /2/
          data uodbox(4,2) /' '/
          data iodbox(4,2) /0/
          data uodbox(5,2) /' '/
          data iodbox(5,2) /0/
c
c       iodb(3):
c
          data uodbtx(3) /'Print Pre-Newton-Raphson Optimization Informa
     $tion'/
          data uodbcx(3) /'EQ3NR, EQ6'/
c
          data uodbox(1,3) /"Don't print"/
          data iodbox(1,3) /0/
          data uodbox(2,3) /'Print summary information'/
          data iodbox(2,3) /1/
          data uodbox(3,3) /'Print detailed information (including the b
     $eta and del vectors)'/
          data iodbox(3,3) /2/
          data uodbox(4,3) /'Print more detailed information (including
     $matrix equations)'/
          data iodbox(4,3) /3/
          data uodbox(5,3) /'Print most detailed information (including
     $activity coefficients)'/
          data iodbox(5,3) /4/
c
c       iodb(4):
c
          data uodbtx(4) /'Print Newton-Raphson Iteration Information'/
          data uodbcx(4) /'EQ3NR, EQ6'/
c
          data uodbox(1,4) /"Don't print"/
          data iodbox(1,4) /0/
          data uodbox(2,4) /'Print summary information'/
          data iodbox(2,4) /1/
          data uodbox(3,4) /'Print detailed information (including the b
     $eta and del vectors)'/
          data iodbox(3,4) /2/
          data uodbox(4,4) /'Print more detailed information (including
     $the Jacobian)'/
          data iodbox(4,4) /3/
          data uodbox(5,4) /'Print most detailed information (including
     $activity coefficients)'/
          data iodbox(5,4) /4/
c
c       iodb(5):
c
          data uodbtx(5) /'Print Step-Size and Order Selection'/
          data uodbcx(5) /'EQ6 only'/
c
          data uodbox(1,5) /"Don't print"/
          data iodbox(1,5) /0/
          data uodbox(2,5) /'Print summary information'/
          data iodbox(2,5) /1/
          data uodbox(3,5) /'Print detailed information'/
          data iodbox(3,5) /2/
          data uodbox(4,5) /' '/
          data iodbox(4,5) /0/
          data uodbox(5,5) /' '/
          data iodbox(5,5) /0/
c
c       iodb(6):
c
          data uodbtx(6) /'Print Details of Hypothetical Affinity Calcul
     $ations'/
          data uodbcx(6) /'EQ3NR, EQ6'/
c
          data uodbox(1,6) /"Don't print"/
          data iodbox(1,6) /0/
          data uodbox(2,6) /'Print summary information'/
          data iodbox(2,6) /1/
          data uodbox(3,6) /'Print detailed information'/
          data iodbox(3,6) /2/
          data uodbox(4,6) /' '/
          data iodbox(4,6) /0/
          data uodbox(5,6) /' '/
          data iodbox(5,6) /0/
c
c       iodb(7):
c
          data uodbtx(7) /'Print General Search (e.g., for a phase bound
     $ary) Information'/
          data uodbcx(7) /'EQ6 only'/
c
          data uodbox(1,7) /"Don't print"/
          data iodbox(1,7) /0/
          data uodbox(2,7) /'Print summary information'/
          data iodbox(2,7) /1/
          data uodbox(3,7) /' '/
          data iodbox(3,7) /0/
          data uodbox(4,7) /' '/
          data iodbox(4,7) /0/
          data uodbox(5,7) /' '/
          data iodbox(5,7) /0/
c
c       iodb(8):
c
          data uodbtx(8) /'Print ODE Corrector Iteration Information'/
          data uodbcx(8) /'EQ6 only'/
c
          data uodbox(1,8) /"Don't print"/
          data iodbox(1,8) /0/
          data uodbox(2,8) /'Print summary information'/
          data iodbox(2,8) /1/
          data uodbox(3,8) /'Print detailed information (including the b
     $etar and delvcr vectors)'/
          data iodbox(3,8) /2/
          data uodbox(4,8) /' '/
          data iodbox(4,8) /0/
          data uodbox(5,8) /' '/
          data iodbox(5,8) /0/
c
c       iodb(9):
c
          data uodbtx(9) /'Not Used'/
          data uodbcx(9) /'None'/
c
          data uodbox(1,9) /' '/
          data iodbox(1,9) /0/
          data uodbox(2,9) /' '/
          data iodbox(2,9) /0/
          data uodbox(3,9) /' '/
          data iodbox(3,9) /0/
          data uodbox(4,9) /' '/
          data iodbox(4,9) /0/
          data uodbox(5,9) /' '/
          data iodbox(5,9) /0/
c
c       iodb(10):
c
          data uodbtx(10) /'Not Used'/
          data uodbcx(10) /'None'/
c
          data uodbox(1,10) /' '/
          data iodbox(1,10) /0/
          data uodbox(2,10) /' '/
          data iodbox(2,10) /0/
          data uodbox(3,10) /' '/
          data iodbox(3,10) /0/
          data uodbox(4,10) /' '/
          data iodbox(4,10) /0/
          data uodbox(5,10) /' '/
          data iodbox(5,10) /0/
c
c       iodb(11):
c
          data uodbtx(11) /'Not Used'/
          data uodbcx(11) /'None'/
c
          data uodbox(1,11) /' '/
          data iodbox(1,11) /0/
          data uodbox(2,11) /' '/
          data iodbox(2,11) /0/
          data uodbox(3,11) /' '/
          data iodbox(3,11) /0/
          data uodbox(4,11) /' '/
          data iodbox(4,11) /0/
          data uodbox(5,11) /' '/
          data iodbox(5,11) /0/
c
c       iodb(12):
c
          data uodbtx(12) /'Not Used'/
          data uodbcx(12) /'None'/
c
          data uodbox(1,12) /' '/
          data iodbox(1,12) /0/
          data uodbox(2,12) /' '/
          data iodbox(2,12) /0/
          data uodbox(3,12) /' '/
          data iodbox(3,12) /0/
          data uodbox(4,12) /' '/
          data iodbox(4,12) /0/
          data uodbox(5,12) /' '/
          data iodbox(5,12) /0/
c
c       iodb(13):
c
          data uodbtx(13) /'Not Used'/
          data uodbcx(13) /'None'/
c
          data uodbox(1,13) /' '/
          data iodbox(1,13) /0/
          data uodbox(2,13) /' '/
          data iodbox(2,13) /0/
          data uodbox(3,13) /' '/
          data iodbox(3,13) /0/
          data uodbox(4,13) /' '/
          data iodbox(4,13) /0/
          data uodbox(5,13) /' '/
          data iodbox(5,13) /0/
c
c       iodb(14):
c
          data uodbtx(14) /'Not Used'/
          data uodbcx(14) /'None'/
c
          data uodbox(1,14) /' '/
          data iodbox(1,14) /0/
          data uodbox(2,14) /' '/
          data iodbox(2,14) /0/
          data uodbox(3,14) /' '/
          data iodbox(3,14) /0/
          data uodbox(4,14) /' '/
          data iodbox(4,14) /0/
          data uodbox(5,14) /' '/
          data iodbox(5,14) /0/
c
c       iodb(15):
c
          data uodbtx(15) /'Not Used'/
          data uodbcx(15) /'None'/
c
          data uodbox(1,15) /' '/
          data iodbox(1,15) /0/
          data uodbox(2,15) /' '/
          data iodbox(2,15) /0/
          data uodbox(3,15) /' '/
          data iodbox(3,15) /0/
          data uodbox(4,15) /' '/
          data iodbox(4,15) /0/
          data uodbox(5,15) /' '/
          data iodbox(5,15) /0/
c
c       iodb(16):
c
          data uodbtx(16) /'Not Used'/
          data uodbcx(16) /'None'/
c
          data uodbox(1,16) /' '/
          data iodbox(1,16) /0/
          data uodbox(2,16) /' '/
          data iodbox(2,16) /0/
          data uodbox(3,16) /' '/
          data iodbox(3,16) /0/
          data uodbox(4,16) /' '/
          data iodbox(4,16) /0/
          data uodbox(5,16) /' '/
          data iodbox(5,16) /0/
c
c       iodb(17):
c
          data uodbtx(17) /'Not Used'/
          data uodbcx(17) /'None'/
c
          data uodbox(1,17) /' '/
          data iodbox(1,17) /0/
          data uodbox(2,17) /' '/
          data iodbox(2,17) /0/
          data uodbox(3,17) /' '/
          data iodbox(3,17) /0/
          data uodbox(4,17) /' '/
          data iodbox(4,17) /0/
          data uodbox(5,17) /' '/
          data iodbox(5,17) /0/
c
c       iodb(18):
c
          data uodbtx(18) /'Not Used'/
          data uodbcx(18) /'None'/
c
          data uodbox(1,18) /' '/
          data iodbox(1,18) /0/
          data uodbox(2,18) /' '/
          data iodbox(2,18) /0/
          data uodbox(3,18) /' '/
          data iodbox(3,18) /0/
          data uodbox(4,18) /' '/
          data iodbox(4,18) /0/
          data uodbox(5,18) /' '/
          data iodbox(5,18) /0/
c
c       iodb(19):
c
          data uodbtx(19) /'Not Used'/
          data uodbcx(19) /'None'/
c
          data uodbox(1,19) /' '/
          data iodbox(1,19) /0/
          data uodbox(2,19) /' '/
          data iodbox(2,19) /0/
          data uodbox(3,19) /' '/
          data iodbox(3,19) /0/
          data uodbox(4,19) /' '/
          data iodbox(4,19) /0/
          data uodbox(5,19) /' '/
          data iodbox(5,19) /0/
c
c       iodb(20):
c
          data uodbtx(20) /'Not Used'/
          data uodbcx(20) /'None'/
c
          data uodbox(1,20) /' '/
          data iodbox(1,20) /0/
          data uodbox(2,20) /' '/
          data iodbox(2,20) /0/
          data uodbox(3,20) /' '/
          data iodbox(3,20) /0/
          data uodbox(4,20) /' '/
          data iodbox(4,20) /0/
          data uodbox(5,20) /' '/
          data iodbox(5,20) /0/
c
c     End of INCLUDE file eqlo8d.h
c-----------------------------------------------------------------------
