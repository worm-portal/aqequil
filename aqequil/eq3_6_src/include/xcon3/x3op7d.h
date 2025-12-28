c x3op7d.h
c
c     Strings, arrays, and such for the menu-style ("D") input
c     format for version 7.
c
c***********************************************************************
c      VALIDATION VALUES FOR THE VARIABLE cspb                         *
c***********************************************************************
c
      data ujflg7(-1) / 'Suppress          ' /
      data ujflg7(0)  / 'Molality          ' /
      data ujflg7(1)  / 'Molarity          ' /
      data ujflg7(2)  / 'mg/L              ' /
      data ujflg7(3)  / 'mg/kg             ' /
      data ujflg7(4)  / 'Free molal        ' /
      data ujflg7(5)  / 'Free molar        ' /
      data ujflg7(6)  / 'Free mg/L         ' /
      data ujflg7(7)  / 'Free mg/kg        ' /
      data ujflg7(8)  / 'Free cm3/cm3      ' /
      data ujflg7(9)  / 'ERROR             ' /
      data ujflg7(10) / 'ERROR             ' /
      data ujflg7(11) / 'ERROR             ' /
      data ujflg7(12) / 'ERROR             ' /
      data ujflg7(13) / 'ERROR             ' /
      data ujflg7(14) / 'ERROR             ' /
      data ujflg7(15) / 'ERROR             ' /
      data ujflg7(16) / 'Log activity      ' /
      data ujflg7(17) / 'Log activity combo' /
      data ujflg7(18) / 'Log mean activity ' /
      data ujflg7(19) / 'Mineral           ' /
      data ujflg7(20) / 'Solid solution    ' /
      data ujflg7(21) / 'Log fugacity      ' /
      data ujflg7(22) / 'ERROR             ' /
      data ujflg7(23) / 'ERROR             ' /
      data ujflg7(24) / 'ERROR             ' /
      data ujflg7(25) / 'ERROR             ' /
      data ujflg7(26) / 'ERROR             ' /
      data ujflg7(27) / 'Dependent         ' /
      data ujflg7(28) / 'ERROR             ' /
      data ujflg7(29) / 'ERROR             ' /
      data ujflg7(30) / 'Eliminated        ' /
      data ujflg7(31) / 'pH                ' /
      data ujflg7(32) / 'pHCl              ' /
c
c***********************************************************************
c                         O  P  T  I  O  N  S                          *
c***********************************************************************
c
      data uopt3(1)
     $   / 'SOLID SOLUTIONS                                 ' /
      data uvar3(1)     / 'iopt' /
      data index3(1)   / 4 /
      data udesc3(1)    / 'ignore solid solutions' /
      data iopti3(1)   / 1 /
      data ivalu3(1)   / 0 /
      data udesc3(2)    / 'process hypothetical solid solutions' /
      data iopti3(2)   / 1 /
      data ivalu3(2)   / 1 /
      data udesc3(3)
     $        / 'process input and hypothetical solid solutions'/
      data iopti3(3)  / 1 /
      data ivalu3(3)  / 2 /
c
      data uopt3(2)
     $   / 'LOADING OF SPECIES INTO MEMORY                  ' /
      data uvar3(2)     / 'iopr'  /
      data index3(2)   / 1 /
      data udesc3(4)    / 'does nothing' /
      data iopti3(4)   / 2 /
      data ivalu3(4)   / 0 /
      data udesc3(5)    / 'lists species loaded into memory' /
      data iopti3(5)   / 2 /
      data ivalu3(5)   / 1 /
c
      data uopt3(3)   / 'ECHO DATABASE INFORMATION' /
      data uvar3(3)     / 'iopr' /
      data index3(3)   / 2 /
      data udesc3(6)    / 'does nothing' /
      data iopti3(6)   / 3 /
      data ivalu3(6)   / 0 /
      data udesc3(7)    / 'lists all reactions' /
      data iopti3(7)   / 3 /
      data ivalu3(7)   / 1 /
      data udesc3(8)    / 'lists reactions and log K values' /
      data iopti3(8)   / 3 /
      data ivalu3(8)   / 2 /
      data udesc3(9)
     $           / 'lists reactions, log K values and polynomial coef.'/
      data iopti3(9)   / 3 /
      data ivalu3(9)   / 3 /
c
      data uopt3(4)   / 'LIST OF AQUEOUS SPECIES (ordering)' /
      data uvar3(4)     / 'iopr' /
      data index3(4)   / 3 /
      data udesc3(10)    / 'in order of decreasing concentration' /
      data iopti3(10)   / 4 /
      data ivalu3(10)   / 0 /
      data udesc3(11)   / 'in same order as input file' /
      data iopti3(11)  / 4 /
      data ivalu3(11)  / 1 /
c
      data uopt3(5)
     $            / 'LIST OF AQUEOUS SPECIES (concentration limit)' /
      data uvar3(5)      / 'iopr' /
      data index3(5)    / 4 /
      data udesc3(12)    / 'all species' /
      data iopti3(12)   / 5 /
      data ivalu3(12)   / 0 /
      data udesc3(13)    / 'only species > 10**-20 molal' /
      data iopti3(13)   / 5 /
      data ivalu3(13)   / -1 /
      data udesc3(14)    / 'only species > 10**-12 molal' /
      data iopti3(14)   / 5 /
      data ivalu3(14)   / -2 /
      data udesc3(15)    / 'not printed' /
      data iopti3(15)   / 5 /
      data ivalu3(15)   / -3 /
c
      data uopt3(6)    / 'LIST OF AQUEOUS SPECIES (by element)    ' /
      data uvar3(6)      / 'iopr' /
      data index3(6)    / 5 /
      data udesc3(16)    / 'print major species' /
      data iopti3(16)   / 6 /
      data ivalu3(16)   / 0 /
      data udesc3(17)    / 'print all species' /
      data iopti3(17)   / 6 /
      data ivalu3(17)   / 1 /
      data udesc3(18)    / 'don''t print' /
      data iopti3(18)   / 6 /
      data ivalu3(18)   / -1 /
c
      data uopt3(7)    / 'MINERAL SATURATION STATES               ' /
      data uvar3(7)      / 'iopr' /
      data index3(7)    / 7 /
      data udesc3(19)    / 'print if affinity > -10 kcals' /
      data iopti3(19)   / 7 /
      data ivalu3(19)   / 0 /
      data udesc3(20)    / 'print all' /
      data iopti3(20)   / 7 /
      data ivalu3(20)   / 1 /
      data udesc3(21)    / 'don''t print' /
      data iopti3(21)   / 7 /
      data ivalu3(21)   / -1 /
c
      data uopt3(8)    / 'pH SCALE CONVENTION                     ' /
      data uvar3(8)      / 'iopg' /
      data index3(8)    / 2 /
      data udesc3(22)    / 'modified NBS' /
      data iopti3(22)   / 8 /
      data ivalu3(22)   / 0 /
      data udesc3(23)    / 'internal' /
      data iopti3(23)   / 8 /
      data ivalu3(23)   / -1 /
      data udesc3(24)    / 'rational' /
      data iopti3(24)   / 8 /
      data ivalu3(24)   / 1 /
c
      data uopt3(9)    / 'ACTIVITY COEFFICIENT OPTIONS' /
      data uvar3(9)      / 'iopg' /
      data index3(9)    / 1 /
      data udesc3(25)    / 'use B-dot equation' /
      data iopti3(25)   / 9 /
      data ivalu3(25)   / 0 /
      data udesc3(26)    / 'Davies'' equation' /
      data iopti3(26)   / 9 /
      data ivalu3(26)   / -1 /
      data udesc3(27)    / 'Pitzer''s equations' /
      data iopti3(27)   / 9 /
      data ivalu3(27)   / 1 /
c
      data uopt3(10)   / 'AUTO BASIS SWITCHING' /
      data uvar3(10)     / 'iopt' /
      data index3(10)   / 2 /
      data udesc3(28)    / 'off' /
      data iopti3(28)   / 10 /
      data ivalu3(28)   / 0 /
      data udesc3(29)    / 'on' /
      data iopti3(29)   / 10 /
      data ivalu3(29)   / 1 /
c
      data uopt3(11) / 'PITZER DATABASE INFORMATION' /
      data uvar3(11)     / 'iopr' /
      data index3(11)   / 9 /
      data udesc3(30)    / 'print only warnings' /
      data iopti3(30)   / 11 /
      data ivalu3(30)   / 0 /
      data udesc3(31)
     $  / 'print species in model and number of Pitzer coefficients' /
      data iopti3(31)   / 11 /
      data ivalu3(31)   / 1 /
      data udesc3(32)
     $  / 'print species in model and names of Pitzer coefficients' /
      data iopti3(32)   / 11 /
      data ivalu3(32)   / 2 /
c
      data uopt3(12)   / 'PICKUP FILE' /
      data uvar3(12)    / 'iopt' /
      data index3(12)   / 3 /
      data udesc3(33)   / 'write pickup file' /
      data iopti3(33)   / 12 /
      data ivalu3(33)   / 0 /
      data udesc3(34)   / 'don''t write pickup file' /
      data iopti3(34)   / 12 /
      data ivalu3(34)   / -1 /
c
      data uopt3(13) / 'LIST MEAN IONIC PROPERTIES' /
      data uvar3(13)     / 'iopr' /
      data index3(13)   / 6 /
      data udesc3(35)    / 'don''t print' /
      data iopti3(35)   / 13 /
      data ivalu3(35)   / 0  /
      data udesc3(36)    / 'print' /
      data iopti3(36)   / 13 /
      data ivalu3(36)   / 1 /
c
      data uopt3(14)
     $  / 'LIST AQUEOUS SPECIES, ION SIZES, AND HYDRATION NUMBERS'/
      data uvar3(14)    / 'iopr' /
      data index3(14)   /  8  /
      data udesc3(37)   / 'print' /
      data iopti3(37)   /  14     /
      data ivalu3(37)   /  0      /
      data udesc3(38)   / 'don''t print'  /
      data iopti3(38)   / 14     /
      data ivalu3(38)   / -1     /
c
      data uopt3(15)  / 'CONVERGENCE CRITERIA' /
      data uvar3(15)    / 'iopt' /
      data index3(15)   / 6      /
      data udesc3(39)
     $   / 'test both residual functions and correction terms' /
      data iopti3(39)   / 15     /
      data ivalu3(39)   / 0      /
      data udesc3(40)   / 'test only residual functions' /
      data iopti3(40)   / 15     /
      data ivalu3(40)   / -1     /
c
c*********************************************************************
c     D E B U G      O P T I O N S
c*********************************************************************
c
      data udebug(1)  / 'generic debugging information' /
      data idbugi(1)  / 1 /
c
      data udebug(2)
     $   / 'print details of pre-Newton-Raphson iteration' /
      data idbugi(2)  / 2 /
c
      data udebug(3)  /'print details of Newton-Raphson iteration'/
      data idbugi(3)  / 4 /
c
      data udebug(4)  /'print details of stoichiometric factors'/
      data idbugi(4)  / 5 /
c
      data udebug(5)
     $   /'print details of stoichiometric factors calculation' /
      data idbugi(5)  / 6 /
c
      data udebug(6)  /'write reactions on RLIST'   /
      data idbugi(6)  / 7 /
c
      data udebug(7)
     $   /'list stoichiometric concentrations of master species' /
      data idbugi(7)  / 10 /
c
      data udebug(8) /'request iteration variables to be killed' /
      data idbugi(8) / 3 /
c
c**********************************************************************
c          D  E  V  E  L  O  P  M  E  N  T      O  P  T  I  O  N  S   *
c**********************************************************************
c
c     None
c
c**********************************************************************
c                      T  O  L  E  R  A  N  C  E S                    *
c**********************************************************************
c
      data utol3(1)   / '      residual functions' /
      data utol3(2)   / '        correction terms' /
      data utol3(3)   / '        saturation state' /
      data utol3(4)   / 'number of N-R iterations' /
c
c     End of include file x3op7d.h
c-----------------------------------------------------------------------
