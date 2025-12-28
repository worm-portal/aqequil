c x6op7d.h
c
c     Strings, arrays, and such for the menu-style ("D") input
c     format for version 7.
c
c***********************************************************************
c                         O  P  T  I  O  N  S                          *
c***********************************************************************
c
      data uopt6(1)
     $   / 'SOLID SOLUTIONS                                 ' /
      data uvar6(1)     / 'iopt' /
      data index6(1)   / 4 /
      data udesc6(1)    / 'ignore solid solutions' /
      data iopti6(1)   / 1 /
      data ivalu6(1)   / 0 /
      data udesc6(2)    / 'process solid solutions' /
      data iopti6(2)   / 1 /
      data ivalu6(2)   / 1 /
c
      data uopt6(2)
     $   / 'LOADING OF SPECIES INTO MEMORY                  ' /
      data uvar6(2)     / 'iopr'  /
      data index6(2)   / 1 /
      data udesc6(3)    / "don't print" /
      data iopti6(3)   / 2 /
      data ivalu6(3)   / 0 /
      data udesc6(4)    / 'lists species loaded into memory' /
      data iopti6(4)   / 2 /
      data ivalu6(4)   / 1 /
c
      data uopt6(3)
     $  / 'LIST DERIVATIVES OF BASIS ELEMENTS AT EACH PRINT POINT' /
      data uvar6(3)     / 'iopr' /
      data index6(3)   / 2 /
      data udesc6(5)    / "don't print" /
      data iopti6(5)   / 3 /
      data ivalu6(5)   / 0 /
      data udesc6(6)    / 'print' /
      data iopti6(6)   / 3 /
      data ivalu6(6)   / 1 /
c
      data uopt6(4)
     $ / 'LIST ALL SPECIES LOADED INTO MEMORY AND THEIR LOG K VALUES' /
      data uvar6(4)     / 'iopr' /
      data index6(4)   / 3 /
      data udesc6(7)    / "don't print" /
      data iopti6(7)   / 4 /
      data ivalu6(7)   / 0 /
      data udesc6(8)   / 'print' /
      data iopti6(8)  / 4 /
      data ivalu6(8)  / 1 /
c
      data uopt6(5)
     $   / 'LIST DISTRIBUTION OF AQUEOUS SPECIES AT EACH PRINT POINT' /
      data uvar6(5)      / 'iopr' /
      data index6(5)     / 4 /
      data udesc6(9)     / 'only species > 10**-12 molal' /
      data iopti6(9)     / 5 /
      data ivalu6(9)     / 0 /
      data udesc6(10)    / 'all species' /
      data iopti6(10)    / 5 /
      data ivalu6(10)    / 1 /
      data udesc6(11)    / "don't print" /
      data iopti6(11)    / 5 /
      data ivalu6(11)    / -1 /
c
      data uopt6(6)
     $     / 'LIST CATION/H+ ACTIVITY RATIOS AT EACH PRINT POINT' /
      data uvar6(6)      / 'iopr' /
      data index6(6)    / 5 /
      data udesc6(12)    / "don't print" /
      data iopti6(12)   / 6 /
      data ivalu6(12)   / 0 /
      data udesc6(13)    / 'print' /
      data iopti6(13)   / 6 /
      data ivalu6(13)   / 1 /
c
      data uopt6(7)
     $ / 'LIST BULK ELEMENT AND OXIDE COMPOSITION AT EACH PRINT POINT'/
      data uvar6(7)      / 'iopr' /
      data index6(7)    / 6 /
      data udesc6(14)    / "don't print" /
      data iopti6(14)   / 7 /
      data ivalu6(14)   / 0 /
      data udesc6(15)    / 'print' /
      data iopti6(15)   / 7 /
      data ivalu6(15)   / 1 /
c
      data uopt6(8)    / 'MINERAL SATURATION STATES               ' /
      data uvar6(8)      / 'iopr' /
      data index6(8)    / 7 /
      data udesc6(16)    / 'print if affinity > -10 kcals' /
      data iopti6(16)   / 8 /
      data ivalu6(16)   / 0 /
      data udesc6(17)    / 'print all' /
      data iopti6(17)   / 8 /
      data ivalu6(17)   / 1 /
      data udesc6(18)    / "don't print" /
      data iopti6(18)   / 8 /
      data ivalu6(18)   / -1 /
c
      data uopt6(9)
     $   / 'LIST GAS SPECIES SUMMARY AT EACH PRINT POINT' /
      data uvar6(9)      / 'iopr' /
      data index6(9)    / 8 /
      data udesc6(19)    / "don't print" /
      data iopti6(19)   / 9 /
      data ivalu6(19)   / 0 /
      data udesc6(20)    / 'print' /
      data iopti6(20)   / 9 /
      data ivalu6(20)   / 1 /
c
      data uopt6(10)
     $   / 'PRINT AQUEOUS MASS AND CONCENTRATION TOTALS' /
      data uvar6(10)    / 'iopr' /
      data index6(10)   / 11 /
      data udesc6(21)   / "don't print" /
      data iopti6(21)   / 10 /
      data ivalu6(21)   / 0 /
      data udesc6(22)   / 'print' /
      data iopti6(22)   / 10 /
      data ivalu6(22)   / 1 /
c
      data uopt6(11)   / 'TAB FILES' /
      data uvar6(11)     / 'iopt' /
      data index6(11)   / 13 /
      data udesc6(23)    / 'write' /
      data iopti6(23)   / 11 /
      data ivalu6(23)   / 0 /
      data udesc6(24)    / 'append to previous tabx file' /
      data iopti6(24)   / 11 /
      data ivalu6(24)   / 1 /
      data udesc6(25)    / "don't write" /
      data iopti6(25)   / 11 /
      data ivalu6(25)   / -1 /
c
      data uopt6(12)   / 'WRITE PICKUP FILE' /
      data uvar6(12)     / 'iopt' /
      data index6(12)   / 3 /
      data udesc6(26)    / 'write pickup file at end of run' /
      data iopti6(26)   / 12 /
      data ivalu6(26)   / 0 /
      data udesc6(27)    / "don't write pickup file" /
      data iopti6(27)   / 12 /
      data ivalu6(27)   / -1 /
      data udesc6(28)
     $   / 'write pickup file for each print point' /
      data iopti6(28)   / 12 /
      data ivalu6(28)   / 1  /
c
      data uopt6(13)   / 'PHYSICALLY REMOVED SUBSYSTEM' /
      data uvar6(13)     / 'iopt' /
      data index6(13)   / 5 /
      data udesc6(29)    / 'does nothing' /
      data iopti6(29)   / 13 /
      data ivalu6(29)   / 0 /
      data udesc6(30)
     $  / 'transfer minerals but leave trivial mass in the system' /
      data iopti6(30)   / 13 /
      data ivalu6(30)   / 1 /
      data udesc6(31)
     $  / 'transfer minerals' /
      data iopti6(31)   / 13 /
      data ivalu6(31)   / 2 /
c
      data uopt6(14)   / 'CLEAR INITIAL PHYSICALLY REMOVED SUBSYSTEM' /
      data uvar6(14)     / 'iopt' /
      data index6(14)   / 6 /
      data udesc6(32)    / 'does nothing' /
      data iopti6(32)   / 14 /
      data ivalu6(32)   / 0 /
      data udesc6(33)
     $       / 'clear p.r.s. before first reaction progress advance' /
      data iopti6(33)   / 14 /
      data ivalu6(33)   / 1 /
c
      data uopt6(15) / 'PHASE BOUNDARY SEARCH' /
      data uvar6(15)     / 'iopt' /
      data index6(15)   / 2 /
      data udesc6(34)
     $  / 'step size constrained by predicted phase boundaries' /
      data iopti6(34)   / 15 /
      data ivalu6(34)   / 0  /
      data udesc6(35)
     $ /"phase boundaries estimated from Taylor's series and printed"/
      data iopti6(35)   / 15 /
      data ivalu6(35)   / 1 /
      data udesc6(36)    / 'locations of phase boundaries ignored'/
      data iopti6(36)   / 15 /
      data ivalu6(36)   / 2 /
c
      data uopt6(16)   / 'AUTO BASIS SWITCHING' /
      data uvar6(16)     / 'iopt' /
      data index6(16)   / 7 /
      data udesc6(37)    / 'off' /
      data iopti6(37)   / 16 /
      data ivalu6(37)   / 0  /
      data udesc6(38)    / 'on' /
      data iopti6(38)   / 16 /
      data ivalu6(38)   / 1 /
c
      data uopt6(17)   / 'SUPPRESS REDOX REACTIONS'/
      data uvar6(17)     / 'iopt' /
      data index6(17)   /  11 /
      data udesc6(39)    / 'does nothing' /
      data iopti6(39)   /  17     /
      data ivalu6(39)   /  0      /
      data udesc6(40)    / 'suppress all redox reactions'  /
      data iopti6(40)   /  17     /
      data ivalu6(40)   /  1     /
c
      data uopt6(18)   / "LINEAR OR LOGARITHMIC TAYLOR'S SERIES" /
      data uvar6(18)     / 'iopt' /
      data index6(18)   / 8      /
      data udesc6(41)
     $   / 'linear for kcol = 1,kdim, logarithmic for kcol = 1,kbt' /
      data iopti6(41)   / 18    /
      data ivalu6(41)   / 0      /
      data udesc6(42)    / 'logarithmic for kcol = 1,kbt' /
      data iopti6(42)   / 18     /
      data ivalu6(42)   / 1      /
      data udesc6(43)    / 'linear for kcol = 1,kdim' /
      data iopti6(43)   / 18     /
      data ivalu6(43)   / -1      /
c
      data uopt6(19)   / 'AZERO AND HYDRATION NUMBERS' /
      data uvar6(19)     / 'iopt' /
      data index6(19)   / 14      /
      data udesc6(44)    / 'no change' /
      data iopti6(44)   / 19    /
      data ivalu6(44)   / 0      /
      data udesc6(45)    / 'read in new azero and hydration numbers' /
      data iopti6(45)   / 19     /
      data ivalu6(45)   / 1      /
c
      data uopt6(20)
     $/ 'PRINT MEAN MOLAL ACTIVITY COEFFICIENTS FOR DISSOLVED SPECIES'/
      data uvar6(20)     / 'iopr' /
      data index6(20)   / 9      /
      data udesc6(46)    / 'does nothing' /
      data iopti6(46)   / 20     /
      data ivalu6(46)   / 0      /
      data udesc6(47)    / 'print'/
      data iopti6(47)   / 20     /
      data ivalu6(47)   / 1      /
c
      data uopt6(21)   / 'PITZER DATABASE INFORMATION'/
      data uvar6(21)     / 'iopr' /
      data index6(21)   / 10     /
      data udesc6(48)    / 'print only warnings' /
      data iopti6(48)   / 21     /
      data ivalu6(48)   / 0      /
      data udesc6(49)
     $   / 'print species in model and number of Pitzer coefficients' /
      data iopti6(49)   / 21     /
      data ivalu6(49)   / 1      /
      data udesc6(50)
     $   / 'print species in model and names of Pitzer coefficients' /
      data iopti6(50)   / 21     /
      data ivalu6(50)   / 2      /
c
      data uopt6(22)   / 'PRINT DIAGNOSTIC MESSAGES'/
      data uvar6(22)     / 'iodb' /
      data index6(22)   / 1      /
      data udesc6(51)    / "don't print" /
      data iopti6(51)   / 22    /
      data ivalu6(51)   / 0      /
      data udesc6(52)    / 'print level 1 messages' /
      data iopti6(52)   / 22     /
      data ivalu6(52)   / 1      /
      data udesc6(53)    / 'print level 2 messages' /
      data iopti6(53)   / 22     /
      data ivalu6(53)   / 2      /
c
      data uopt6(23)   / 'PRINT PRE-NEWTON-RAPHSON OPTIMIZATION'/
      data uvar6(23)     / 'iodb' /
      data index6(23)   / 2     /
      data udesc6(54)    / "don't print" /
      data iopti6(54)   / 23    /
      data ivalu6(54)   / 0      /
      data udesc6(55)    / 'print summary information' /
      data iopti6(55)   / 23     /
      data ivalu6(55)   / 1      /
      data udesc6(56)    / 'print detailed information' /
      data iopti6(56)   / 23     /
      data ivalu6(56)   / 2      /
c
      data uopt6(24)   / 'PRINT STEP SIZE AND ORDER'/
      data uvar6(24)     / 'iodb' /
      data index6(24)   / 3     /
      data udesc6(57)    / "don't print" /
      data iopti6(57)   / 24    /
      data ivalu6(57)   / 0      /
      data udesc6(58)    / 'print scale factor' /
      data iopti6(58)   / 24     /
      data ivalu6(58)   / 1      /
      data udesc6(59)    / 'print orders and step size scaling factors'/
      data iopti6(59)   / 24     /
      data ivalu6(59)   / 2      /
c
      data uopt6(25)   / 'CONTROL STEP SIZE AND ORDER PRINT'/
      data uvar6(25)     / 'iodb' /
      data index6(25)   / 8     /
      data udesc6(60)    / 'does nothing' /
      data iopti6(60)   / 25    /
      data ivalu6(60)   / 0      /
      data udesc6(61)
     $          / 'print step size and order when delzi .le. dlzmx1' /
      data iopti6(61)   / 25     /
      data ivalu6(61)   / 1      /
c
      data uopt6(26)   / 'NEWTON ITERATIONS'/
      data uvar6(26)     / 'iodb' /
      data index6(26)   / 4     /
      data udesc6(62)    / "don't print" /
      data iopti6(62)   / 26    /
      data ivalu6(62)   / 0      /
      data udesc6(63)    / 'print summary of newton iterations' /
      data iopti6(63)   / 26     /
      data ivalu6(63)   / 1      /
      data udesc6(64)
     $   / 'print summary, residual functions and correction terms' /
      data iopti6(64)   / 26    /
      data ivalu6(64)   / 2     /
      data udesc6(65)
     $/'print summary, residual functions, correction terms and matrix'/
      data iopti6(65)   / 26     /
      data ivalu6(65)   / 4      /
c
      data uopt6(27)   / 'PRINT SEARCH ITERATIONS'/
      data uvar6(27)     / 'iodb' /
      data index6(27)   / 5     /
      data udesc6(66)    / "don't print" /
      data iopti6(66)   / 27    /
      data ivalu6(66)   / 0      /
      data udesc6(67)    / 'print' /
      data iopti6(67)   / 27     /
      data ivalu6(67)   / 1      /
c
      data uopt6(28)   / 'PRINT HPSAT ITERATIONS'/
      data uvar6(28)     / 'iodb' /
      data index6(28)   / 6     /
      data udesc6(68)    / "don't print" /
      data iopti6(68)   / 28    /
      data ivalu6(68)   / 0      /
      data udesc6(69)    / 'print' /
      data iopti6(69)   / 28     /
      data ivalu6(69)   / 1      /
c
      data uopt6(29)   / 'PRINT FINITE DIFFERENCE AND DERIVATIVE DATA'/
      data uvar6(29)     / 'iodb' /
      data index6(29)   / 7     /
      data udesc6(70)    / "don't print" /
      data iopti6(70)   / 29    /
      data ivalu6(70)   / 0      /
      data udesc6(71)   / 'print computations from RDERIV, and RTAYLR' /
      data iopti6(71)   / 29     /
      data ivalu6(71)   / 1      /
      data udesc6(72)
     $   / 'print computations from RDERIV, RTAYLR, DERIV and TAYLOR' /
      data iopti6(72)   / 29     /
      data ivalu6(72)   / 2      /
c
      data uopt6(30)   / 'PRINT KINETICS DIAGNOSTIC MESSAGES'/
      data uvar6(30)     / 'iodb' /
      data index6(30)   / 9     /
      data udesc6(73)    / "don't print" /
      data iopti6(73)   / 30    /
      data ivalu6(73)   / 0      /
      data udesc6(74)    / 'print level 1 diagnostics' /
      data iopti6(74)   / 30     /
      data ivalu6(74)   / 1      /
      data udesc6(75)    / 'print level 1 and level 2 diagnostics' /
      data iopti6(75)   / 30     /
      data ivalu6(75)   / 2      /
c
      data uopt6(31)   / 'PRINT AKMATR'/
      data uvar6(31)     / 'iodb' /
      data index6(31)   / 16     /
      data udesc6(76)    / "don't print" /
      data iopti6(76)   / 31     /
      data ivalu6(76)   / 0      /
      data udesc6(77)    / 'print level 1 diagnostics' /
      data iopti6(77)   / 31     /
      data ivalu6(77)   / 1      /
c
      data uopt6(32)
     $      / 'KILL ITERATION VARIABLES'/
      data uvar6(32)     / 'iodb' /
      data index6(32)    / 12     /
      data udesc6(78)    / 'does nothing' /
      data iopti6(78)    / 32     /
      data ivalu6(78)    / 0      /
      data udesc6(79)    / 'allow selection of variables to remove' /
      data iopti6(79)    / 32     /
      data ivalu6(79)    / 1      /
c
      data uopt6(33)    / 'pH SCALE CONVENTION                     ' /
      data uvar6(33)      / 'iopg' /
      data index6(33)    / 2 /
      data udesc6(80)    / 'modified NBS' /
      data iopti6(80)   / 33 /
      data ivalu6(80)   / 0 /
      data udesc6(81)    / 'internal' /
      data iopti6(81)   / 33 /
      data ivalu6(81)   / -1 /
      data udesc6(82)    / 'rational' /
      data iopti6(82)   / 33 /
      data ivalu6(82)   / 1 /
c
      data uopt6(34)    / 'ACTIVITY COEFFICIENT OPTIONS' /
      data uvar6(34)      / 'iopg' /
      data index6(34)    / 1 /
      data udesc6(83)    / 'use B-dot equation' /
      data iopti6(83)   / 34 /
      data ivalu6(83)   / 0 /
      data udesc6(84)    / 'Davies'' equation' /
      data iopti6(84)   / 34 /
      data ivalu6(84)   / -1 /
      data udesc6(85)    / 'Pitzer''s equations' /
      data iopti6(85)   / 34 /
      data ivalu6(85)   / 1 /
c
c**********************************************************************
c          D  E  V  E  L  O  P  M  E  N  T      O  P  T  I  O  N  S   *
c**********************************************************************
c
      data udevl6(1)
     $  / 'check finite difference and Taylor series expression' /
      data udv6vr(1)  / 'iodb' /
      data idev6i(1)  / 10 /
c
      data udevl6(2)
     $  /'check reaction rate finite difference and Taylor series' /
      data udv6vr(2)  / 'iodb' /
      data idev6i(2)  / 11 /
c
c**********************************************************************
c                      T  O  L  E  R  A  N  C  E S                    *
c**********************************************************************
c
      data utol6(1) / 'number of N-R iterations'         /
      data utol6(2) / 'p.r.s. transfer interval'         /
      data utol6(3) / 'residual magnitude'               /
      data utol6(4) / 'correction magnitude'             /
      data utol6(5) / 'search/find tolerance'            /
      data utol6(6) / 'supersaturation '                 /
      data utol6(7) / 'supersaturation set size '        /
      data utol6(8) / "max. size Taylor's series term "  /
      data utol6(9) / 'max. initial value betamx '       /
      data utol6(10) /"max. Taylor's series term (kin.)" /
      data utol6(11) /'corrector iteration'              /
      data utol6(12) /'max. size of N-R correction term' /
      data utol6(13) /'step size (economy mode)'         /
      data utol6(14) /'log mass of phases'               /
      data utol6(15) /'decrement mass (p.r.s.)'          /
      data utol6(16) /'min. left after p.r.s.'           /
      data utol6(17) /'initial step size '               /
      data utol6(18) /'upper limit step size'            /
      data utol6(19) /'maximum order'                    /
      data utol6(20) /'num. attempted assemblages'       /
      data utol6(21) /'slide -> over phase bound.'       /
      data utol6(22) /'slide -> over redox insta.'       /
      data utol6(23) /'fo2 scan control'                 /
c
c end include file x6op7d
c-----------------------------------------------------------------------
