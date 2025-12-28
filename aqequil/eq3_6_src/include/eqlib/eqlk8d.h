c eqlk8d.h
c
c     Strings, arrays, and such for the menu-style ("D") input format
c     for INPUT file options exclusive to EQ6 in Version 8. This INCLUDE
c     INCLUDE file is referenced by EQ6, EQ3NR, and XCON6.
c
c     The relevant arrays are declared in the EQLIB INCLUDE file
c     eqlk8.h.
c
c-----------------------------------------------------------------------
c
c     Reactant type strings.
c
        data urcjco(0) /'Pure mineral            '/
        data urcjco(1) /'Solid solution          '/
        data urcjco(2) /'Special reactant        '/
        data urcjco(3) /'Aqueous species         '/
        data urcjco(4) /'Gas species             '/
        data urcjco(5) /'Generic ion exchanger   '/
c
c-----------------------------------------------------------------------
c
c     Reactant status strings.
c
        data urcjre(-1) /'Saturated, reacting     '/
        data urcjre(0)  /'Reacting                '/
        data urcjre(1)  /'Exhausted               '/
        data urcjre(2)  /'Saturated, not reacting '/
c
c-----------------------------------------------------------------------
c
c     Rate law strings.
c
        data urcnrk(-1,1) /'Use backward rate law   '/
        data urcnrk(0,1)  /'Illegal value           '/
        data urcnrk(1,1)  /'Relative rate equation  '/
        data urcnrk(2,1)  /'TST rate equation       '/
        data urcnrk(3,1)  /'Linear rate equation    '/
        data urcnrk(-1,2) /'Use forward rate law    '/
        data urcnrk(0,2)  /'Partial equilibrium     '/
        data urcnrk(1,2)  /'Relative rate equation  '/
        data urcnrk(2,2)  /'TST rate equation       '/
        data urcnrk(3,2)  /'Linear rate equation    '/
c
c-----------------------------------------------------------------------
c
c     Rate direction header strings.
c
        data urcrld(1) /'Forward rate law        '/
        data urcrld(2) /'Backward rate law       '/
c
c-----------------------------------------------------------------------
c
c     Rate direction sign strings.
c
        data urcrls(1) /'+'/
        data urcrls(2) /'-'/
c
c-----------------------------------------------------------------------
c
c     Relative rate header strings.
c
        data urcrel(1,1) /'dXi(n)/dXi (mol/mol)    '/
        data urcrel(2,1) /'d2Xi(n)/dXi2 (mol/mol2) '/
        data urcrel(3,1) /'d3Xi(n)/dXi3 (mol/mol3) '/
        data urcrel(1,2) /'-dXi(n)/dXi (mol/mol)   '/
        data urcrel(2,2) /'-d2Xi(n)/dXi2 (mol/mol2)'/
        data urcrel(3,2) /'-d3Xi(n)/dXi3 (mol/mol3)'/
c
c-----------------------------------------------------------------------
c
c       Nxopt mineral subset-selection suppression option strings:
c
        data uxopti(1) /'None    '/
        data uxopti(2) /'All     '/
        data uxopti(3) /'Alwith  '/
        data uxopti(4) /'Allwith '/
c
c-----------------------------------------------------------------------
c
c       Jflgi strings for EQ6.
c
        data ujf6(1) /'Moles           '/
        data ujf6(2) /'Make non-basis  '/
c
c       Note: 'Moles" -> jflgi = 0, 'Make non-basis' -> jflgi = 30.
c       A distinct, larger set of jflgi option strings (ujf3) is used
c       by EQ3NR.
c
c     End of INCLUDE file eqlk8d.h
c-----------------------------------------------------------------------
