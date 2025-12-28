c eqlj8d.h
c
c     Strings, arrays, and such for the menu-style ("D") input format
c     for INPUT file options exclusive to EQ3NR in Version 8. This
c     INCLUDE file is referenced by EQ3NR and XCON3.
c
c      The relevant arrays are declared in the EQLIB INCLUDE file
c      eqlj8.h.
c
c-----------------------------------------------------------------------
c
c     Jflgi options for EQ3NR:
c
        data ujf3(-1) /'Suppressed      '/
        data ujf3(0)  /'Molality        '/
        data ujf3(1)  /'Molarity        '/
        data ujf3(2)  /'mg/L            '/
        data ujf3(3)  /'mg/kg.sol       '/
        data ujf3(4)  /'ERROR           '/
        data ujf3(5)  /'ERROR           '/
        data ujf3(6)  /'ERROR           '/
        data ujf3(7)  /'Alk., eq/kg.H2O '/
        data ujf3(8)  /'Alk., eq/L      '/
        data ujf3(9)  /'Alk., eq/kg.sol '/
        data ujf3(10) /'Alk., mg/L CaCO3'/
        data ujf3(11) /'Alk., mg/L HCO3-'/
        data ujf3(12) /'ERROR           '/
        data ujf3(13) /'ERROR           '/
        data ujf3(14) /'ERROR           '/
        data ujf3(15) /'ERROR           '/
        data ujf3(16) /'Log activity    '/
        data ujf3(17) /'Log act combo   '/
        data ujf3(18) /'Log mean act    '/
        data ujf3(19) /'pX              '/
        data ujf3(20) /'pH              '/
        data ujf3(21) /'pHCl            '/
        data ujf3(22) /'pmH             '/
        data ujf3(23) /'pmX             '/
        data ujf3(24) /'ERROR           '/
        data ujf3(25) /'Hetero. equil.  '/
        data ujf3(26) /'ERROR           '/
        data ujf3(27) /'Homo. equil.    '/
        data ujf3(28) /'ERROR           '/
        data ujf3(29) /'ERROR           '/
        data ujf3(30) /'Make non-basis  '/
c
c       Note: a distinct, smaller set of jflgi option strings (ujf6)
c       is used by EQ6.
c
c     End of INCLUDE file eqlj8d.h
c-----------------------------------------------------------------------
