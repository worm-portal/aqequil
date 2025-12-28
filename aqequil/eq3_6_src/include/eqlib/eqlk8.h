c eqlk8.h
c
c     Strings, arrays, and such for the menu-style ("D") input format
c     for INPUT file options exclusive to EQ6 in Version 8. This INCLUDE
c     file is referenced by EQ6, EQ3NR, and XCON6.
c
c       Reactant type strings:
c
c         urcjco = array of strings corresponding to jcode reactant
c                    type flags
c
c       Reactant status strings:
c
c         urcjre = array of strings corresponding to jreac reactant
c                    status flags
c
c       Rate law strings:
c
c         urcnrk = array of strings corresponding to nrk rate law codes
c
c       Rate direction header strings:
c
c         urcrld = array of strings corresponding to rate law directions
c                    (forward and backward)
c
c       Rate direction sign strings:
c
c         urcrls = array of strings corresponding to signs of rate law
c                    directions (+ and -)
c
c       Relative rate header strings:
c
c         urcrel = array of strings corresponding to relative rates of
c                    various orders
c
c       Nxopt mineral subset-selection suppression option strings:
c
c         uxopti = array of strings corresponding to valid options
c
c       Jflgi strings for EQ6:
c
c         ujf6   = array of jflgi option strings for EQ6
c
c       Note: a distinct, larger set of jflgi option strings (ujf3)
c       is used by EQ3NR.
c
c     The strings are defined in the EQLIB include file eqlk8d.h.
c
c-----------------------------------------------------------------------
c
      character*24 urcrld,urcjco,urcjre,urcnrk,urcrel
      character*16 ujf6
      character*16 uxopti
      character*2 urcrls
c
      common /eqlk8c/ urcrld(2),urcjco(0:5),urcjre(-1:2),
     $ urcnrk(-1:4,2),urcrel(3,2),ujf6(2),uxopti(4),urcrls(2)
c
c     End of INCLUDE file eqlk8.h
c-----------------------------------------------------------------------
