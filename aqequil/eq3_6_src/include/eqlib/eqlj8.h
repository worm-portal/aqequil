c eqlj8.h
c
c     Strings, arrays, and such for the menu-style ("D") input format
c     for INPUT file options exclusive to EQ3NR in Version 8. This
c     INCLUDE file is referenced by EQ3NR and XCON3.
c
c       Jflgi options:
c
c         njfxpa = the largest jflgi value
c         ujf3   = array of strings corresponding to jflgi values
c
c       Note: a distinct, smaller set of jflgi option strings (ujf6)
c       is used by EQ6.
c
c     The strings are defined in the EQLIB include file eqlj8d.h.
c
c-----------------------------------------------------------------------
c
      integer njfxpa
      parameter (njfxpa = 30)
c
      character*16 ujf3
c
      common /eqlj8c/ ujf3(-1:njfxpa)
c
c     End of INCLUDE file eqlj8.h
c-----------------------------------------------------------------------
