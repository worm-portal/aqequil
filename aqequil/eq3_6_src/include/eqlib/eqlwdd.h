c eqlwdd.h
c
c     This is a pseudo-data file for constants needed for the
c     WIPP brine density model:
c
c       density (g/L) = a + b x TDS (g/L)
c
c     This is based on data for pure NaCl solutions at 20C.
c
c     The adwipp (a) and bdwipp (b) variables are declared in
c     EQLIB/eqlwd.h.
c
c     This pseudo-data file is referenced by:
c
c       EQ3NR/bkdeq3.f
c       EQ6/bkdeq6.f
c
      data adwipp /1000.96/, bdwipp /0.639963/
c
c     End of INCLUDE file eqlwdd.h
c-----------------------------------------------------------------------
