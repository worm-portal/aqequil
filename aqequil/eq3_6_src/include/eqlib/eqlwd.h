c eqlwd.h
c
c     This constants needed for the WIPP brine density model:
c
c       density (g/L) = a + b x TDS (g/L)
c
c     This is based on data for pure NaCl solutions at 20C.
c
c     The respective variables corresponding to a and b are adwipp
c     and bdwipp. The values of the coefficients are stored in a
c     DATA statement in the EQLIB INCLUDE file eqlwdd.h.
c
      real(8) adwipp,bdwipp
c
      common /eqlwda/ adwipp,bdwipp
c
c     End of INCLUDE file eqlwd.h
c-----------------------------------------------------------------------
