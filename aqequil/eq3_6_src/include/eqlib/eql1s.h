c eql1s.h
c
c     The presssure on the 1.013-bar/steam-saturation curve. The
c     variable is presh. The coefficients for computing it are in the
c     apresh array, which has fixed dimensions (apresh(5,2)). The values
c     of the coefficients are stored in a DATA statement in the EQLIB
c     INCLUDE file eql1sd.h.
c
c     Note that presg and apresg (which are in EQLIB/eqlxf.h) are,
c     respectively, the pressure and the corresponding coefficients) for
c     the data file reference pressure curve. This may or may not be the
c     1.013-bar/steam-saturation curve. Even if it is, the two
c     representations may not exactly coincide.
c
      integer narxth
c
      common /eql1sa/ narxth(2)
c
      real*8 apresh
c
      common /eql1sb/ apresh(5,2)
c
      real*8 presh
c
      common /eql1sc/ presh
c
c     End of INCLUDE file eql1s.h
c-----------------------------------------------------------------------
