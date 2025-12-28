c eql1sd.h
c
c     This is a pseudo-data file for polynomial coefficients for
c     computing the pressure on the 1.013-bar/steam-saturation
c     curve. These data were obtained using the old standard EQ3/6
c     pressure grid. They therefore are for the case of two temperature
c     ranges (0-100C and 100-300C), with a maximum of five coefficients
c     per range.
c
c     The apresh array is declared in EQLIB/eql1s.h.
c
c     This pseudo-data file is referenced by:
c
c       EQ3NR/bkdeq3.f
c       EQ6/bkdeq6.f
c
c-----------------------------------------------------------------------
c
      data narxth(1) /1/, narxth(2) /5/
c
      data apresh(1,1) / 1.013200000E+00/,
     $     apresh(2,1) / 0./,apresh(3,1) / 0./,
     $     apresh(4,1) / 0./,apresh(5,1) / 0./,
     $     apresh(1,2) /-4.345000000E-01/,
     $     apresh(2,2) / 7.632333333E-03/,
     $     apresh(3,2) / 5.514000000E-05/,
     $     apresh(4,2) /-1.263733333E-06/,
     $     apresh(5,2) / 1.396800000E-08/
c
c     End of INCLUDE file eql1sd.h
c-----------------------------------------------------------------------
