c eqltf1.h
c
c     This is a pseudo-data file for the titration factors which are
c     required to compute the HCO3-CO3-OH total alkalinity in normal
c     (e.g., dilute) groundwaters. The species considered here are
c     the three anions and their 1:1 ion pairs with Na+, K+, Mg++,
c     and Ca++. Addition of higher order ion pairs, such as
c     Mg4(OH)4++++ or Ca(OH)2(aq), is not recommended, as these
c     are not signficant in dilute groundwaters. Do not add other
c     kinds of species which may contribute to total alkalinity,
c     such as acetate, other organics, borate, silicate, phosphate,
c     or transition metal-hydroxy complexes. Instead, add any such
c     needed species to the pseudo-data file in the companion EQLIB
c     INCLUDE file eqltf2.h.
c
c     Note that the list of species terminates with an "endit."
c
c     This pseudo-data file is referenced by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
      integer, parameter :: ntf1a_par = 15
c
      character(len=32), dimension(ntf1a_par) :: utf1xd = (/
     $ '"CaCO3(aq)           ", 2.  ',
     $ '"CaHCO3+             ", 1.  ',
     $ '"CaOH+               ", 1.  ',
     $ '"CO3--               ", 2.  ',
     $ '"HCO3-               ", 1.  ',
     $ '"KCO3-               ", 2.  ',
     $ '"KHCO3(aq)           ", 1.  ',
     $ '"KOH(aq)             ", 1.  ',
     $ '"MgCO3(aq)           ", 2.  ',
     $ '"MgHCO3+             ", 1.  ',
     $ '"MgOH+               ", 1.  ',
     $ '"NaCO3-              ", 2.  ',
     $ '"NaHCO3(aq)          ", 1.  ',
     $ '"NaOH(aq)            ", 1.  ',
     $ '"OH-                 ", 1.  '/)
c
c     End of INCLUDE file eqltf1.h
c-----------------------------------------------------------------------
