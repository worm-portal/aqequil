c eqltf2.h
c
c     This is a pseudo-data file for the titration factors which are
c     required to compute the extended total alkalinity, which is an
c     attempt to match the measured total alkalinity of waters in a
c     range of waters which is more extensive than normal (dilute)
c     ground or surface waters. The set of species considered here is
c     more extensive than the limited set in the pseudo-data file
c     in the EQLIB INCLUDE file eqltf1.h. That set assumes that the
c     the total alkalinity is essentially just HCO3-CO3-OH alkalinity.
c     This set is intended to account for contributions to measured
c     total alkalinity from species such as such as acetate, other
c     organics, borate, silicate, phosphate, and transition metal-
c     hydroxy complexes.
c
c     This pseudo-data file may not be sufficiently complete to
c     represent the non-HCO3-CO3-OH contributions in all waters, nor
c     even all the HCO3-CO3-OH contributions in some kinds of waters.
c     This treatment is intended to be approximate. Not every minor
c     species which might contribute to the measured alkalinity is or
c     should be represented here. The user may wish to extend the range
c     of species considered in this pseudo-data file.
c
c     A species may be known by more than one name in a set of data
c     files. To insure matching against the names on any data file,
c     all likely names of a given species should appear here. For
c     example, the acetate ion may be "Acetate" or "CH3COO-". The
c     dihydrogen borate ion may be "H2BO3-" or its hydration
c     doppelganger "BO2-". The list terminates with an "endit."
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
      integer, parameter :: ntf2a_par = 59
c
      character(len=32), dimension(ntf2a_par) :: utf2xd = (/
     $ '"Acetate             ", 1.  ','"Al(OH)3(aq)         ", 1.  ',
     $ '"Al(OH)4-            ", 2.  ','"BaOH+               ", 1.  ',
     $ '"BO2-                ", 1.  ','"BO3---              ", 3.  ',
     $ '"CaCO3(aq)           ", 2.  ','"CaHCO3+             ", 1.  ',
     $ '"CaHPO4(aq)          ", 1.  ','"CaOH+               ", 1.  ',
     $ '"CaPO4--             ", 2.  ','"CH3CH2COO-          ", 1.  ',
     $ '"CH3COO-             ", 1.  ','"CO3--               ", 2.  ',
     $ '"FeHPO4+             ", 1.  ','"FeHPO4(aq)          ", 1.  ',
     $ '"FeOH+               ", 1.  ','"FePO4--             ", 2.  ',
     $ '"Fe(HS)2(aq)         ", 2.  ','"Fe(HS)3-            ", 3.  ',
     $ '"Fe(OH)2(aq)         ", 2.  ','"Fe(OH)3-            ", 3.  ',
     $ '"Fe(OH)3(aq)         ", 1.  ','"Fe(OH)4-            ", 2.  ',
     $ '"HCO3-               ", 1.  ','"HBO3--              ", 2.  ',
     $ '"HPO4--              ", 1.  ','"HS-                 ", 1.  ',
     $ '"H2BO3-              ", 1.  ','"H2SiO4--            ", 2.  ',
     $ '"H3SiO4-             ", 1.  ','"KCO3-               ", 2.  ',
     $ '"KHCO3(aq)           ", 1.  ','"KHPO4-              ", 1.  ',
     $ '"KOH(aq)             ", 1.  ','"MgCO3(aq)           ", 2.  ',
     $ '"MgHCO3+             ", 1.  ','"MgHPO4(aq)          ", 1.  ',
     $ '"MgOH+               ", 1.  ','"MgPO4--             ", 2.  ',
     $ '"Mg4(OH)4++++        ", 4.  ','"MnCO3(aq)           ", 2.  ',
     $ '"MnHCO3+             ", 1.  ','"MnOH+               ", 1.  ',
     $ '"Mn(OH)3-            ", 3.  ','"NaCO3-              ", 2.  ',
     $ '"NaHCO3(aq)          ", 1.  ','"NaHPO4-             ", 1.  ',
     $ '"NaOH(aq)            ", 1.  ','"NH3(aq)             ", 1.  ',
     $ '"OH-                 ", 1.  ','"PO4---              ", 2.  ',
     $ '"Propanoate          ", 1.  ','"SrCO3(aq)           ", 2.  ',
     $ '"SrHCO3+             ", 1.  ','"SrHPO4(aq)          ", 1.  ',
     $ '"SrOH+               ", 1.  ','"SrPO4--             ", 2.  ',
     $ '"S--                 ", 2.  '/)
c
c     End of INCLUDE file eqltf2.h
c-----------------------------------------------------------------------
