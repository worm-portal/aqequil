c x6sbsp.h
c
c     This is a pseudo-data file for arrays which can be used to map
c     the names of chemical element names to strict basis species.
c     It is intended to be used in mapping certain data between
c     EQ6 input files that are version level less than 8.0 to those
c     that are 8.0 or above, and vice versa. One function is to
c     construct reactions for special reactants for input files that
c     are version level 8.0 or above, as ony the compositions of such
c     reactants are given on input files of version level less than 8.0.
c
c     The element name and corresponding basis species name are
c     followed by the number of moles of the element per mole of
c     the species, the number of moles of H per mole of species
c     (zero for H/H+), the number of moles of O per mole of species
c     (zero for O/H2O), and the electrical charge number of the
c     species.
c
c     Note that the list of species terminates with an "endit."
c
c     This pseudo-data file is referenced by:
c
c       XCON6/xcon6.f
c
c-----------------------------------------------------------------------
c
      data (uelnam(n),uelspn(n),celspe(n),celsph(n),celspo(n),zelsp(n),
     $  n = 1,10) /
     $ "Ag      ", "Ag+                     ",1.,0.,0.,+1.,
     $ "Al      ", "Al+++                   ",1.,0.,0.,+3.,
     $ "Am      ", "Am+++                   ",1.,0.,0.,+3.,
     $ "Ar      ", "Ar(aq)                  ",1.,0.,0.,0.,
     $ "As      ", "H2AsO4-                 ",1.,2.,4.,-1.,
     $ "Au      ", "Au+                     ",1.,0.,0.,+1.,
     $ "B       ", "B(OH)3(aq)              ",1.,3.,3.,0.,
     $ "Ba      ", "Ba++                    ",1.,0.,0.,+2.,
     $ "Be      ", "Be++                    ",1.,0.,0.,+2.,
     $ "Br      ", "Br-                     ",1.,0.,0.,-1./
c
      data (uelnam(n),uelspn(n),celspe(n),celsph(n),celspo(n),zelsp(n),
     $  n = 11,20) /
     $ "C       ", "HCO3-                   ",1.,1.,3.,-1.,
     $ "Ca      ", "Ca++                    ",1.,0.,0.,+2.,
     $ "Cd      ", "Cd++                    ",1.,0.,0.,+2.,
     $ "Ce      ", "Ce+++                   ",1.,0.,0.,+3.,
     $ "Cl      ", "Cl-                     ",1.,0.,0.,-1.,
     $ "Co      ", "Co++                    ",1.,0.,0.,+2.,
     $ "Cr      ", "CrO4--                  ",1.,0.,4.,-2.,
     $ "Cs      ", "Cs+                     ",1.,0.,0.,+1.,
     $ "Cu      ", "Cu++                    ",1.,0.,0.,+2.,
     $ "Dy      ", "Dy+++                   ",1.,0.,0.,+3./
c
      data (uelnam(n),uelspn(n),celspe(n),celsph(n),celspo(n),zelsp(n),
     $  n = 21,30) /
     $ "Er      ", "Er+++                   ",1.,0.,0.,+3.,
     $ "Eu      ", "Eu+++                   ",1.,0.,0.,+3.,
     $ "F       ", "F-                      ",1.,0.,0.,-1.,
     $ "Fe      ", "Fe++                    ",1.,0.,0.,+2.,
     $ "Ga      ", "Ga+++                   ",1.,0.,0.,+3.,
     $ "Gd      ", "Gd+++                   ",1.,0.,0.,+3.,
     $ "H       ", "H+                      ",1.,0.,0.,+1.,
     $ "He      ", "He(aq)                  ",1.,0.,0.,0.,
     $ "Hg      ", "Hg++                    ",1.,0.,0.,+2.,
     $ "Ho      ", "Ho+++                   ",1.,0.,0.,+3./
c
      data (uelnam(n),uelspn(n),celspe(n),celsph(n),celspo(n),zelsp(n),
     $  n = 31,40) /
     $ "I       ", "I-                      ",1.,0.,0.,-1.,
     $ "In      ", "In+++                   ",1.,0.,0.,+3.,
     $ "K       ", "K+                      ",1.,0.,0.,+1.,
     $ "Kr      ", "Kr(aq)                  ",1.,0.,0.,0.,
     $ "La      ", "La+++                   ",1.,0.,0.,+3.,
     $ "Li      ", "Li+                     ",1.,0.,0.,+1.,
     $ "Lu      ", "Lu+++                   ",1.,0.,0.,+3.,
     $ "Mg      ", "Mg++                    ",1.,0.,0.,+2.,
     $ "Mn      ", "Mn++                    ",1.,0.,0.,+2.,
     $ "Mo      ", "MoO4--                  ",1.,0.,4.,-2./
c
      data (uelnam(n),uelspn(n),celspe(n),celsph(n),celspo(n),zelsp(n),
     $  n = 41,50) /
cXX
cXX  $ "N       ", "NO3-                    ",1.,0.,3.,-1.,
     $ "N       ", "NO3-                    ",1.,3.,0.,0.,
cXX
     $ "Na      ", "Na+                     ",1.,0.,0.,+1.,
     $ "Nd      ", "Nd+++                   ",1.,0.,0.,+3.,
     $ "Ne      ", "Ne(aq)                  ",1.,0.,0.,0.,
     $ "Ni      ", "Ni++                    ",1.,0.,0.,+2.,
     $ "Np      ", "Np++++                  ",1.,0.,0.,+4.,
     $ "O       ", "H2O                     ",1.,2.,0.,0.,
     $ "P       ", "HPO4--                  ",1.,1.,4.,-2.,
     $ "Pb      ", "Pb++                    ",1.,0.,0.,+2.,
     $ "Pd      ", "Pd++                    ",1.,0.,0.,+2./
c
      data (uelnam(n),uelspn(n),celspe(n),celsph(n),celspo(n),zelsp(n),
     $  n = 51,60) /
     $ "Pr      ", "Pr+++                   ",1.,0.,0.,+3.,
     $ "Pu      ", "Pu++++                  ",1.,0.,0.,+4.,
     $ "Ra      ", "Ra++                    ",1.,0.,0.,+2.,
     $ "Rb      ", "Rb+                     ",1.,0.,0.,+1.,
     $ "Re      ", "ReO4-                   ",1.,0.,4.,-1.,
     $ "Rn      ", "Rn(aq)                  ",1.,0.,0.,0.,
     $ "Ru      ", "RuO4--                  ",1.,0.,4.,-2.,
     $ "S       ", "SO4--                   ",1.,0.,4.,-2.,
     $ "Sb      ", "Sb(OH)3(aq)             ",1.,3.,3.,0.,
     $ "Sc      ", "Sc+++                   ",1.,0.,0.,+3./
c
      data (uelnam(n),uelspn(n),celspe(n),celsph(n),celspo(n),zelsp(n),
     $  n = 61,70) /
     $ "Se      ", "SeO3--                  ",1.,0.,3.,-2.,
     $ "Si      ", "SiO2(aq)                ",1.,0.,2.,0.,
     $ "Sm      ", "Sm+++                   ",1.,0.,0.,+3.,
     $ "Sn      ", "Sn++                    ",1.,0.,0.,+2.,
     $ "Sr      ", "Sr++                    ",1.,0.,0.,+2.,
     $ "Tb      ", "Tb+++                   ",1.,0.,0.,+3.,
     $ "Tc      ", "TcO4-                   ",1.,0.,4.,-1.,
     $ "Th      ", "Th++++                  ",1.,0.,0.,+4.,
     $ "Ti      ", "Ti(OH)4(aq)             ",1.,4.,4.,0.,
     $ "Tl      ", "Tl+                     ",1.,0.,0.,+1./
c
      data (uelnam(n),uelspn(n),celspe(n),celsph(n),celspo(n),zelsp(n),
     $  n = 71,80) /
     $ "Tm      ", "Tm+++                   ",1.,0.,0.,+3.,
     $ "U       ", "UO2++                   ",1.,0.,2.,+2.,
     $ "V       ", "VO++                    ",1.,0.,1.,+2.,
     $ "W       ", "WO4--                   ",1.,0.,4.,-2.,
     $ "Xe      ", "Xe(aq)                  ",1.,0.,0.,0.,
     $ "Y       ", "Y+++                    ",1.,0.,0.,+3.,
     $ "Yb      ", "Yb+++                   ",1.,0.,0.,+3.,
     $ "Zn      ", "Zn++                    ",1.,0.,0.,+2.,
     $ "Zr      ", "Zr(OH)2++               ",1.,2.,2.,+2.,
     $ "endit.  ", "                        ",0.,0.,0.,0. /
c
c end include file x6sbsp
c-----------------------------------------------------------------------
