c eqlo8.h
c
c     Strings, arrays, and such for the menu-style ("D") input format
c     for INPUT file options common to EQ3NR and EQ6 in Version 8.
c     INCLUDE file is referenced by EQ6, EQ3NR, XCON3, and XCON6.
c
c       Nxmod alter/supress options:
c
c         ukxm   = array of valid option strings for nxmod options
c
c       Exchange model options for generic ion exchangers:
c
c         ugexmv = array of valid exchange model strings for generic
c                    ion exchangers
c
c       Units options for thermodynamic parameters of generic ion
c         exchangers:
c
c         uxfuni = array of strings for units options for
c                    thermodynamic functions of generic ion
c                    exchangers
c
c       Major option switch (iopt, iopr, iodb, and iopg) arrays:
c
c         jptxpa = maximum number of choices per iopt option switch
c         jprxpa = maximum number of choices per iopr print  option
c                    switch
c         jdbxpa = maximum number of choices per iodb debugging print
c                    option switch
c         jpgxpa = maximum number of choices per iopg activity
c           coefficient option switch
c
c         nptxpa = maximum number of iopt option switches for which
c                  strings are defined
c         nprxpa = maximum number of print option switches for which
c                    strings are defined
c         ndbxpa = maximum number of debugging print option switches
c                    for which strings are defined
c         npgxpa = maximum number of activity coefficient option
c                    switches for which strings are defined
c
c         uopttx = array of title strings for iopt option switches
c         uoprtx = array of title strings for iopr option switches
c         uodbtx = array of title strings for iodb option switches
c         uopgtx = array of title strings for iopg option switches
c
c         uoptcx = array of code relevancy strings for iopt option
c                    switches
c         uoprcx = array of code relevancy strings for iopr option
c                    switches
c         uodbcx = array of code relevancy strings for iodb option
c                    switches
c         uopgcx = array of code relevancy strings for iopg option
c                    switches
c
c         uoptox = array of option strings for iopt option switches
c         uoprox = array of option strings for iopr option switches
c         uodbox = array of option strings for iodb option switches
c         uopgox = array of option strings for iopg option switches
c
c         ioptox = array of option values for iopt option switches
c         ioprox = array of option values for iopr option switches
c         iodbox = array of option values for iodb option switches
c         iopgox = array of option values for iopg option switches
c
c     The strings are defined in the EQLIB include file eqlo8d.h.
c
c-----------------------------------------------------------------------
c
      integer jptxpa,jprxpa,jdbxpa,jpgxpa,nptxpa,nprxpa,ndbxpa,npgxpa,
     $ nvet_par
      parameter (jptxpa = 5,jprxpa = 5,jdbxpa = 5,jpgxpa = 5,
     $ nptxpa = 20,nprxpa = 20,ndbxpa = 20,npgxpa = 20,nvet_par = 10)
c
      integer ioptox,ioprox,iodbox,iopgox
c
      common /eqlo8i/ ioptox(jptxpa,nptxpa),ioprox(jprxpa,nprxpa),
     $ iodbox(jdbxpa,ndbxpa),iopgox(jpgxpa,npgxpa)
c
      character*72 uopttx,uoprtx,uodbtx,uopgtx
      character*72 uoptox,uoprox,uodbox,uopgox
      character*24 uoptcx,uoprcx,uodbcx,uopgcx
      character*24 ugexmv
      character*16 ukxm
      character*8 uxfuni
c
      common /eqlo8c/ uopttx(nptxpa),uoprtx(nprxpa),uodbtx(ndbxpa),
     $ uopgtx(npgxpa),uoptox(jptxpa,nptxpa),uoprox(jprxpa,nprxpa),
     $ uodbox(jdbxpa,ndbxpa),uopgox(jpgxpa,npgxpa),uoptcx(nptxpa),
     $ uoprcx(nprxpa),uodbcx(ndbxpa),uopgcx(npgxpa),ugexmv(nvet_par),
     $ ukxm(-1:2),uxfuni(4)
c
c     End of INCLUDE file eqlo8.h
c-----------------------------------------------------------------------
