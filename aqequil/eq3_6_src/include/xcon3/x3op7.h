c x3op7.h
c
c     Strings, arrays, and such for the menu-style ("D") input
c     format for version 7.
c
c     Jflag strings:
c
c       njf7pa = the largest jflag value, including "hidden" values
c                  (> 30) that map to standard ones (0-30)
c       ujflg7 = array of strings corresponding to jflag values
c
c     Options (iopt, iopg, and iopr):
c
c       nop3pa = maximum number of options (the sum of all iopt,
c                  iopg, and iopr options)
c       nod3pa = maximum number of choices for all options in
c                  the above group
c       uopt3  = array of strings describing the above options
c       udesc3 = array of strings describing option choices
c       uvar3  = array of names of option arrays (iopt, iopg, or
c                  iopr) corresponding to the options
c       iopti3 = index of the option array element corresponding to
c                  a string describing an option choice
c       ivalu3 = array of values of option array elements corresponding
c                  to strings describing option choices
c       index3 = pointer array giving the index of the option
c                  corresponding to a string describing an option choice
c
c     Debug options (iodb):
c
c       ndb3pa = maximum number of debug options
c       udebug = array of strings describing debug options
c       idbugi = array of indices of the iodb array corresponding to
c                  strings describing debug options
c
c     Development options (iopt, iopg, iopr, and iodb):
c
c       None
c
c     Tolerances:
c
c       nto3pa = maximum number of tolerance strings
c       utol3  = array of strings describing tolerance variables
c
c-----------------------------------------------------------------------
c
      integer njf7pa,nop3pa,nod3pa,ndb3pa,nto3pa
      parameter (njf7pa = 32,nop3pa = 15,nod3pa = 40,ndb3pa = 8,
     $ nto3pa = 4)
c
      integer iopti3,ivalu3,index3,idbugi
c
      common /x3op7/ iopti3(nod3pa),ivalu3(nod3pa),index3(nop3pa),
     $ idbugi(ndb3pa)
c
      character*80 uopt3,udesc3,udebug
      character*40 ujflg7
      character*32 utol3
      character*8 uvar3
c
      common /x3op7c/ uopt3(nop3pa),udesc3(nod3pa),udebug(ndb3pa),
     $ ujflg7(-1:njf7pa),utol3(nto3pa),uvar3(nop3pa)
c
c     End of include file x3op7.h
c-----------------------------------------------------------------------
