c x6op7.h
c
c     Strings, arrays, and such for the menu-style ("D") input
c     format for version 7.
c
c     Options (iopt, iopg, iopr, and iodb):
c
c       nop6pa = maximum number of options (the sum of all iopt,
c                  iopg, iopr, and iodb options)
c       nod6pa = maximum number of choices for all options
c                  in the above group
c       uopt6  = array of strings describing the above options
c       udesc6 = array of strings describing option choices
c       uvar6  = array of names of option arrays (iopt, iopg, iopr, or
c                  iodb) corresponding to the options
c       iopti6 = array of indices of option array elements corresponding
c                  to strings describing option choices
c       ivalu6 = array of values of option array elements corresponding
c                  to strings describing option choices
c       index6 = pointer array giving the index of the option
c                  corresponding to a string describing an option choice
c
c     Development options (iopt, iopg, iopr, and iodb):
c
c       ndv6pa = maximum number of development options
c       udevl6 = array of strings describing development options
c       udv6vr = array of names of option arrays (iopt, iopg, iopr, or
c                  iodb) corresponding to the development options
c       idev6i = array of values of option array elements corresponding
c                  to strings describing development option choices
c
c     Tolerances:
c
c       nto6pa = maximum number of tolerance variables
c       utol6  = array of tolerance descriptor strings
c
c-----------------------------------------------------------------------
c
      integer nop6pa,nod6pa,ndv6pa,nto6pa
      parameter (nop6pa = 34,nod6pa = 85,ndv6pa = 2,nto6pa = 23)
c
      integer iopti6,ivalu6,index6,idev6i
c
      common /x6op7/ iopti6(nod6pa),ivalu6(nod6pa),index6(nop6pa),
     $ idev6i(ndv6pa)
c
      character*80 udesc6,uopt6,udevl6
      character*32 utol6
      character*8 uvar6,udv6vr
c
      common /x6op7c/ udesc6(nod6pa),udevl6(ndv6pa),udv6vr(ndv6pa),
     $ utol6(nto6pa),uopt6(nop6pa),uvar6(nop6pa)
c
c end include file x6op7
c-----------------------------------------------------------------------
