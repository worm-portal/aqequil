      subroutine cdakey(iopg,nopgmx,noutpt,nttyo,udakey,udatfi)
c
c     This subroutine checks the flag string "udakey" read from the data
c     file and checks it against a key-list to test whether or not
c     the data file used is consistent with the specified model for
c     aqueous species activity coefficients (determined by iopg(1)).
c     The string "udatfi" identifies the data file. The legal
c     combinations are:
c
c           Option                 iopg(1)      Legal key
c
c        Davies equation             -1          SEDH
c        B-dot equation               0          SEDH
c        Pitzer's equations           1          Pitzer
c        HC + DHC equations           2          SEDH
c
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c      Principal input:
c
c
c      Principal output:
c
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nopgmx
c
      integer noutpt,nttyo
c
      integer iopg(nopgmx)
c
      character*8 udakey,udatfi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3,j4,j5
c
      integer ilnobl
c
      character*8 uiakey,usedh,upitz,ux8
c
c-----------------------------------------------------------------------
c
      data usedh  /'SEDH    '/,upitz/'Pitzer  '/
c
c-----------------------------------------------------------------
c
c     Set the correct key.
c
      if (iopg(1) .eq. -1) then
        uiakey = usedh
      elseif (iopg(1) .eq. 0) then
        uiakey = usedh
      elseif (iopg(1) .eq. 1) then
        uiakey = upitz
      elseif (iopg(1) .eq. 2) then
        uiakey = usedh
      else
        write (ux8,'(i5)') iopg(1)
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1000) ux8(1:j2)
        write (nttyo,1000) ux8(1:j2)
 1000   format(/' * Error - (EQLIB/cdakey) The iopg(1) value of ',a,
     $  ' specified on',/7x,"the input file doesn't correspond to an",
     $  ' available model for the',/7x,'activity coefficients of',
     $  ' the aqueous species.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (udakey(1:8) .ne. uiakey(1:8)) then
c
c       Error, bad combination.
c
        write (ux8,'(i5)') iopg(1)
        call lejust(ux8)
        j2 = ilnobl(ux8)
        j3 = ilnobl(udatfi)
        j4 = ilnobl(uiakey)
        j5 = ilnobl(udakey)
        write (noutpt,1010) ux8(1:j2),udatfi(1:j3),uiakey(1:j4),
     $  udakey(1:j5)
        write (nttyo,1010) ux8(1:j2),udatfi(1:j3),uiakey(1:j4),
     $  udakey(1:j5)
 1010   format(/' * Error - (EQLIB/cdakey) The iopg(1) value of ',a,
     $  ' specified on',/7x,"the input file doesn't correspond to the",
     $  ' type of model for the',/7x,'activity coefficients of the',
     $  ' aqueous species required to use the',/7x,a,' data file.',
     $  ' The iopg1(1) value calls for using a ',a,' type model,',
     $  /7x,'while the data file is consistent with a ',a,' type',
     $  ' model. Change',/7x,'the iopg(1) value on the input file',
     $  ' or specify a compatible data file.')
        stop
      endif
c
      end
