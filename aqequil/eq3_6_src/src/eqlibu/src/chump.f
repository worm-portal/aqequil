      subroutine chump(eps,irang,noutpt,nttyo)
c
c     This subroutine tests the floating point precision and exponent
c     range for adequacy for EQ3/6 calculations. This is essentially
c     a trap to insure that most floating point variables are real*8
c     and, if the computer is a VAX, that the code has been compiled
c     with the 'G_FLOATING' option (which increases the exponent range).
c
c     The machine epsilon and exponent range are respectively tested
c     against epstst and irgtst. The former should be no greater than
c     1.e-14. The latter should be at least 100, and the recommended
c     value is 300.
c
c     This subroutine is called by:
c
c        EQLIBU/flpars.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       noutpt = the unit number of the output file
c       nttyo  = the unit number of the screen file
c
c     Output:
c
c       eps    = the real*8 machine epsilon
c       irang  = the real*8 exponent range
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer irang,noutpt,nttyo
c
      real*8 eps
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer irgtst
c
      logical qwarn
c
      real*8 epstst
c
c-----------------------------------------------------------------------
c
c     WARNING: Under no circumstances should the value of irgtst be
c     changed to a value less than 100. If decreased to from the
c     recommended value of 300, make a corresponding adjustment of the
c     overflow truncation limit in EQLIBU/texp.f or you will lose the
c     floating-point overflow protection built into this software.
c
      data irgtst /300/
c
      data epstst /1.e-14/
c
c-----------------------------------------------------------------------
c
      qwarn = .false.
      if (eps .gt. epstst) then
        write (noutpt,1007) eps,epstst
        write (nttyo,1007) eps,epstst
 1007   format(/' * Warning - (EQLIBU/chump) Have insufficient',
     $  ' floating point',/7x,'epsilon= ',1pe12.5,'. This parameter',
     $  ' should be no larger',/7x,'than ',1pe12.5,'. You need a',
     $  ' computer with better floating',/7x,'point characteristics',
     $  ' to run this software.')
        qwarn=.true.
      endif
c
      if (irang .lt. irgtst) then
        write (noutpt,1008) irang,irgtst
        write (nttyo,1008) irang,irgtst
 1008   format(/' * Warning - (EQLIBU/chump) Have insufficient',
     $  ' floating point',/7x,'exponent range= +/- ',i4,'. This should',
     $  ' be at least +/- ',i4,/7x,'to run this software. You need',
     $  ' a computer with better',/7x,'floating point characteristics.',
     $  ' Alternatively, you need to recompile',/7x,'with a special',
     $  ' compiler option.')
        if (irang.ge.37 .and. irang.le.39) then
          write (noutpt,1010)
          write (nttyo,1010)
 1010     format(/' * Note - (EQLIBU/chump) This computer appears to',
     $    ' be a VAX.',/7x,'If you run under VMS, try recompiling with',
     $    /7x,'the G_FLOATING option.',/7x,'This will increase the',
     $    ' exponent range to +/- 308. Some older',/7x,'VAXes do not',
     $    ' have the necessary hardware to implement this option.',
     $    /7x,'Also, there appears to be no equivalent to the',
     $    ' G_FLOATING option',/7x,'on VAXes running the ULTRIX',
     $    ' operating system. If you are',/7x,'running a VAX with',
     $    ' ULTRIX, activate the call to the system subroutine TRPFPE',
     $    /7x,'in this subroutine (EQLIBU/chump.f). This call is',
     $    ' normally commented.',/7x,'out. You might wish to consider',
     $    ' getting a computer which has better',/7x,'floating point',
     $    ' characteristics.')
        endif
        qwarn=.true.
      endif
c
      if (qwarn) then
        write (noutpt,1020)
        write (nttyo,1020)
 1020   format(/' * Warning - (EQLIBU/chump) Will continue execution.',
     $  ' You may encounter',/7x,'floating point overflows, especially',
     $  ' if you are running EQ3NR',/7x,'or EQ6. You may be able to',
     $  ' make some EQ3NR and EQ6 problems run',/7x,'by specifying',
     $  ' judicious basis switches on the input file.',/)
      endif
c
c-----------------------------------------------------------------------
c
com   BEGIN_VAX_DEPENDENT_CODING
com
com     The following VAX system call traps and repairs floating point
com     overflows when running under the ULTRIX operating system.
com     This call may have be to be activated (by de-commenting it) if
com     you are running on a VAX with ULTRIX. On some older VAXES, you
com     may have to substitute "TRAPOV" for "TRPFPE".
com
c       call trpfpe(0,1.d+30)
com
com   END_VAX_DEPENDENT_CODING
c
c-----------------------------------------------------------------------
c
      end
