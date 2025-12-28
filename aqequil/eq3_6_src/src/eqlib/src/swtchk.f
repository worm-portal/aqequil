      subroutine swtchk(cdrs,jflag,jsflag,nbaspx,nbt,nbtmax,ndrs,ndrsmx,
     $ ndrsr,noutpt,ns1,ns2,nstmax,nttyo,uspec)
c
c     This subroutine checks a proposed basis switch for EQLIB/switch.f
c     to see if doing the switch is ok or not.
c
c     This subroutine is called by:
c
c       EQLIB/switch.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
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
      integer nbtmax,ndrsmx,nstmax
c
      integer noutpt,nttyo
c
      integer jflag(nstmax),jsflag(nstmax),nbaspx(nbtmax),
     $ ndrs(ndrsmx),ndrsr(2,nstmax)
      integer nbt,ns1,ns2
c
      character*48 uspec(nstmax)
c
      real*8 cdrs(ndrsmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen1,jlen2,nb2,ncmpt2,nerr
c
      integer nbasis
c
      character*56 usp156,usp256
c
      real*8 cx12
c
      real*8 coefdr
c
c-----------------------------------------------------------------------
c
      nerr = 0
c
c     Calling sequence substitutions:
c       jlen1 for jlen
c       uspec(ns1) for unam48
c       usp156 for uspn56
c
      call fmspnx(jlen1,uspec(ns1),usp156)
c
c     Calling sequence substitutions:
c       jlen2 for jlen
c       uspec(ns2) for unam48
c       usp256 for uspn56
c
      call fmspnx(jlen2,uspec(ns2),usp256)
c
      if (ns1 .eq. ns2) then
        write (noutpt,1000) usp156(1:jlen1)
        write (nttyo,1000) usp156(1:jlen1)
 1000   format(/" * Error - (EQLIB/swtchk) Can't replace the species ",
     $  a,/7x,'with itself in the basis set.')
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make sure that the first species is not suppressed
c
      if (jsflag(ns1) .gt. 0) then
        write (noutpt,1010) usp156(1:jlen1),usp256(1:jlen2)
        write (nttyo,1010) usp156(1:jlen1),usp256(1:jlen2)
 1010   format(/" * Error - (EQLIB/swtchk) Can't replace the species ",
     $  a,/7x,'in the basis set with ',a,', because the former is',
     $  ' suppressed.')
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make sure that the second species is not suppressed
c
      if (jsflag(ns2) .gt. 0) then
        write (noutpt,1020) usp156(1:jlen1),usp256(1:jlen2)
        write (nttyo,1020) usp156(1:jlen1),usp256(1:jlen2)
 1020   format(/" * Error - (EQLIB/swtchk) Can't replace the species ",
     $  a,/7x,'in the basis set with ',a,', because the latter is',
     $  ' suppressed.')
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the linking reaction.
c
      ncmpt2 = ndrsr(2,ns2) - ndrsr(1,ns2) + 1
      if (ncmpt2 .lt. 2) then
c
c       The second species must not be in the strict basis, because
c       there is then no possibility of a linking reaction.
c
        write (noutpt,1030) usp156(1:jlen1),usp256(1:jlen2)
        write (nttyo,1030) usp156(1:jlen1),usp256(1:jlen2)
 1030   format(/" * Error - (EQLIB/swtchk) Can't replace the species ",
     $  a,/7x,'in the basis set with ',a,', because the latter is',
     $  ' already in',/7x,'the strict basis set, so there is no',
     $  ' possibility of a valid',/7x,'linking reaction.')
      else
c
c       Make sure that the first species appears as a product in the
c       reaction belonging to the second species.
c
c       Calling sequence substitutions:
c         ns1 for nse
c         ns2 for ns
c
        cx12 = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,ns2,nstmax)
        if (cx12 .eq. 0.) then
        write (noutpt,1040) usp156(1:jlen1),usp256(1:jlen2)
        write (nttyo,1040) usp156(1:jlen1),usp256(1:jlen2)
 1040   format(/" * Error - (EQLIB/swtchk) Can't replace the species ",
     $  a,/7x,'in the basis set with ',a,', because the former does',
     $  ' not appear',/7x,'in the linking reaction.')
c
c         Calling sequence substitutions:
c           noutpt for nf
c           ns2 for ns
c
          call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns2,nstmax,uspec)
          nerr = nerr + 1
        endif
      endif
c
      if (jflag(ns2) .ne. 30) then
c
c       Calling sequence substitutions:
c         nbaspx for nbasp
c         ns2 for ns
c
        nb2 = nbasis(nbaspx,nbt,nbtmax,ns2)
        if (nb2 .gt. 0) then
          write (noutpt,1050) usp156(1:jlen1),usp256(1:jlen2)
          write (nttyo,1050) usp156(1:jlen1),usp256(1:jlen2)
 1050     format(/" * Error - (EQLIB/swtchk) Can't replace the",
     $    ' species ',a,/7x,'in the basis set with ',a,', because',
     $    ' the latter is already in',/7x,'the active basis set.')
          nerr = nerr + 1
        endif
      endif
      if (nerr .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
