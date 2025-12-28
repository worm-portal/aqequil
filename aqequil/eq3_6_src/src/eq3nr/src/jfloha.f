      subroutine jfloha(jflgi,nbtd,nbti,nbtmax,noutpt,nttyo,
     $ nxmdmx,nxmod,ubasp,uspeci,uxmod)
c
c     This subroutine modifies the jflgi array to insure that jflag
c     defaults to 27, not 30, for the species H2(aq) and O2(aq).
c     These defaults are desirable because they provide some
c     poising (buffering of the fO2 or Eh) in EQ6 mode calculations.
c     However, the jflag is defaulted to -1 if a species is on the
c     suppression list (nxmod option).
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
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
      integer nbtmax,nxmdmx
c
      integer noutpt,nttyo
c
      integer jflgi(nbtmax)
      integer nbtd,nbti,nxmod
c
      character(len=48) ubasp(nbtmax),uspeci(nbtmax),uxmod(nxmdmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,nb,nbi
c
      character(len=24) ux24
c
c-----------------------------------------------------------------------
c
c     Process H2(aq).
c
      ux24 = 'H2(aq)'
c
c     Is it on the input file?
c
      do nbi = 1,nbti
        if (uspeci(nbi)(1:24) .eq. ux24(1:24)) go to 110
      enddo
c
c     Is it on the data file?
c
      do nb = 1,nbtd
        if (ubasp(nb)(1:24) .eq. ux24(1:24)) go to 100
      enddo
      go to 110
c
c     See if the input can be extended.
c
  100 if (nbti .ge. nbtmax) then
        write (noutpt,1000) nbtmax,ux24(1:6)
        write (nttyo,1000) nbtmax,ux24(1:6)
 1000   format(/' * Error - (EQ3NR/jfloha) Have exceeded the maximum',
     $  /7x,i4,' basis species trying to add ',a,
     $  /7x,'to the species read from the input file.',
     $  /7x,'Increase the dimensioning parameter nbtpar.')
        stop
      endif
c
c     Extend the input with the desired default.
c
      nbti = nbti + 1
      uspeci(nbti) = ux24
c
c     Is it on the suppression list?
c
      do n = 1,nxmod
        if (uxmod(n)(1:24) .eq. ux24(1:24)) then
          jflgi(nbti) = -1
          go to 110
        endif
      enddo
c
      jflgi(nbti) = 27
c
  110 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Process O2(aq).
c
      ux24 = 'O2(aq)'
c
c     Is it on the input file?
c
      do nbi = 1,nbti
        if (uspeci(nbi)(1:24) .eq. ux24(1:24)) go to 210
      enddo
c
c     Is it on the data file?
c
      do nb = 1,nbtd
        if (ubasp(nb)(1:24) .eq. ux24(1:24)) go to 200
      enddo
      go to 210
c
c     See if the input can be extended.
c
  200 if (nbti .ge. nbtmax) then
        write (noutpt,1000) nbtmax,ux24(1:6)
        write (nttyo,1000) nbtmax,ux24(1:6)
        stop
      endif
c
c     Extend the input with the desired default.
c
      nbti = nbti + 1
      uspeci(nbti) = ux24
c
c     Is it on the suppression list?
c
      do n = 1,nxmod
        if (uxmod(n)(1:24) .eq. ux24(1:24)) then
          jflgi(nbti) = -1
          go to 210
        endif
      enddo
c
      jflgi(nbti) = 27
c
  210 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
