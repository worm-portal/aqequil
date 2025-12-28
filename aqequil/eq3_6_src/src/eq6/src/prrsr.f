      subroutine prrsr(cbsr,jcode,nbaspd,nbtd,nbtmax,nbt1mx,
     $ nf,noutpt,nrc,nrctmx,nrndex,nsrtmx,nstmax,nttyo,
     $ ureac,uspec)
c
c     This subroutine writes the reaction for the nsr-th special
c     reactant on the file whose unit number is nf.
c
c     This subroutine is similar in function to EQLIB/prreac.f, which
c     writes an ordinary reaction.
c
c     This subroutine is called by:
c
c       EQ6/chzrsr.f
c       EQ6/ckfrsr.f
c       EQ6/echoz.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       cbsr   = array of coefficients for reactions for special
c                  reactants
c       nbtmax = the maximum number of basis species
c       nbt1mx = the maximum number of basis species plus 1
c       nf     = the unit number of the file on which to write the
c                  reaction
c       nrc    = the reactant index of the special reactant whose
c                  reaction is to be printed
c       nsrtmx = the maximum number of special reactants
c       ureac  = array of names of reactants
c       uspec  = array of names of species
c
c     Principal output:
c
c       None
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nf,noutpt,nttyo
c
      integer nbtmax,nbt1mx,nrctmx,nsrtmx,nstmax
c
      integer jcode(nrctmx),nbaspd(nbtmax),nrndex(nrctmx)
c
      integer nbtd,nrc
c
      character*48 uspec(nstmax)
      character*24 ureac(nrctmx)
c
      real*8 cbsr(nbt1mx,nsrtmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen,j2,nb,ns,nsr
c
      integer ilnobl
c
      logical qfirst
c
      character*56 uspn56
c
      real*8 cx
c
c-----------------------------------------------------------------------
c
c     Make sure that the reactant is a special reactant.
c
      if (jcode(nrc) .ne. 2) then
        j2 = ilnobl(ureac(nrc))
        write (noutpt,1000) ureac(nrc)(1:j2)
        write (nttyo,1000) ureac(nrc)(1:j2)
 1000   format(/' * Error - (EQ6/prrsr.f) Programming error trap: this',
     $  ' subroutine was called',/7x,'to print the reaction for the',
     $  ' special reactant ',a,'. However,',/7x,'this reactant is not',
     $  ' a special reactant.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Get the special reactant index.
c
      nsr = nrndex(nrc)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (nf,1020)
 1020 format(1x)
c
c     The reactant.
c
      cx = -cbsr(nbt1mx,nsr)
      j2 = ilnobl(ureac(nrc))
      write (nf,1030) cx,ureac(nrc)(1:j2)
 1030 format(6x,1pe22.15,2x,a)
c
c     Reactants exclusive of the special reactant.
c
      do nb = 1,nbtd
        cx = cbsr(nb,nsr)
        if (cx .lt. 0.) then
          cx = -cx
          ns = nbaspd(nb)
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnx(jlen,uspec(ns),uspn56)
          write (nf,1040) cx,uspn56(1:jlen)
 1040     format(4x,'+ ',1pe22.15,2x,a)
        endif
      enddo
c
      write (nf,1050)
 1050 format(18x,'==')
c
c     Products.
c
      qfirst = .true.
      do nb = 1,nbtd
        cx = cbsr(nb,nsr)
        if (cx .gt. 0.) then
          ns = nbaspd(nb)
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnx(jlen,uspec(ns),uspn56)
          if (qfirst) then
            write (nf,1030) cx,uspn56(1:jlen)
            qfirst = .false.
          else
            write (nf,1040) cx,uspn56(1:jlen)
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
