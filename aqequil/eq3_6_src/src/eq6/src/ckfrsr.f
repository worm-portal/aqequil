      subroutine ckfrsr(cbsr,csts,jcode,jflag,nbaspd,nbtd,nbtmax,
     $ nbt1mx,noutpt,nrct,nrctmx,nrndex,nsrtmx,nstmax,nsts,nstsmx,
     $ nstsr,nttyo,ureac,uspec)
c
c     This subroutine checks the reactions for special reactants. As
c     needed, it rewrites these reactions to eliminate any basis
c     species for which jflag = 30. Such species are not in the active
c     basis set, and there are no corresonding mass balances.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       cbsr   = array of coefficients for reactions for special
c                  reactants
c       nbtmax = the maximum number of basis species
c       nbt1mx = the maximum number of basis species plus 1
c       nsrtmx = the maximum number of special reactants
c       ureac  = array of names of reactants
c       uspec  = array of names of species
c
c     Principal output:
c
c       cbsr   = array of coefficients for reactions for special
c                  reactants (modified if necessary)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbtmax,nbt1mx,nrctmx,nsrtmx,nstmax,nstsmx
c
      integer noutpt,nttyo
c
      integer jcode(nrctmx),jflag(nstmax),nbaspd(nbtmax),nrndex(nrctmx),
     $ nsts(nstsmx),nstsr(2,nstmax)
c
      integer nbtd,nrct
c
      character*48 uspec(nstmax)
      character*24 ureac(nrctmx)
c
      real*8 cbsr(nbt1mx,nsrtmx),csts(nstsmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nf
c
      integer jlen,j2,n,nb,nbb,nrc,nr1,nr2,ns,nsr
c
      integer ilnobl
c
      logical qcaught
c
      character*56 uspn56
c
      real*8 cx
c
c-----------------------------------------------------------------------
c
      do nrc = 1,nrct
        if (jcode(nrc) .eq. 2) then
          nsr = nrndex(nrc)
          qcaught = .false.
          do nb = 1,nbtd
            ns = nbaspd(nb)
            cx = cbsr(nb,nsr)
            if (cx.ne.0. .and. jflag(ns).eq.30) then
              if (.not.qcaught) then
c
c               Write the original reaction.
c
                j2 = ilnobl(ureac(nrc))
                write (noutpt,1000) ureac(nrc)(1:j2)
                write (nttyo,1000) ureac(nrc)(1:j2)
 1000           format(/1x,'The reaction for special reactant ',a,
     $          ' is written in terms of',/1x,'one or more species',
     $          ' which are not in the active basis set. The',/1x,
     $          'currently written reaction is:')
c
                nf = noutpt
                call prrsr(cbsr,jcode,nbaspd,nbtd,nbtmax,nbt1mx,
     $          nf,noutpt,nrc,nrctmx,nrndex,nsrtmx,nstmax,nttyo,
     $          ureac,uspec)
                write (noutpt,1010)
 1010           format(/1x,'The following species are not in the',
     $          ' active basis set and will be eliminated',
     $          /1x,'from the reaction:',/)
c
                nf = nttyo
                call prrsr(cbsr,jcode,nbaspd,nbtd,nbtmax,nbt1mx,
     $          nf,noutpt,nrc,nrctmx,nrndex,nsrtmx,nstmax,nttyo,
     $          ureac,uspec)
                write (nttyo,1010)
c
                qcaught = .true.
              endif
c
c             Calling sequence substitutions:
c               uspec(ns) for unam48
c
              call fmspnx(jlen,uspec(ns),uspn56)
c
              write (noutpt,1020) uspn56(1:jlen)
              write (nttyo,1020) uspn56(1:jlen)
 1020         format(3x,a)
c
c             The reaction contains a non-zero reaction coefficient for
c             a species which is not in the active basis set. Rewrite
c             the reaction so that this species does not appear.
c
              nr1 = nstsr(1,ns)
              nr2 = nstsr(2,ns)
              do n = nr1,nr2
                nbb = nsts(n)
                cbsr(nbb,nsr) = cbsr(nbb,nsr) + csts(n)*cx
              enddo
              cbsr(nb,nsr) = 0.
            endif
          enddo
c
          if (qcaught) then
c
c           Write the modified reaction.
c
            write (noutpt,1070)
            write (nttyo,1070)
 1070       format(//1x,'The modified reaction is:')
c
            nf = noutpt
            call prrsr(cbsr,jcode,nbaspd,nbtd,nbtmax,nbt1mx,
     $      nf,noutpt,nrc,nrctmx,nrndex,nsrtmx,nstmax,nttyo,
     $      ureac,uspec)
            write (noutpt,1030)
 1030       format(/1x)
c
            nf = nttyo
            call prrsr(cbsr,jcode,nbaspd,nbtd,nbtmax,nbt1mx,
     $      nf,noutpt,nrc,nrctmx,nrndex,nsrtmx,nstmax,nttyo,
     $      ureac,uspec)
            write (nttyo,1030)
          endif
c
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
