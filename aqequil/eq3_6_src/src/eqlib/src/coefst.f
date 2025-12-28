      real*8 function coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
c
c     This subroutine finds the stoichiometric coefficient of the nb-th
c     basis species for the ns-th species. Note that 'nb' denotes
c     position in the list of basis species, not the general list of
c     species. The corresponding index on the latter list would be
c     given by nbasp(nb).
c
c     This subroutine is called by:
c
c       EQLIB/prtpct.f
c       EQ3NR/betas.f
c       EQ3NR/matrix.f
c       EQ3NR/scripx.f
c       EQ6/betaz.f
c       EQ6/matrxz.f
c       EQ6/pabssw.f
c       EQ6/setffg.f
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
      integer nstsmx,nstmax
c
      integer nsts(nstsmx),nstsr(2,nstmax)
      integer nb,ns
c
      real*8 csts(nstsmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,n1,n2
c
c-----------------------------------------------------------------------
c
      coefst = 0.
      n1 = nstsr(1,ns)
      n2 = nstsr(2,ns)
      do n = n1,n2
        if (nb .eq. nsts(n)) then
          coefst = csts(n)
          go to 999
        endif
      enddo
c
  999 continue
      end
