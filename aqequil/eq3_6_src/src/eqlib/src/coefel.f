      real*8 function coefel(cess,ness,nessmx,nessr,nc,ns,nstmax)
c
c     This subroutine finds the coefficient of the nc-th element
c     in the composition of the ns-th species.
c
c     This subroutine is called by:
c
c       None
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
      integer nessmx,nstmax
c
      integer ness(nessmx),nessr(2,nstmax)
      integer nc,ns
c
      real*8 cess(nessmx)
c
c----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer n,n1,n2
c
c-----------------------------------------------------------------------
c
      coefel = 0.
      n1 = nessr(1,ns)
      n2 = nessr(2,ns)
      do n = n1,n2
        if (nc .eq. ness(n)) then
          coefel = cess(n)
          go to 999
        endif
      enddo
c
  999 continue
      end
