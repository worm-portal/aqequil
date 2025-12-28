      subroutine lambda(acflgc,afcnst,bpx,ibpxmx,ibpxt,iktmax,ixrn1,
     $ ixrn2,jsol,ncmpr,noutpt,np,nptmax,nstmax,nttyo,nxtmax,wfac,
     $ xbar,xbarlg,uphase,uspec)
c
c     This subroutine computes activity coefficients for the components
c     of the np-th phase, which is a solid solution.
c
c     This subroutine is called by:
c
c       EQLIB/ngcadv.f
c       EQ3NR/eq3nr.f
c       EQ6/raff.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       afcnst = affinity constant, 2.303RT/1000
c       bpx    = array of site-mixing parameters for solid solution
c                  activity coefficients
c       ibpxt  = array of the numbers of non-zero site-mixing
c                  paramters for computing activity coefficients in
c                  solid solutions
c       jsol   = array of activity coefficient model flags
c       ncmpr  = species range pointer array for phases
c       noutpt = unit number of the output file
c       nttyo  = unity number of the tty output file
c       nx     = solid solution index
c       nxtmax = array dimension, maximum number of solid solutions
c       wfac   = array of non-ideal mixing parameters calculated from
c                   the apx array
c       xbar   = array of mole fractions
c       xbarlg = array of log mole fraction values
c       wterm  = array of temperature-dependent coefficients for
c                  solid solution activity coefficients
c
c     Principal output:
c
c       acflgc  = array of activity coefficients
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ibpxmx,iktmax,nptmax,nstmax,nxtmax
c
      integer noutpt,nttyo
c
      integer ibpxt(nxtmax),jsol(nxtmax),ncmpr(2,nptmax)
      integer ixrn1,ixrn2,np
c
      character*48 uspec(nstmax)
      character*24 uphase(nptmax)
c
      real*8 acflgc(nstmax),bpx(ibpxmx,nxtmax),xbar(nstmax),
     $ xbarlg(nstmax),wfac(iktmax,nxtmax)
      real*8 afcnst
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ix,iy,j2,k,nr1,nr2,ns,nt,nx
c
      integer ilnobl
c
      character*48 ux48,uy48
c
      real*8 stxm1,xx,xy
c
c-----------------------------------------------------------------------
c
c     Note: the following statements don't really do anything except
c     cause the compiler not to complain that ixrn2, uspec, afcnst,
c     xbar, and wfac are not used.
c
      ix = ixrn2
      iy = ix
      ixrn2 = iy
c
      ux48 = uspec(1)
      uy48 = ux48
      uspec(1) = uy48
c
      xx = afcnst
      xy = xx
      afcnst = xy
c
      xx = xbar(1)
      xy = xx
      xbar(1) = xy
c
      xx = wfac(1,1)
      xy = xx
      wfac(1,1) = xy
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      nx = np - ixrn1 + 1
cXX   Temporarily make all solid solutions ideal (for version 8.0).
cXX   Non-ideal solid solutions will be reintroduced in some later
cXX   version.
      if (jsol(nx) .ne. 1) jsol(nx) = 1
cXX
      k = jsol(nx)
c
      nr1 = ncmpr(1,np)
      nr2 = ncmpr(2,np)
      nt = nr2 - nr1 + 1
      if (nt .le. 0) then
        j2 = ilnobl(uphase(np))
        write (noutpt,1000) uphase(np)(1:j2)
        write (nttyo,1000) uphase(np)(1:j2)
 1000   format(/' * Error - (EQLIBG/lambda) Programming error trap:',
     $  ' Have no',/7x,'components present for solid solution ',a,'.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (k .eq. 1) then
c
        if (ibpxt(nx) .eq. 0.) then
c
c         Molecular mixing ideal solution.
c
          do ns = nr1,nr2
            acflgc(ns) = 0.
          enddo
        else
c
c         Site mixing ideal solution.
c
          stxm1 = bpx(1,nx) - 1.0
          do ns = nr1,nr2
            acflgc(ns) = stxm1*xbarlg(ns)
          enddo
        endif
c
      else
        j2 = ilnobl(uphase(np))
        write (noutpt,1010) jsol(nx),nx,uphase(np)(1:j2)
        write (nttyo,1010) jsol(nx),nx,uphase(np)(1:j2)
 1010   format(/' * Error - (EQLIBG/lambda) Programming error trap:',
     $  ' Have an',/7x,'undefined value of ',i3,' for jsol(',i3,'),',
     $  ' the activity'/7x,'coefficient model flag for solid',
     $  ' solution ',a,'.')
        stop
      endif
c
      end
