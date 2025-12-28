      subroutine wterm(apx,iapxmx,iktmax,ixrn1,ixrn2,jsol,ncmpr,
     $ noutpt,nptmax,nstmax,nttyo,nxt,nxtmax,press,tempk,uphase,
     $ uspec,wfac)
c
c     This subroutine computes the wfac(i,nx) array, which contains
c     the coefficients for the excess free energy function of
c     solid solutions. If non-zero coefficients are lacking, the
c     solid solution reduces to an ideal solution and the jsol value
c     is changed to zero to reflect this.
c
c     This subroutine is called by:
c
c       EQLIB/evdata.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       apx    = array of non-ideal mixing parameters
c       jsol   = array of non-ideal mixing law flags
c       press  = pressure, bars
c       tempk  = temperature, K
c       uphase = array of phase names
c
c     Principal input:
c
c       wfac   = array of non-ideal mixing parameters calculated from
c                  the apx array
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer iapxmx,iktmax,nptmax,nstmax,nxtmax
c
      integer jsol(nxtmax),ncmpr(2,nptmax)
      integer ixrn1,ixrn2,noutpt,nttyo,nxt
c
      character*48 uspec(nstmax)
      character*24 uphase(nptmax)
c
      real*8 apx(iapxmx,nxtmax),wfac(iktmax,nxtmax)
      real*8 press,tempk
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ix,iy,j2,k,np,nx
c
      integer ilnobl
c
      character*48 ux48,uy48
c
      real*8 sum
c
c-----------------------------------------------------------------------
c
c     Note: the following statements don't really do anything except
c     cause the compiler not to complain that uspec, ncmpr, and ixrn2
c     are not used.
c
      ux48 = uspec(1)
      uy48 = ux48
      uspec(1) = uy48
c
      ix = ncmpr(1,1)
      iy = ix
      ncmpr(1,1) = iy
c
      ix = ixrn2
      iy = ix
      ixrn2 = iy
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      np = ixrn1 - 1
      do nx = 1,nxt
        np = np + 1
        sum = 0.
        do i = 1,iapxmx
          sum = sum + apx(i,nx)*apx(i,nx)
        enddo
c
        if (jsol(nx).eq.1 .and. sum.gt. 0.) then
          j2 = ilnobl(uphase(np))
          write (noutpt,1735) uphase(np)(1:j2)
          write (nttyo,1735) uphase(np)(1:j2)
1735      format(/' * Note - (EQLIBG/wterm) The phase ',a,' is',
     $    /7x,'listed as ideal. This is inconsistent with the',
     $    /7x,'presence of non-zero data in the apx array. These',
     $    /7x,'data will be ignored.')
          do i = 1,iapxmx
            apx(i,nx) = 0.
          enddo
        endif
c
        if (jsol(nx).ne.1 .and. sum.eq.0.) then
          j2 = ilnobl(uphase(np))
          write (noutpt,1736) uphase(np)(1:j2)
          write (nttyo,1736) uphase(np)(1:j2)
1736      format(/' * Note - (EQLIBG/wterm) The phase ',a,' is',
     $    /7x,'listed as non-ideal, but no non-zero data are present'
     $    /7x,'in the apx array. This phase will be treated as ideal.')
          jsol(nx) = 1
        endif
c
c       Initialize wfac values to zero
c
        do i = 1,iktmax
          wfac(i,nx) = 0.
        enddo
c
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      np = ixrn1 - 1
      do nx = 1,nxt
        np = np + 1
        k = jsol(nx)
        if (k .eq. 1) then
c
c         Ideal solution.
c
          continue
        elseif (k.eq.2) then
c
c         binary solution, third-order maclaurin expansion
c         original pathi solid solution model.
c
          wfac(1,nx) = apx(1,nx)
          wfac(2,nx) = apx(2,nx)
          wfac(3,nx) = apx(3,nx)
          wfac(1,nx) = -wfac(2,nx)/2. - wfac(3,nx)/6.
c
        elseif (k.eq.3) then
c
c       binary solution, parabolic maclaurin expansion
c
          wfac(1,nx) = apx(1,nx)
c
        elseif (k.eq.4) then
c
c         binary solution, cubic maclaurin (p,t dependent)
c
          wfac(1,nx) = apx(1,nx) + apx(2,nx)*tempk + apx(3,nx)*press
          wfac(2,nx) = apx(4,nx) + apx(5,nx)*tempk + apx(6,nx)*press
c
        elseif (k.eq.5) then
c
c         binary solution, guggenheim polynomial (t dependent)
c
          wfac(1,nx) = apx(1,nx) + apx(2,nx)*tempk
     $    + apx(3,nx)*tempk*tempk
          wfac(2,nx) = apx(4,nx) + apx(5,nx)*tempk
     $    + apx(6,nx)*tempk*tempk
          wfac(3,nx) = apx(7,nx) + apx(8,nx)*tempk
     $    + apx(9,nx)*tempk*tempk
c
        elseif (k.eq.6) then
c
c         ternary regular solution (see prigogine and defay, p. 257)
c
          wfac(1,nx) = apx(1,nx)
          wfac(2,nx) = apx(2,nx)
          wfac(3,nx) = apx(3,nx)
c
        elseif (k.eq.7) then
c
c         newton et al plagioclase model (gca vol 44 p. 933, 1980)
c         1 - albite; 2 - anorthite
c
          wfac(1,nx) = apx(1,nx)
          wfac(2,nx) = apx(2,nx)
        else
          j2 = ilnobl(uphase(np))
          write (noutpt,100) jsol(nx),uphase(np)(1:j2)
          write (nttyo,100) jsol(nx),uphase(np)(1:j2)
 100      format(/' * Error - (EQLIBG/wterm) Have an Illegal jsol',
     $    /7x,'value of ',i2,' for solid solution ',a,'.')
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
