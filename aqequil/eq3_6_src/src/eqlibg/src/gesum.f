      subroutine gesum(conc,delam,elam,elsump,elsums,elsumw,
     $ fxi,narn1,narn2,nazpmx,nstmax,zchar)
c
c     This subroutine calculates three second order sums involving the
c     E-lambda function (elam) and its ionic strength derivatives
c     (delam). These are used in Pitzer's equations. These
c     sums are:
c
c       elsumw: SUM(ij) [E-lambda(ij) + I*E-lambda'(ij)]*m(i)*m(j)
c
c       elsums: SUM(ij) E-lambda'(ij)*m(i)*m(j)
c
c       elsump: SUM(ij) E-lambda''(ij)*m(i)*m(j) [Not used]
c
c     Here ij refers to any species pair, but non-zero
c     contributions only arise from cc' and aa' pairs. Note that:
c
c       elsumw = 2*
c             [SUM(c'>c) {E-theta(cc') + I*E-theta'(cc')}m(c)m(c')
c            + SUM(a'>a) {E-theta(aa') + I*E-theta'(aa')}m(a)m(a')]
c
c       elsums = 2*[SUM(c'>c) E-theta'(cc')m(c)m(c')
c                 + SUM(a'>a) E-theta'(aa')m(a)m(a')]
c
c       elsump = 2*[SUM(c'>c) E-theta''(cc')m(c)m(c')
c                 + SUM(a'>a) E-theta''(aa')m(a)m(a')] [Not used]
c
c     This subroutine is called by:
c
c       EQLIBG/gcoeff.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       conc   = array of concentrations
c       elam   = array of E-lambda functions
c       delam  = array of ionic strength derivatives of E-lambda
c                  functions
c       fxi    = the ionic strength (the 2nd-order electrostatic
c                  moment function I)
c       narn1  = start of species range for aqueous solution
c       narn2  = end of species range for aqueous solution
c       zchar  = array of charges
c
c     Principal output:
c
c       elsumw = sum: SUM(ij) [ E-lambda(ij) + I*E-lambda'(ij) ]*mi*mj
c       elsums = sum: SUM(ij) E-lambda'(ij)*mi*mj
c       elsump = sum: SUM(ij) E-lambda''(ij)*mi*mj
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nazpmx,nstmax
c
      integer narn1,narn2
c
      real*8 conc(nstmax),delam(2,nazpmx,nazpmx),elam(nazpmx,nazpmx),
     $ zchar(nstmax)
      real*8 elsump,elsums,elsumw,fxi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ijz,iz,jz,nsi,nsj
c
      real*8 cij,el,elp,elpp,zi,zj
c
c-----------------------------------------------------------------------
c
      elsumw = 0.
      elsums = 0.
      elsump = 0.
c
      do nsi = narn1,narn2
        zi = zchar(nsi)
        iz = nint(abs(zi))
        if (iz .ne. 0) then
          do nsj = narn1,narn2
            if (nsi .ne. nsj) then
              zj = zchar(nsj)
              jz = nint(abs(zj))
              if (jz .ne. 0) then
                ijz = nint(zchar(nsi)*zchar(nsj))
                if ((ijz) .gt. 0) then
                  el = elam(iz,jz)
                  elp = delam(1,iz,jz)
                  elpp = delam(2,iz,jz)
                  cij = conc(nsi)*conc(nsj)
                  elsumw = elsumw + (el + fxi*elp)*cij
                  elsums = elsums + elp*cij
                  elsump = elsump + elpp*cij
                endif
              endif
            endif
          enddo
        endif
      enddo
c
      end
