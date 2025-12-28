      subroutine rxnchk(cdrs,cess,mtotr,nbt,nbtmx1,nbtmx2,nco,
     $ nct,nctmax,nerr,noutpt,ns,nsb,nttyo,uelem,uspec,zchar)
c
c     This suboutine checks the reaction associated with the ns-th
c     species for mass and charge balance. If an imbalance is found,
c     a message is written to the screen and output files.
c
c     This suboutine is called by:
c
c       EQPT/pcraq.f
c       EQPT/pcrsg.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       cdrs   = array of reaction coefficients
c       cess   = array of elemental composition coefficients
c       nbt    = the number of basis species
c       nct    = the number of chemical elements
c       nerr   = cumulative error counter
c       uelem  = array of chemical element names
c       uspec  = array of species names
c       zchar  = array of electrical charge numbers
c
c     Principal output:
c
c       nerr   = cumulative error counter
c
c     Workspace:
c
c       mtotr  = array of mass balance residuals
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbtmx1,nbtmx2
c
      integer noutpt,nttyo
c
      integer nbt,nco,nct,nctmax,nerr,ns,nsb
c
      character*24 uspec(nbtmx1)
      character*8 uelem(nctmax)
      character*8 ux8
c
      real*8 cdrs(nbtmx2,nbtmx1),cess(nctmax,nbtmx1),mtotr(nctmax),
     $ zchar(nbtmx1)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      logical qbadr,qo2gx
c
      integer j2,j3,nbt1,nc,nfile,nlim,nse
c
      integer ilnobl
c
      real*8 mx,rxntol,ztotr
c
c-----------------------------------------------------------------------
c
      data rxntol /1.e-6/
c
c-----------------------------------------------------------------------
c
      qbadr = .false.
c
      qo2gx = uspec(nsb) .eq. 'O2(g)'
      if (qo2gx) cess(nco,nsb) = 2.
c
      nlim = ns - 1
      nlim = min(nlim,nbt)
      nbt1 = nbt + 1
c
      ztotr = cdrs(nbt1,ns)*zchar(ns)
      do nse = 1,nlim
        ztotr = ztotr + cdrs(nse,ns)*zchar(nse)
      enddo
      if (abs(ztotr) .gt. rxntol) qbadr = .true.
c
      do nc = 1,nct
        mx = cdrs(nbt1,ns)*cess(nc,ns)
        do nse = 1,nlim
          mx = mx + cdrs(nse,ns)*cess(nc,nse)
        enddo
        mtotr(nc) = mx
        if (abs(mx) .gt. rxntol) qbadr = .true.
      enddo
c
      if (qo2gx) cess(nco,nsb) = 0.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qbadr) then
c
c       Have a reaction with one or more imbalances.
c
        j2 = ilnobl(uspec(ns))
        write (noutpt,1000) uspec(ns)(1:j2)
        write (nttyo,1000) uspec(ns)(1:j2)
 1000   format(/' * Error - (EQPT/rxnchk) The reaction for the',
     $  ' destruction of',/7x,a,' has the following imbalances:',/)
        if (abs(ztotr) .gt. rxntol) then
          ux8 = 'charge'
          j3 = ilnobl(ux8)
          write (noutpt,1010) ztotr,ux8(1:j3)
          write (nttyo,1010) ztotr,ux8(1:j3)
 1010     format(9x,1pe12.5,' in ',a)
        endif
        do nc = 1,nct
          if (abs(mtotr(nc)) .gt. rxntol) then
            ux8 = uelem(nc)
            j3 = ilnobl(ux8)
            write (noutpt,1010) mtotr(nc),ux8(1:j3)
            write (nttyo,1010) mtotr(nc),ux8(1:j3)
          endif
        enddo
        write (noutpt,1020)
        write (nttyo,1020)
 1020   format(/9x,'The reaction is:')
        nfile = noutpt
        call prrecy(cdrs,nbtmx1,nbtmx2,nbt,ns,nfile,uspec)
        nfile = nttyo
        call prrecy(cdrs,nbtmx1,nbtmx2,nbt,ns,nfile,uspec)
        nerr = nerr + 1
      endif
c
      end
