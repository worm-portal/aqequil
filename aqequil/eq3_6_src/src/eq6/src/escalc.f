      subroutine escalc(csts,iindx1,jcsort,kbt,kmax,moph,mosp,mtb,
     $ mtb0,nbaspd,nbt,nbtmax,ncmpr,noutpt,npt,nptmax,nstmax,nsts,
     $ nstsmx,nstsr,qprflg,uspec)
c
c     This subroutine recomputes the mass balance totals for the
c     Equilibrium System (ES). Normally, this is done immediately
c     after a shift of material from the ES to the Physically Removed
c     System (PRS). The calculation is performed by summing over what
c     remains in the ES, so it is not sensitive to the nature of
c     the material that was transferred.
c
c     This subroutine is called by:
c
c       EQ6/dumpdp.f
c       EQ6/pshfta.f
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
      integer kmax,nbtmax,nptmax,nstmax,nstsmx
c
      integer noutpt
c
      integer iindx1(kmax),jcsort(nstmax),nbaspd(nbtmax),
     $ ncmpr(2,nptmax),nsts(nstsmx),nstsr(2,nstmax)
c
      integer kbt,nbt,npt
c
      logical qprflg
c
      character*48 uspec(nstmax)
c
      real*8 csts(nstsmx),moph(nptmax),mosp(nstmax),mtb(nbtmax),
     $ mtb0(nbtmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer kcol,n,nb,np,nrn1,nrn2,nr1,nr2,ns,nss
c
c-----------------------------------------------------------------------
c
      call copyaa(mtb,mtb0,nbt)
c
      do nb = 1,nbt
        mtb(nb) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do np = 1,npt
        if (moph(np) .gt. 0.) then
          nrn1 = ncmpr(1,np)
          nrn2 = ncmpr(2,np)
          do nss = nrn1,nrn2
            ns = jcsort(nss)
            if (mosp(ns) .gt. 0.) then
              nr1 = nstsr(1,ns)
              nr2 = nstsr(2,ns)
              do n = nr1,nr2
                nb = nsts(n)
                mtb(nb) = mtb(nb) + csts(n)*mosp(ns)
              enddo
            endif
          enddo
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qprflg) then
        write (noutpt,1000)
 1000   format(/'   --- Mass Balance Totals ---',
     $  //3x,'Species',20x,'Old',10x,'New',/)
        do kcol = 1,kbt
          nb = iindx1(kcol)
          ns = nbaspd(nb)
          write (noutpt,1010) uspec(ns),mtb0(nb),mtb(nb)
 1010     format(1x,a24,2(2x,1pe12.5))
        enddo
        write (noutpt,1020)
 1020   format(1x)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      call copyaa(mtb,mtb0,nbt)
c
      end
