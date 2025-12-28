      subroutine clones(axlksa,cdrsa,cessa,mwtspa,narx_asv,
     $ ndrsa,ndrsa_asv,ndrsn,ndrsra,nessa,nessa_asv,nessn,
     $ nessra,np,npta_asv,ns,nsc,nsta_asv,ntpr_asv,uphasa,
     $ uspeca,zchara)
c
c     This subroutine clones the nsc-th species into the ns-th. The
c     ns-th species belongs to the np-th phase. Typically, the nsc-th
c     species is a pure mineral or liquid and the ns-th species is the
c     corresponding component species of a solid or liquid solution.
c
c     This subroutine is called by:
c
c       EQLIB/indata.f
c       EQLIB/indatp.f
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
      integer narx_asv,ndrsa_asv,nessa_asv,npta_asv,nsta_asv,ntpr_asv
c
      integer ndrsa(ndrsa_asv),ndrsra(2,nsta_asv),nessa(nessa_asv),
     $ nessra(2,nsta_asv)
c
      integer ndrsn,nessn,np,ns,nsc
c
      character*24 uphasa(npta_asv)
      character*48 uspeca(nsta_asv)
c
      real*8 axlksa(narx_asv,ntpr_asv,nsta_asv),cessa(nessa_asv),
     $ cdrsa(ndrsa_asv),mwtspa(nsta_asv),zchara(nsta_asv)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j,ndrsn1,nessn1,nr1,nr2,nse,nt
c
c----------------------------------------------------------------------
c
c     Copy the name, molecular weight, and electrical charge.
c     Note that the phase part of the name is different.
c
      uspeca(ns)(1:24) = uspeca(nsc)(1:24)
      uspeca(ns)(25:48) = uphasa(np)(1:24)
      mwtspa(ns) = mwtspa(nsc)
      zchara(ns) = zchara(nsc)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Copy the elemental composition.
c
      nessra(1,ns) = nessn + 1
      nr1 = nessra(1,nsc)
      nr2 = nessra(2,nsc)
      do nessn1 = nr1,nr2
        nessn = nessn + 1
        cessa(nessn) = cessa(nessn1)
        nessa(nessn) = nessa(nessn1)
      enddo
      nessra(2,ns) = nessn
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Copy the associated reaction.
c
      ndrsra(1,ns) = ndrsn + 1
      nr1 = ndrsra(1,nsc)
      nr2 = ndrsra(2,nsc)
      nt = nr2 - nr1 + 1
      if (nt .lt. 2) then
c
c       The nsc-th species is a data file basis species that
c       has only an identity reaction. Write a linking reaction.
c
        ndrsn = ndrsn + 1
        cdrsa(ndrsn) = -1.
        ndrsa(ndrsn) = ns
        ndrsn = ndrsn + 1
        cdrsa(ndrsn) = 1.
        ndrsa(ndrsn) = nsc
      else
c
c       The nsc-th species has a reaction which can be copied.
c
        do ndrsn1 = nr1,nr2
          ndrsn = ndrsn + 1
          cdrsa(ndrsn) = cdrsa(ndrsn1)
          nse = ndrsa(ndrsn1)
          if (nse .eq. nsc) nse = ns
          ndrsa(ndrsn) = nse
        enddo
      endif
      ndrsra(2,ns) = ndrsn
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Copy the coefficients of the polynomials which represent
c     the thermodynamic functions of the reaction. Note that
c     this is valid whether the reaction is an identity reaction
c     or a real one.
c
      do j = 1,ntpr_asv
        do i = 1,narx_asv
          axlksa(i,j,ns) = axlksa(i,j,nsc)
        enddo
      enddo
c
      end
