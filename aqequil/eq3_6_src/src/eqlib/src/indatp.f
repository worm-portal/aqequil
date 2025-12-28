      subroutine indatp(apxa,axlksa,bpxa,cdrsa,cessa,iapxa_asv,iapxta,
     $ ibpxa_asv,ibpxta,ikta_asv,jsola,mwtspa,nad1,narx_asv,ncmpra,
     $ ndrsa,ndrsa_asv,ndrsn,ndrsra,nerr,nessa,nessa_asv,nessn,
     $ nessra,nmrn1a,nmrn2a,noutpt,np,npta_asv,ns,nsta_asv,ntpr_asv,
     $ nttyo,nxta,nxta_asv,qclnsa,uendit,uspeca,uphasa,uptsld,
     $ uptypa,zchara)
c
c     This subroutine reads the solid solution superblock from the
c     supporting data file "data1".
c
c     This subroutine is called by:
c
c       EQLIB/indata.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nad1   = unit number of the data file
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
      integer iapxa_asv,ibpxa_asv,ikta_asv,narx_asv,ndrsa_asv,nessa_asv,
     $ npta_asv,nsta_asv,ntpr_asv,nxta_asv
c
      integer noutpt,nttyo
c
      integer iapxta(nxta_asv),ibpxta(nxta_asv),jsola(nxta_asv),
     $ ncmpra(2,npta_asv),ndrsa(ndrsa_asv),ndrsra(2,nsta_asv),
     $ nessa(nessa_asv),nessra(2,nsta_asv)
c
      integer nad1,ndrsn,nerr,nessn,nmrn1a,nmrn2a,np,ns,nxta
c
      logical qclnsa(nsta_asv)
c
      character*24 uphasa(npta_asv),uptypa(npta_asv)
      character*24 uptsld
      character*48 uspeca(nsta_asv)
      character*8 uendit
c
      real*8 apxa(iapxa_asv,nxta_asv),
     $ axlksa(narx_asv,ntpr_asv,nsta_asv),
     $ bpxa(ibpxa_asv,nxta_asv),cdrsa(ndrsa_asv),cessa(nessa_asv),
     $ mwtspa(nsta_asv),zchara(nsta_asv)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3,jsolv,n,ncmptv,nsc,nx
c
      integer ilnobl
c
      character*24 uphasv
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with global dimensioning.
c     The following need not be SAVEd.
c
      character(len=24), dimension(:), allocatable :: ucompv
c
c-----------------------------------------------------------------------
c
c     Allocate  a scratch array for the names of the end-member
c     components.
c
      ALLOCATE(ucompv(ikta_asv))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The label below defines a loop for reading the blocks in the
c     superblock.
c
  100 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read a species block.
c
c     uphasv  = name of solid solution
c     ncmptv  = number of component species
c     jsolv   = code for solid solution activity coefficients
c
      read (nad1) uphasv,ncmptv,jsolv
      j3 = ilnobl(uphasv)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for end of this current superblock.
c
      if (uphasv(1:8) .eq. uendit(1:8)) go to 200
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      np = np + 1
      uphasa(np) = uphasv
      uptypa(np) = uptsld
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      nxta = nxta + 1
      nx = nxta
      jsola(nx) = jsolv
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ncmptv .gt. 0) then
c
c       Read the names of component species.
c
        read (nad1) (ucompv(n), n = 1,ncmptv)
      else
        write (noutpt,1020) uphasv(1:j3)
        write (noutpt,1020) uphasv(1:j3)
 1020   format(/' * Error - (EQLIB/indatp) The solid solution ',a,
     $  /7x,'has no components specified on the data file.')
       nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Read the activity coefficient parameters, first the ordinary
c     parameters (apx), then the site-mixing parameters (bpx).
c
      read (nad1) iapxta(nx)
      if (iapxta(nx) .gt. 0) read (nad1) (apxa(n,nx), n = 1,iapxta(nx))
      read (nad1) ibpxta(nx)
      if (ibpxta(nx) .gt. 0) read (nad1) (bpxa(n,nx), n = 1,ibpxta(nx))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Clone the component species.
c
      ncmpra(1,np) = ns + 1
      ncmpra(2,np) = ns + ncmptv
      do n = 1,ncmptv
        ns = ns + 1
        do nsc = nmrn1a,nmrn2a
          if (ucompv(n)(1:24) .eq. uspeca(nsc)(1:24)) then
            call clones(axlksa,cdrsa,cessa,mwtspa,narx_asv,
     $      ndrsa,ndrsa_asv,ndrsn,ndrsra,nessa,nessa_asv,nessn,
     $      nessra,np,npta_asv,ns,nsc,nsta_asv,ntpr_asv,uphasa,
     $      uspeca,zchara)
            qclnsa(ns) = .true.
            go to 130
          endif
        enddo
c
c       The listed component was not found as a pure mineral species.
c
        j2 = ilnobl(ucompv(n))
        write (noutpt,1040) ucompv(n)(1:j2),uphasv(1:j3)
        write (nttyo,1040) ucompv(n)(1:j2),uphasv(1:j3)
 1040   format(/" * Error - (EQLIB/indatp) Couldn't find the pure",
     $  ' phase equivalent',/7x,'of the species',/7x,a,' (',a,'),',
     $  ' therefore',/7x,'could not create it by cloning.')
        nerr = nerr + 1
  130   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Loop back to read another species block.
c
      go to 100
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Deallocate  the scratch array.
c
  200 DEALLOCATE(ucompv)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
