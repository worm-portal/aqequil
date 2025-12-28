      subroutine chsgex(adhfs,adhfsd,adhfsx,advfs,advfsd,advfsx,
     $ axhfs,axhfsd,axhfsx,axlks,axlksd,axlksx,axvfs,axvfsd,axvfsx,
     $ cdrs,cdrsd,cdrsx,eps100,iern1,ipch,ipchmx,ipcv,ipcvmx,jern1,
     $ jetmax,jflag,jgext,jsflag,narn1,narxmx,narxt,nbasp,nbaspd,
     $ nbaspx,nbt,nbtmax,nbw,ndrs,ndrsd,ndrsmx,ndrsr,ndrsrd,ndrsrx,
     $ ndrsx,nern1,nern2,net,netmax,ngext,noutpt,nphasx,nst,nstmax,
     $ ntprmx,ntprt,nttyo,qbassw,qbswok,ugexmo,uspec)
c
c     This subroutine changes the original setup of component species
c     of generic ion exchanger phases for certain exchange models
c     (e.g., Gapon, Vanselow). It switches bare site species out of
c     the active basis set; e.g., it replaces __-Z with Na-Z. The
c     thermodynamic properties and the jsflag status array element
c     for each such base site species are changed to effectively take
c     this species out of the model. It persists only as a reference
c     for the mass balance of the corresponding exchanger substrate.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       narxmx = number of coefficient elements of axlks per species
c                  per temperature range
c       nbtmax = maximum number of basis species
c       nstmax = maximum number of species
c       ntprmx = maximum number of temperature ranges
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
      integer ipchmx,ipcvmx,jetmax,narxmx,nbtmax,ndrsmx,netmax,nstmax,
     $ ntprmx
c
      integer noutpt,nttyo
c
      integer jern1(jetmax,netmax),jflag(nstmax),jgext(netmax),
     $ jsflag(nstmax),narxt(ntprmx),nbasp(nbtmax),nbaspd(nbtmax),
     $ nbaspx(nbtmax),ndrs(ndrsmx),ndrsd(ndrsmx),ndrsr(2,nstmax),
     $ ndrsrd(2,nstmax),ndrsx(ndrsmx),ndrsrx(2,nstmax),
     $ ngext(jetmax,netmax),nphasx(nstmax)
c
      integer iern1,ipch,ipcv,narn1,nbt,nbw,nern1,nern2,net,nst,ntprt
c
      logical qbassw,qbswok
c
      character*24 ugexmo(netmax)
      character*48 uspec(nstmax)
c
      real*8 adhfs(narxmx,ntprmx,ipchmx,nstmax),
     $ adhfsd(narxmx,ntprmx,ipchmx,nstmax),
     $ adhfsx(narxmx,ntprmx,ipchmx,nstmax),
     $ advfs(narxmx,ntprmx,ipcvmx,nstmax),
     $ advfsd(narxmx,ntprmx,ipcvmx,nstmax),
     $ advfsx(narxmx,ntprmx,ipcvmx,nstmax),
     $ axhfs(narxmx,ntprmx,nstmax),axhfsd(narxmx,ntprmx,nstmax),
     $ axhfsx(narxmx,ntprmx,nstmax),
     $ axlks(narxmx,ntprmx,nstmax),axlksd(narxmx,ntprmx,nstmax),
     $ axlksx(narxmx,ntprmx,nstmax),axvfs(narxmx,ntprmx,nstmax),
     $ axvfsd(narxmx,ntprmx,nstmax),axvfsx(narxmx,ntprmx,nstmax),
     $ cdrs(ndrsmx),cdrsd(ndrsmx),cdrsx(ndrsmx)
c
      real*8 eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ie,j,je,j2,nb,ne,np,ns,ns1,ns2
c
      integer ilnobl
c
c-----------------------------------------------------------------------
c
c     Copy the existing nbasp set into the nbaspx array.
c
      do nb = 1,nbt
        nbaspx(nb) = nbasp(nb)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Switch bare site species of generic ion exchange phases with
c     certain exchange models (e.g., Gapon, Vanselow) out of the
c     basis set. E.g., use Na-Z instead of __-Z.
c
      do nb = 1,nbt
        ns1 = nbaspx(nb)
        if (ns1.ge.nern1 .and. ns1.le.nern2) then
          np = nphasx(ns1)
          ne = np - iern1 + 1
          j2 = ilnobl(ugexmo(ne))
          if (ugexmo(ne)(1:j2).eq.'Gapon' .or.
     $      ugexmo(ne)(1:6).eq.'Gapon-' .or.
     $      ugexmo(ne)(1:j2).eq.'Vanselow' .or.
     $      ugexmo(ne)(1:9).eq.'Vanselow-') then
c
c           Theoretically, the test below should be unnecessary,
c           as the only exchanger species presently in the basis
c           set should be bare site species.
c
            if (uspec(ns1)(1:3) .eq. '__ ') then
c
c             Set up to switch the bare exchanger species with the
c             next species in the same site.
c
              ns2 = ns1 + 1
              nbasp(nb) = ns2
            endif
c
          endif
        endif
      enddo
c
c     Make the switches.
c
      do nb = 1,nbt
        ns1 = nbaspx(nb)
        ns2 = nbasp(nb)
        if (ns1 .ne. ns2) then
          call switch(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,
     $    axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,
     $    ipcv,ipcvmx,jflag,jsflag,narn1,narxmx,nbasp,nbaspd,nbaspx,
     $    nb,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,noutpt,
     $    ns2,nst,nstmax,ntprmx,nttyo,qbassw,qbswok,uspec)
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Now for certain models (e.g., Gapon, Vanselow) adjust the
c     thermodynamic data for the bare site species so that these
c     species are effectively taken out of the model.
c
      do ne = 1,net
        j2 = ilnobl(ugexmo(ne))
        if (ugexmo(ne)(1:j2).eq.'Gapon' .or.
     $    ugexmo(ne)(1:6).eq.'Gapon-' .or.
     $    ugexmo(ne)(1:j2).eq.'Vanselow' .or.
     $    ugexmo(ne)(1:9).eq.'Vanselow-') then
c
          do je = 1,jgext(ne)
            ns = jern1(je,ne) - 1
            do ie = 1,ngext(je,ne)
              ns = ns + 1
c
              if (uspec(ns)(1:3) .eq. '__ ') then
c
c               Null the thermodynamic properties of the
c               bare site species.
c
                do j = 1,ntprt
                  axlks(1,j,ns) = +9999999.
                  axhfs(1,j,ns) = 0.
                  axvfs(1,j,ns) = 0.
                  do i = 2,narxt(j)
                    axlks(i,j,ns) = 0.
                    axhfs(i,j,ns) = 0.
                    axvfs(i,j,ns) = 0.
                  enddo
                enddo
                jsflag(ns) = 2
              endif
c
            enddo
          enddo
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make the 'd' set consistent with the ion exchanger part of the
c     ordinary set as it now exists. Do not modify other parts of the
c     'd' set. This action will leave, say, Na-Z in the 'd' set in
c     place of __-Z, which was the basis species on which the
c     exchanger substrate mass balance was defined. The only memory
c     of this definition is contained in the nbaspi array. Note
c     that the 'd' set reactions for exchangers will now include
c     the effects of eliminations of other basis species from the
c     active basis set.
c
c     The effects of eliminations on exchanger reactions could be
c     overcome by saving the exchanger reactions and associated data
c     prior to making the eliminations. However, new arrays would be
c     required to do that. It doesn't seem worthwhile at the present
c     time.
c
      call mdrgex(adhfs,adhfsd,adhfsx,advfs,advfsd,advfsx,
     $ axhfs,axhfsd,axhfsx,axlks,axlksd,axlksx,axvfs,axvfsd,axvfsx,
     $ cdrs,cdrsd,cdrsx,ipch,ipchmx,ipcv,ipcvmx,narxmx,nbasp,nbaspd,
     $ nbtmax,ndrs,ndrsd,ndrsx,ndrsmx,ndrsr,ndrsrd,ndrsrx,nern1,
     $ nern2,noutpt,nst,nstmax,ntprmx,nttyo)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
