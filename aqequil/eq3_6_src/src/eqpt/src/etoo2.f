      subroutine etoo2(cdrsi,dhfe,dhfs,dvfe,dvfs,ipch,ipchmx,ipcv,
     $ ipcvmx,narxmx,narxt,nbtmx1,ndrsts,ns,ntprmx,ntprt,udrsi,
     $ xhfe,xhfs,xlke,xlks,xvfe,xvfs)
c
c     This subroutine converts the reaction for the ns-th species from
c     one written in terms of e- to one written in terms of O2(g).
c     The "Eh" reaction is used to do this.
c
c     This subroutine is called by:
c       EQPT/pcraq.f
c       EQPT/pcrsg.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       cdrsi  = array of coefficients for the original reaction
c                  associated with the ns-th species
c       nbtmx1 = the maximum number of basis species plus 1
c       ndrsts = number of species in the reaction
c       ns     = index of the species whose reaction is to be rewritten
c       udrsi  = names of species in the original reaction
c       xlke   = the array of log K values for the "Eh reaction"
c       xlks   = the array of log K values for the original reaction
c
c     Output:
c
c       cdrsi  = array of coefficients for a reaction, rewritten
c                  in terms of e- instead of O2(g)
c       udrsi  = names of species in the rewritten reaction
c       xlks   = the array of log K values for the rewritten reaction
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipchmx,ipcvmx,narxmx,nbtmx1,ntprmx
c
      integer narxt(ntprmx)
c
      integer ipch,ipcv,ndrsts,ns,ntprt
c
      character*24 udrsi(nbtmx1)
c
      real*8 cdrsi(nbtmx1),dhfe(narxmx,ntprmx,ipchmx),
     $ dhfs(narxmx,ntprmx,ipchmx,nbtmx1),dvfe(narxmx,ntprmx,ipcvmx),
     $ dvfs(narxmx,ntprmx,ipcvmx,nbtmx1),xhfe(narxmx,ntprmx),
     $ xhfs(narxmx,ntprmx,nbtmx1),xlke(narxmx,ntprmx),
     $ xlks(narxmx,ntprmx,nbtmx1),xvfe(narxmx,ntprmx),
     $ xvfs(narxmx,ntprmx,nbtmx1)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ipc,kdrsts,n,nelect,nhydr,ntpr,nwater
c
      real*8 factor
c
c-----------------------------------------------------------------------
c
c     Search for e-, H+, and H2O in the reaction.
c
      nelect = 0
      nhydr = 0
      nwater = 0
c
      do n = 2,ndrsts
        if (udrsi(n)(1:3) .eq. 'e- ') nelect = n
        if (udrsi(n)(1:3) .eq. 'H+ ') nhydr = n
        if (udrsi(n)(1:4) .eq. 'H2O ') nwater = n
      enddo
c
      if (nelect .eq. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Re-write the reaction, using the Eh reaction:
c
c       H2O = 4H+ + 4e- + O2(g)
c
      factor = -cdrsi(nelect)/4.
c
      udrsi(nelect) = 'O2(g)'
      cdrsi(nelect) = factor
c
      if (nhydr .gt. 0) then
        cdrsi(nhydr) = cdrsi(nhydr) + 4.*factor
      else
        ndrsts = ndrsts + 1
        nhydr = ndrsts
        udrsi(nhydr) = 'H+'
        cdrsi(nhydr) = 4.*factor
      endif
c
      if (nwater .gt. 0) then
        cdrsi(nwater) = cdrsi(nwater) - 2.*factor
      else
        ndrsts = ndrsts + 1
        nwater = ndrsts
        udrsi(nwater) = 'H2O'
        cdrsi(nwater) = -2.*factor
      endif
c
c     Clear any zeroes in the reaction coefficients.
c
      kdrsts = ndrsts
      do n = 1,ndrsts
        if (abs(cdrsi(n)) .le. 1.e-6) then
          if (n .le. kdrsts - 1) then
            do i = n,kdrsts - 1
              cdrsi(i) = cdrsi(i + 1)
              udrsi(i) = udrsi(i + 1)
            enddo
          endif
          kdrsts = kdrsts - 1
        endif
        if (n .ge. kdrsts) go to 150
      enddo
  150 ndrsts = kdrsts
c
c     Modify the log K values.
c
      do ntpr = 1,ntprt
        do n = 1,narxt(ntpr)
          if (xlks(n,ntpr,ns).lt.9999999. .and.
     $      xlke(n,ntpr).lt.9999999.) then
            xlks(n,ntpr,ns) = xlks(n,ntpr,ns) + factor*xlke(n,ntpr)
          else
            xlks(n,ntpr,ns) = 9999999.
          endif
        enddo
      enddo
c
      if (ipch .ge. 0) then
        do ntpr = 1,ntprt
          do n = 1,narxt(ntpr)
            if (xhfs(n,ntpr,ns).lt.9999999. .and.
     $        xhfe(n,ntpr).lt.9999999.) then
              xhfs(n,ntpr,ns) = xhfs(n,ntpr,ns) + factor*xhfe(n,ntpr)
            else
              xhfs(n,ntpr,ns) = 9999999.
            endif
          enddo
        enddo
c
        do ipc = 1,ipch
          do ntpr = 1,ntprt
            do n = 1,narxt(ntpr)
              if (dhfs(n,ntpr,ipc,ns).lt.9999999. .and.
     $          dhfe(n,ntpr,ipc).lt.9999999.) then
                dhfs(n,ntpr,ipc,ns) = dhfs(n,ntpr,ipc,ns)
     $          + factor*dhfe(n,ntpr,ipc)
              else
                dhfs(n,ntpr,ipc,ns) = 9999999.
              endif
            enddo
          enddo
        enddo
      endif
c
      if (ipcv .ge. 0) then
        do ntpr = 1,ntprt
          do n = 1,narxt(ntpr)
            if (xvfs(n,ntpr,ns).lt.9999999. .and.
     $        xvfe(n,ntpr).lt.9999999.) then
              xvfs(n,ntpr,ns) = xvfs(n,ntpr,ns) + factor*xvfe(n,ntpr)
            else
              xvfs(n,ntpr,ns) = 9999999.
            endif
          enddo
        enddo
c
        do ipc = 1,ipcv
          do ntpr = 1,ntprt
            do n = 1,narxt(ntpr)
              if (dvfs(n,ntpr,ipc,ns).lt.9999999. .and.
     $          dvfe(n,ntpr,ipc).lt.9999999.) then
                dvfs(n,ntpr,ipc,ns) = dvfs(n,ntpr,ipc,ns)
     $          + factor*dvfe(n,ntpr,ipc)
              else
                dvfs(n,ntpr,ipc,ns) = 9999999.
              endif
            enddo
          enddo
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
