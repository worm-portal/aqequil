c eqldv.h
c
c     Dimensioning variables.
c
c-----------------------------------------------------------------------
c
c     Those used in both EQ3NR and EQ6 modes (fixed/variable mass of
c     solvent water).
c
      integer iapxmx,ibpxmx,ietmax,iktmax,ipbtmx,ipchmx,ipcvmx,
     $ jetmax,jpfcmx,jsomax,ketmax,kmax,napmax,narxmx,natmax,
     $ nazmmx,nazpmx,nbtmax,nbt1mx,nctmax,netmax,ndrsmx,nessmx,
     $ ngtmax,nltmax,nmtmax,nmutmx,nmxmax,nodbmx,nopgmx,noprmx,
     $ noptmx,nptmax,nsltmx,nstmax,nstsmx,nsxmax,ntidmx,ntitmx,
     $ ntf1mx,ntf2mx,ntfxmx,ntprmx,nvetmx,nxmdmx,nxtmax
c
      common /eqldv/ iapxmx,ibpxmx,ietmax,iktmax,ipbtmx,ipchmx,
     $ ipcvmx,jetmax,jpfcmx,jsomax,ketmax,kmax,napmax,narxmx,
     $ natmax,nazmmx,nazpmx,nbtmax,nbt1mx,nctmax,netmax,ndrsmx,
     $ nessmx,ngtmax,nltmax,nmtmax,nmutmx,nmxmax,nodbmx,nopgmx,
     $ noprmx,noptmx,nptmax,nsltmx,nstmax,nstsmx,nsxmax,ntidmx,
     $ ntitmx,ntf1mx,ntf2mx,ntfxmx,ntprmx,nvetmx,nxmdmx,nxtmax
c
c-----------------------------------------------------------------------
c
c     Those used in only EQ3NR mode (fixed mass of solvent water).
c
      integer njfmax,nxicmx,nxtimx
c
      common /eqldv1/ njfmax,nxicmx,nxtimx
c
c-----------------------------------------------------------------------
c
c     Those used in only EQ6 mode (variable mass of solvent water).
c
      integer imchmx,ndctmx,nertmx,nffgmx,nordmx,npetmx,nprpmx,nprsmx,
     $ nptkmx,nrctmx,nrd1mx,nsetmx,nttkmx,nsrtmx,nxopmx,nxpemx,nxrtmx,
     $ nsscmx
c
      common /eqldv2/ imchmx,ndctmx,nertmx,nffgmx,nordmx,npetmx,nprpmx,
     $ nprsmx,nptkmx,nrctmx,nrd1mx,nsetmx,nttkmx,nsrtmx,nxopmx,nxpemx,
     $ nxrtmx,nsscmx
c
c     End of INCLUDE file eqldv.h
c-----------------------------------------------------------------------
