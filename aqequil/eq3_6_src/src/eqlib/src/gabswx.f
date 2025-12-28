      subroutine gabswx(beta,ibswx,iindx1,kbt,kmax,nbt,nbtmax)
c
c     This subroutine supports automatic basis switching as a
c     pre-Newton-Raphson optimization technique. It resolves conflicts
c     in the ibswx array, which contains candidate switches. This array
c     could call for the same non-basis species to be switched with
c     more than one basis species. The approach here is to use the size
c     of the associated relative mass balance residual to any resolve
c     conflicts. This subroutine is somewhat similar to EQLIB/gbfac.f,
c     which resolves conflicts affecting continued fraction
c     calculations.
c
c     This subroutine is called by:
c
c       EQLIB/absswa.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ibswx = array of indices of species not in the active
c                 master set that are candidates for switching into
c                 the basis set
c       beta  = array of normalized Newton-Raphson residual functions
c       kbt   = number of species in the active basis set
c
c     Principal output:
c
c       ibswx = input array modified to remove conflicts.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer kmax,nbtmax
c
      integer ibswx(nbtmax),iindx1(kmax)
c
      integer kbt,nbt
c
      real(8) beta(kmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer kcol,kk,krow,nb,nbig,nb1,nb2,ns1,ns2
c
      real(8) bbig,frac
c
c-----------------------------------------------------------------------
c
c     Find the largest residual among potential basis switching cases.
c
      nbig = 0
      bbig = 0.
      do kcol = 1,kbt
        nb = iindx1(kcol)
        if (ibswx(nb) .gt. 0) then
          if (beta(kcol) .gt. bbig) then
            nbig = nb
            bbig = beta(kcol)
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Eliminate candidates not involving largest residuals.
c
      if (bbig .le. 0.5) then
c
c       If the largest residual is not very large, wipe out all other
c       proposed switches.
c
        do nb = 1,nbt
          if (nb .ne. nbig) ibswx(nb) = 0
        enddo
        go to 999
      else
c
c       If the largest residual is rather large, wipe out all other
c       proposed switches with residuals that are not within two
c       orders of magnitude of the largest residual.
c
        do kcol = 1,kbt
          nb = iindx1(kcol)
          if (ibswx(nb) .gt. 0) then
            frac = beta(kcol)/bbig
            if (frac .le. 1.e-2) ibswx(nb) = 0
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Eliminate any remaining conflicts.
c
      do krow = 1,kbt - 1
        nb1 = iindx1(krow)
        ns1 = ibswx(nb1)
        if (ns1 .gt. 0) then
          kk = krow + 1
          do kcol = kk,kbt
            nb2 = iindx1(kcol)
            ns2 = ibswx(nb2)
            if (ns2 .eq. ns1) then
              if (beta(krow) .gt. beta(kcol)) then
                ibswx(nb2) = 0
              else
                ibswx(nb1) = 0
                go to 25
              endif
            endif
          enddo
        endif
   25   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
