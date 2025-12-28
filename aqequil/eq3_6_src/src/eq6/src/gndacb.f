      subroutine gndacb(cdac,cdacb,cdrs,eps100,imech,imchmx,
     $ jflag,nbasp,nbt,nbtmax,ndac,ndacb,ndact,ndactb,ndctmx,
     $ ndrs,ndrsmx,ndrsr,nrct,nrctmx,nstmax)
c
c     This subroutine computes the ndacb, ndactb, and cdacb arrays,
c     which are used to support the higher-order stiff ODE corrector.
c     Basically, these are the respective analogs of the ndac, ndact,
c     and cdac arrays. The latter arrays partially describe the
c     kinetic activity products of terms in rate expressions. The
c     ndac array contains the species indices of the species whose
c     activities appear in an activity product, the ndact array
c     contains the total numbers of such species for the kinetic
c     activity products, and the cdac array contains the respective
c     exponents for the species whose indices are given in the
c     ndac array. The species whose indices appear in the ndac
c     array need not be basis species. The ndacb, ndactb, and
c     cdacb arrays are transforms in which the species in ndacb
c     are all members of the active basis set. These arrays must
c     be adjusted if this set is changed by basis switching.
c
c     This subroutine is called by:
c
c       EQ6/path.f
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
      integer imchmx,nbtmax,ndctmx,ndrsmx,nrctmx,nstmax
c
      integer nbt,nrct
c
      integer imech(2,nrctmx),jflag(nstmax),nbasp(nbtmax),
     $ ndac(ndctmx,imchmx,2,nrctmx),ndacb(nbt,imchmx,2,nrct),
     $ ndact(imchmx,2,nrctmx),ndactb(imchmx,2,nrct),ndrs(ndrsmx),
     $ ndrsr(2,nstmax)
c
      real(8) cdac(ndctmx,imchmx,2,nrctmx),cdacb(nbt,imchmx,2,nrct),
     $ cdrs(ndrsmx)
      real(8) eps100
c
c----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j,k,n,nb,nn,nnn,ns,nse,nr1,nr2,nsi
c
      real(8) cxe
c
c----------------------------------------------------------------------
c
      do k = 1,nrct
        do j = 1,2
          do i = 1,imech(j,k)
c
c           Process the current kinetic activity product.
c
            nn = 0
            do nb = 1,nbt
              nse = nbasp(nb)
              if (jflag(nse) .eq. 0) then
                cxe = 0.
                do n = 1,ndact(i,j,k)
                  ns = ndac(n,i,j,k)
                  if (ns .eq. nse) then
c
c                   The current species in the kinetic activity product
c                   is the basis species presently being examined.
c
                    cxe = cxe + cdac(n,i,j,k)
                  elseif (jflag(ns) .eq. 30) then
c
c                   The current species in the activity product
c                   is a dependent species. Look for a dependence
c                   on the basis species presently being examined.
c
                    nr1 = ndrsr(1,ns)
                    nr2 = ndrsr(2,ns)
                    do nnn = nr1 + 1,nr2
                      nsi = ndrs(nnn)
                      if (nsi .eq. nse) then
                        cxe = cxe - cdac(n,i,j,k)*(cdrs(nnn)/cdrs(nr1))
                      endif
                    enddo
                  endif
c
                  if (abs(cxe) .ge. eps100) then
c
c                   The basis species presently being examined
c                   belongs in the current transformed kinetic
c                   activity product.
c
                    nn = nn + 1
                    ndacb(nn,i,j,k) = nse
                    cdacb(nn,i,j,k) = cxe
                  endif
                enddo
              endif
            enddo
            ndactb(i,j,k) = nn
            do nn = ndactb(i,j,k) + 1,nbt
              ndacb(nn,i,j,k) = 0
              cdacb(nn,i,j,k) = 0.
            enddo
c
          enddo
        enddo
      enddo
c
      end
