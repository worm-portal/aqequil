      subroutine gdlgxw(cdrs,cjbasp,cnufac,conc,dlogxw,eps100,
     $ ixbasp,jcsort,jflag,narn1,narn2,nbasp,nbt,nbtmax,nbw,
     $ ndrs,ndrsmx,ndrsr,nern1,nern2,noutpt,nstmax,nttyo,omega,
     $ xbar,xbarw)
c
c     This subroutine computes the array dlogxw, which contains the
c     partial derivatives:
c
c       d log x(w)/d log m(s')
c
c     where s' denotes a basis species other than water. Allowed
c     basis species may be aqueous or generic ion exchanger species.
c
c     This subroutine is called by:
c
c       EQLIB/matrix.f
c       EQ3NR/arrsim.f
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
      integer nbtmax,ndrsmx,nstmax
c
      integer noutpt,nttyo
c
      integer ixbasp(nbtmax),jcsort(nstmax),jflag(nstmax),nbasp(nbtmax),
     $ ndrs(ndrsmx),ndrsr(2,nstmax)
      integer narn1,narn2,nbt,nbw,nern1,nern2
c
      real*8 cdrs(ndrsmx),cjbasp(nbtmax),cnufac(nstmax),conc(nstmax),
     $ dlogxw(nbtmax),xbar(nstmax)
      real*8 eps100,omega,xbarw
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nb,nsc,nse,nss
c
      real*8 cw,cxec,dx,sumw,xx
c
      real*8 coefdr
c
c-----------------------------------------------------------------------
c
      if (nbw .eq. 0) then
        write (noutpt,1000)
        write (nttyo,1000)
 1000   format(/' * Error - (EQLIB/dxlgxw) Programming error trap:',
     $  /7x,'Water is not in the basis set, so the dlogxw array',
     $  /7x,'should not be calculated.')
        stop
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make sure that the "conc" of water is zero. This will be restored
c     at the end of this subroutine. This practice allows avoiding "if"
c     testing in the following loops. Note that there is no actual
c     element of the dlogxw array for water (basis index nbw). However,
c     there is a sum used in the denominator of the other array elements
c     that is based on water, and this resembles the sums accumulated
c     for the other elements. Therefore, dlogxw(nbw) will be used here
c     to accumulate that sum. It will be reset to zero after the sum
c     is used.
c
      cw = conc(narn1)
      conc(narn1) = 0.
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the sums needed to calculate the dlogxw array.
c     Store them in this array.
c
      do nb = 1,nbt
        nse = nbasp(nb)
        dlogxw(nb) = conc(nse)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do nss = narn1,narn2
        nsc = jcsort(nss)
        if (jflag(nsc) .eq. 30) then
c
          do nb = 1,nbt
            nse = nbasp(nb)
c
c           Calling sequence substitutions:
c             nsc for ns
c
            cxec = coefdr(cdrs,ndrs,ndrsmx,ndrsr,nse,nsc,nstmax)
c
c           Here cnufac(nsc) = -conc(nsc)/cdrs(ndrsr(1,nsc)), or
c           the negative of the molality of a dependent species
c           divided by its (intrinsically negative) reaction
c           coefficient in its associated reaction.
c
            dlogxw(nb) = dlogxw(nb) + cxec*cnufac(nsc)
          enddo
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the dlogxw array from the sums. The common denominator
c     (dx) for the elements corresponding to basis species other than
c     water is computed from the sum for water.
c
c     Do not protect here against a value of xbarw that is ridiculously
c     small, say on the order of 1.e-30. Such a value seems to be
c     necessary to obtaining good results in pre-Newton-Raphson
c     optimization. Only protect against a zero divide below, where
c     dx (derived from xbarw) appears in a denominator.
c
      xx = xbarw/omega
      sumw = dlogxw(nbw)
      dlogxw(nbw) = 0.
      dx = 1. + xx*sumw
c
c     Protect against the case in which dx is zero. This is designed
c     to prevent a fatal floating point exception in pre-Newton-Raphson
c     optimization.
c
      if (dx.ge.0. .and.  dx.lt.eps100) dx = eps100
      if (dx.lt.0. .and.  dx.gt.-eps100) dx = -eps100
c
      do nb = 1,nbt
        if (nb .ne. nbw) then
          dlogxw(nb) = -xx*dlogxw(nb)/dx
          if (ixbasp(nb) .ge. 1) then
c
c           Here nb refers to a basis species whose thermodynamic
c           activity is defined by its mole fraction.
c
            nse = nbasp(nb)
            dlogxw(nb) = dlogxw(nb)*cjbasp(nb)*(1.0 - xbar(nse))
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Restore the "conc" of water.
c
      conc(narn1) = cw
c
      end
