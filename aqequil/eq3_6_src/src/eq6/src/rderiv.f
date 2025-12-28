      subroutine rderiv(akmat0,drer0,drir0,fdri0,fdrr0,jreac,nord,
     $ nrct,nrctmx,nrd1mx)
c
c     This subroutine transforms finite differences of the rates of
c     irreversible reactions into the corresponding derivatives.
c     Note that (drer0) = (akmat0)(fdrr0). The inverse rate is treated
c     similarly; i.e., (drir0) = (akmat0)(fdri0).
c
c     Compare with:
c
c       EQ6/aderiv.f
c       EQ6/bderiv.f
c       EQ6/pderiv.f
c       EQ6/zderiv.f
c
c     This subroutine is called by:
c
c       EQ6/chksti.f
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
      integer nrctmx,nrd1mx
c
      integer jreac(nrctmx)
c
      integer nord,nrct
c
      real*8 akmat0(nrd1mx,nrd1mx),drer0(nrd1mx,nrctmx),
     $ drir0(nrd1mx),fdri0(nrd1mx),fdrr0(nrd1mx,nrctmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer k,n,nrc
c
      real*8 dxx
c
c-----------------------------------------------------------------------
c
c     Calculate the derivatives of the inverse rate.
c
      do n = 1,nord
        dxx = 0.
        do k = n,nord
          dxx = dxx + fdri0(k)*akmat0(n,k)
        enddo
        drir0(n) = dxx
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Calculate the derivatives of the rates of irreversible reactions.
c
      do nrc = 1,nrct
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1)  then
          do n = 1,nord
            dxx = 0.
            do k = n,nord
              dxx = dxx + fdrr0(k,nrc)*akmat0(n,k)
            enddo
            drer0(n,nrc) = dxx
          enddo
        endif
      enddo
c
      end
