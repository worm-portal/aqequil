      subroutine intinx(ier,ixrn1,ixrn2,ncmpr,ncmpri,noutpt,npnxp,
     $ nptmax,nstmax,nttyo,nxicmx,nxti,nxtimx,umemi,uphase,usoli,
     $ uspec,xbar,xbari,xbarlg)
c
c     This subroutine interpets the solid solution compositions
c     specified on the input file.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
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
      integer nptmax,nstmax,nxicmx,nxtimx
c
      integer noutpt,nttyo
c
      integer ncmpr(2,nptmax),ncmpri(2,nxtimx),npnxp(nxtimx)
      integer ier,ixrn1,ixrn2,nxti
c
      character*48 uspec(nstmax)
      character*24 umemi(nxicmx),uphase(nptmax),usoli(nxtimx)
c
      real(8) xbar(nstmax),xbari(nxicmx),xbarlg(nstmax)
c
c----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3,np,nr1,nr1i,nr2,nr2i,ns,nxi,nxic
c
      integer ilnobl
c
      character*24 uxp,uxs
c
      real(8) xx
c
      real(8) tlg
c
c----------------------------------------------------------------------
c
      ier = 0
      if (nxti .le. 0) go to 999
c
      do nxi = 1,nxti
        npnxp(nxi) = 0
        uxp = usoli(nxi)
c
        do np = ixrn1,ixrn2
          if (uxp(1:24) .eq. uphase(np)(1:24)) go to 100
        enddo
c
        ier = ier + 1
        j2 = ilnobl(uxp)
        write (noutpt,1000) uxp(1:j2)
        write (nttyo,1000) uxp(1:j2)
 1000   format(/' * Error - (EQ3NR/intinx) The solid solution ',a,
     $   /7x,'is on the input file, but not on the data file.')
        go to 130
c
  100   npnxp(nxi) = np
        nr1i = ncmpri(1,nxi)
        nr2i = ncmpri(2,nxi)
        nr1 = ncmpr(1,np)
        nr2 = ncmpr(2,np)
c
        do nxic = nr1i,nr2i
          uxs = umemi(nxic)
c
          do ns = nr1,nr2
            if (uxs(1:24) .eq. uspec(ns)(1:24)) go to 110
          enddo
c
          ier = ier + 1
          j2 = ilnobl(uxs)
          j3 = ilnobl(uxp)
          write (noutpt,1010) uxs(1:j2),uxp(1:j3)
          write (nttyo,1010) uxs(1:j2),uxp(1:j3)
 1010     format(/' * Error - (EQ3NR/intinx) The solid solution',
     $    /7x,'species ',a,' (',a,') is on the'
     $    /7x,'input file, but not on the data file.')
          go to 120
c
  110     continue
          xx = xbari(nxic)
          xbar(ns) = xx
          xbarlg(ns) = tlg(xx)
  120     continue
        enddo
      enddo
c
  130 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
