      subroutine prtmip(acflg,actlg,conclg,ctb,nbaspd,nbt,nbtmax,
     $ nelect,nhydr,nhydx,noutpt,nstmax,uspec,zchar)
c
c     This subroutine prints  tables of the mean ionic activities and
c     activity coefficients.
c
c     This subroutine is called by:
c
c       EQ3NR/scripx.f
c       EQ6/scripz.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       acflg  = array of log activity coefficients of species
c       actlg  = array of log activities of species
c       conclg = array of log concentrations of species
c       ctb    = array of total molalities of basis species
c       nbt    = the number of species in the basis set
c       nelect = index of the fictive species aqueous e-
c       nhydr  = index of the species aqueous H+
c       nhydx  = index of the species aqueous OH-
c       uspec  = array of species names
c       zchar  = array of species electrical charge numbers
c
c     Principal output:
c
c       None
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbtmax,nstmax
c
      integer noutpt
c
      integer nbaspd(nbtmax)
      integer nbt,nelect,nhydr,nhydx
c
      character*48 uspec(nstmax)
c
      real*8 acflg(nstmax),actlg(nstmax),conclg(nstmax),ctb(nbtmax),
     $ zchar(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3,nb1,nb2,ns1,ns2
c
      integer ilnobl
c
      real*8 apmlg,cpmlg,cpmslg,ctb1,ctb2,gpm,gpmlg,gpms,gpmslg,zxa,
     $ zx1,zx2
c
      real*8 texp,tlg
c
c-----------------------------------------------------------------------
c
c     Compute and print a table of stoichiometric ionic properties.
c
      write (noutpt,1000)
 1000 format(//16x,'--- Stoichiometric Mean Ionic Properties ---')
      write (noutpt,1010)
 1010 format(/3x,'Species',18x,'Log a(+/-)',2x,'Log m(+/-)',2x,
     $ 'Log gamma(+/-)',2x,'gamma(+/-)',/)
c
      do nb1 = 1,nbt
        ns1 = nbaspd(nb1)
        zx1 = zchar(ns1)
c
        ctb1 = ctb(nb1)
        if (ctb1 .lt. 0.) then
          if (ns1 .eq. nhydx) then
            ns1 = nhydr
            zx1 = zchar(nhydr)
            ctb1 = -ctb1
          endif
        endif
c
        if (zx1 .le. 0.) go to 120
        if (ctb1 .le. 0.) go to 120
c
        do nb2 = 1,nbt
          ns2 = nbaspd(nb2)
          zx2 = zchar(ns2)
c
          ctb2 = ctb(nb2)
          if (ctb2 .lt. 0.) then
            if (ns2 .eq. nhydr) then
              ns2 = nhydx
              zx2 = zchar(nhydx)
              ctb2 = -ctb2
            endif
          endif
c
          if (zx2 .ge. 0.) go to 110
          if (ctb2 .le. 0.) go to 110
          if (ns2 .eq. nelect) go to 110
c
          zxa = zx1 - zx2
          apmlg = (zx1*actlg(ns2) - zx2*actlg(ns1))/zxa
c
          cpmslg = (zx1*tlg(ctb2) - zx2*tlg(ctb1))/zxa
c
          gpmslg = apmlg - cpmslg
          gpms = texp(gpmslg)
c
          j2 = ilnobl(uspec(ns1)(1:24))
          j3 = ilnobl(uspec(ns2)(1:24))
c
  100     if ((j2 + j3) .gt. 24) then
            if (j2 .gt. j3) then
              j2 = j2 - 1
            else
              j3 = j3 - 1
            endif
            go to 100
          endif
c
          write (noutpt,1020) uspec(ns1)(1:j2),uspec(ns2)(1:j3),
     $    apmlg,cpmslg,gpmslg,gpms
 1020     format(1x,a,'/',a,t29,f9.4,3x,f9.4,4x,f9.4,5x,1pe11.4)
  110     continue
        enddo
c
  120   continue
      enddo
c
      write (noutpt,1030)
 1030 format(/3x,'The stoichiometric mean molalities and activity',
     $ ' coefficients given',/3x,'above are consistent with the',
     $ ' sensible composition of the',/3x,'aqueous solution.')
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute and print a table of ionic properties.
c
      write (noutpt,1100)
 1100 format(//21x,'--- Mean Ionic Properties ---')
      write (noutpt,1010)
c
      do nb1 = 1,nbt
        ns1 = nbaspd(nb1)
        zx1 = zchar(ns1)
c
        if (zx1 .le. 0.) go to 220
        if (conclg(ns1) .le. -99999.) go to 220
c
        do nb2 = 1,nbt
          ns2 = nbaspd(nb2)
          zx2 = zchar(ns2)
c
          if (zx2 .ge. 0.) go to 210
          if (conclg(ns2) .le. -99999.) go to 210
          if (ns2 .eq. nelect) go to 210
c
          zxa = zx1 - zx2
          apmlg = (zx1*actlg(ns2) - zx2*actlg(ns1))/zxa
c
          cpmlg = (zx1*conclg(ns2) - zx2*conclg(ns1))/zxa
c
          gpmlg = (zx1*acflg(ns2) - zx2*acflg(ns1))/zxa
          gpm = texp(gpmlg)
c
          j2 = ilnobl(uspec(ns1)(1:24))
          j3 = ilnobl(uspec(ns2)(1:24))
c
  200     if ((j2 + j3) .gt. 32) then
            if (j2 .gt. j3) then
              j2 = j2 - 1
            else
              j3 = j3 - 1
            endif
            go to 200
          endif
c
          write (noutpt,1020) uspec(ns1)(1:j2),uspec(ns2)(1:j3),
     $    apmlg,cpmlg,gpmlg,gpm
  210     continue
        enddo
c
  220   continue
      enddo
c
      write (noutpt,1110)
 1110 format(/3x,'The mean molalities and activity coefficients given',
     $ ' above are',/3x,'consistent with the speciation in the model',
     $ ' employed.',/)
c
      end
