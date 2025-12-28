      subroutine wrazp(azero,insgf,nazt,naztmx,ndata1,ndat1f,
     $ noutpt,nttyo,uazp)
c
c     This subroutine writes on the DATA1 and DATA1F files the "bdot"
c     data read from the DATA0 file by EQPT/rdazp.f.
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ndata1 = unit number of the DATA1 file
c       ndat1f = unit number of the DATA1F file
c       nazt   = the number of elements in the uazp array
c       uazp   = array containing lines of data
c       azero  = array of corresponding hard core diameters
c       insgf  = array of corresponding neutral species
c                  activity coefficient flags
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
      integer naztmx
c
      integer ndata1,ndat1f,noutpt,nttyo
c
      integer nazt
c
      integer insgf(naztmx)
c
      character*24 uazp(naztmx)
c
      real*8 azero(naztmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ineu,naz
c
      character*72 uterm,utermc
      character*24 uendit
c
      real*8 zero
c
c-----------------------------------------------------------------------
c
      data uendit / 'endit.' /
c
c-----------------------------------------------------------------------
c
      uterm(1:48) = '+-----------------------------------------------'
      uterm(49:72) = '------------------------'
      utermc = uterm
      utermc(1:1) = '*'
c
c     Write the azero ('bdot') data on the DATA1 and DATA1F files.
c
      do naz = 1,nazt
        write (ndata1) uazp(naz),azero(naz),insgf(naz)
        write (ndat1f,1010) uazp(naz),azero(naz),insgf(naz)
 1010   format(a24,2x,f7.1,2x,i2)
      enddo
c
c     Write the block terminator.
c
      zero = 0.
      ineu = 0
      write (ndata1) uendit,zero,ineu
      write (ndat1f,1010) uendit,zero,ineu
      write (ndat1f,1020) utermc(1:72)
 1020 format(a)
c
      end
