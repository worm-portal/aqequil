      subroutine rdazp(azero,insgf,nazt,naztmx,ndat0s,nerr,noutpt,
     $ nttyo,uazp)
c
c     This subroutine reads hard core diamaters and related parameters
c     used in the B-dot equation from the DATA0 file. This data
c     consists of aqueous species names, hard core diameter
c     values, and an integer flag to indicate whether a species that
c     is electrically neutral should be considered polar (set log
c     gamma = 0) or non-polar (log gamma = log gamma for CO2(aq)).
c     The data are returned in the uazp, azero, and insgf arrays. They
c     are later written on the DATA1 and DATA1F files by EQPT/wrazp.f.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ndat0s = unit number of the stripped DATA0 file
c
c     Principal output:
c
c       nazt   = the number of specified hard core diameters
c       uazp   = array of aqueous species names used to specify
c                  hard core diamters on the data file
c       azero  = array of corresponding hard core diameters
c       insgf  = array of corresponding neutral species
c                  activity coefficient flags
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
      integer ndat0s,noutpt,nttyo
c
      integer nazt,nerr
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
      integer j2,j3,j4
c
      integer ilnobl
c
      character*80 uline
      character*24 ux24
      character*8 ubdot,uterm,ux8
c
c-----------------------------------------------------------------------
c
      data ubdot / 'bdot par' /
      data uterm / '+-------' /
c
c-----------------------------------------------------------------------
c
      nazt = 0
c
c     Advance to bdot header.
c
  100 continue
      read (ndat0s,1000,end=990,err=995) uline
 1000 format(a)
      ux8 = uline(1:8)
      if (ux8(1:8) .ne. ubdot(1:8)) go to 100
c
c     Skip the terminator line.
c
      read (ndat0s,1000,end=990,err=995) uline
c
c     Read in the azero ('bdot') data.
c
  110 continue
      read (ndat0s,1000,end=990,err=995) uline
      ux8 = uline(1:8)
c
      if (ux8(1:8) .eq. uterm(1:8)) then
c
c       Have found the block terminator. Skip past the element
c       header.
c
        read (ndat0s,1000,end=990,err=995) uline
c
        go to 999
      endif
c
c     Found another entry. Read the data from the line.
c
      nazt = nazt + 1
      read(uline,1010,err=998) uazp(nazt),azero(nazt),insgf(nazt)
 1010 format(a24,7x,f7.1,4x,i2)
c
      j2 = ilnobl(uazp(nazt))
      if (j2 .le. 0) then
        write (ux8,'(i5)') nazt
        call lejust(ux8)
        j3 = ilnobl(ux8)
        write (noutpt,1030) ux8(1:j3)
        write (nttyo,1030) ux8(1:j3)
 1030   format(/' * Error - (EQPT/rdazp) Have encountered a blank',
     $  ' species name on line ',a,/7x,'of the block of hard core',
     $  ' diameter values for aqueous species.')
        if (nazt .gt. 1) then
          ux24 = uazp(nazt - 1)
          if (ux24(1:7) .ne. '<blank>') then
            j4 = ilnobl(ux24)
            write (noutpt,1040) ux24(1:j4)
            write (nttyo,1040) ux24(1:j4)
 1040       format(7x,'This line follows the one for ',a,'.')
          endif
        endif
        uazp(nazt) = '<blank>'
        nerr = nerr + 1
      endif
c
      go to 110
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 write (noutpt,2000)
      write (nttyo,2000)
 2000 format(/' * Error - (EQPT/rdazp) Unexpectedly encountered',
     $ /7x,'end-of-file while reading the DATA0 file.')
c
      write (noutpt,2010) nazt
      write (nttyo,2010) nazt
 2010 format(7x,'The value of the local block line counter is ',i4,'.')
      j2 = ilnobl(uline)
      if (j2 .gt. 0) then
        j2 = min(j2,70)
        write (noutpt,2030) uline(1:j2)
        write (nttyo,2030) uline(1:j2)
 2030   format(7x,'The last line read was the following:',
     $  /7x,'"',a,'"')
      endif
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  995 write (noutpt,2040)
      write (nttyo,2040)
 2040 format(/' * Error - (EQPT/rdazp) Encountered a read format',
     $ /7x,'error while reading the DATA0 file.')
c
      write (noutpt,2010) nazt
      write (nttyo,2010) nazt
      j2 = ilnobl(uline)
      if (j2 .gt. 0) then
        j2 = min(j2,70)
        write (noutpt,2030) uline(1:j2)
        write (nttyo,2030) uline(1:j2)
      endif
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  998 write (noutpt,2050)
      write (nttyo,2050)
 2050 format(/' * Error - (EQPT/rdazp) Encountered a read format error',
     $ /7x,'while reading data from lines read from the DATA0 file.')
      j2 = ilnobl(uline)
      if (j2 .gt. 0) then
        j2 = min(j2,70)
        write (noutpt,2060) uline(1:j2)
        write (nttyo,2060) uline(1:j2)
 2060   format(7x,'The line with the problem is:',
     $  /7x,'"',a,'"')
      endif
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
