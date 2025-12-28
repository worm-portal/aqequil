      subroutine gakey(ndat0s,noutpt,nttyo,uakey)
c
c     This suboutine scans the DATA0 file to determine the aqueous
c     species activity coefficient model (e.g., Pitzer, Simple
c     Extended Debye-Huckel) associated with this file.
c
c     This suboutine is called by:
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
c       uakey  = the type of data file ("SEDH" or "Pitzer")
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ndat0s,noutpt,nttyo
c
      character(len=8) uakey
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,icount,n
c
      character(len=80) uline
c
c-----------------------------------------------------------------------
c
c     Determine whether or not this is a Pitzer-type data file.
c     Look for patterns indicating the presence of Pitzer interaction
c     coefficients.
c
      icount = 0
      do n = 1,500
        read (ndat0s,1000,end=100,err=995) uline
 1000   format(a)
        i = index(uline,'beta0 =')
        icount = icount + i
        i = index(uline,'beta1 =')
        icount = icount + i
        i = index(uline,'beta2 =')
        icount = icount + i
        i = index(uline,'cphi =')
        icount = icount + i
        i = index(uline,'Beta0 =')
        icount = icount + i
        i = index(uline,'Beta1 =')
        icount = icount + i
        i = index(uline,'Beta2 =')
        icount = icount + i
        i = index(uline,'Cphi =')
        icount = icount + i
        i = index(uline,'beta(0)')
        icount = icount + i
        i = index(uline,'beta(1)')
        icount = icount + i
        i = index(uline,'beta(2)')
        icount = icount + i
        i = index(uline,'Cphi')
        icount = icount + i
        i = index(uline,'C(phi)')
        icount = icount + i
        if (icount .ge. 24) go to 100
      enddo
  100 continue
c
      if (icount .ge. 8) then
        uakey = 'Pitzer'
      else
        uakey = 'SEDH'
      endif
c
      rewind(ndat0s)
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write a message for any read error other than end-of-file.
c
  995 write (noutpt,2010)
      write (nttyo,2010)
 2010 format(/' * Error - (EQPT/gakey) Encountered a read format',
     $ ' error while',/7x,'scanning the DATA0 file to determine',
     $ ' the associated aqueous',/7x,'activity coefficient model.')
      stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
