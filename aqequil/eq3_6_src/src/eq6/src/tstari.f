      subroutine tstari(afrc1,alphar,betar,delxi,iodb,modr,nodbmx,
     $ noutpt,nrct,nrctmx,nrct1,nsscmx,qodeok,rirec1,rirecp,rrelr1,
     $ rrelrp,sscrew,time1,tistrt,ureac)
c
c     This subroutine tests the accuracy of the rate law integration
c     by comparing computed rates with predicted values.
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
c       qodeok = logical flag, .true. if the rate law integration
c                  tolerances have been satisfied.
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nodbmx,nrctmx,nsscmx
c
      integer noutpt
c
      integer iodb(nodbmx)
c
      integer nrct,nrct1
c
      logical qodeok
c
      character*24 ureac(nrctmx)
c
      real(8) afrc1(nrctmx),alphar(nrct1),betar(nrct1),modr(nrctmx),
     $ rrelr1(nrctmx),rrelrp(nrctmx),sscrew(nsscmx)
c
      real(8) delxi,rirec1,rirecp,time1,tistrt
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,nrc
c
      integer ilnobl
c
      logical qbad
c
      real(8) atx,etx,rdelmd,rdelti,ss4
c
c-----------------------------------------------------------------------
c
      qodeok = .true.
      qbad = .false.
      ss4 = sscrew(4)**2
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do nrc = 1,nrct
c
c       The error in a relative rate is acceptable is any of three
c       tests is satisfied.
c
        qbad = .false.
c
c       Test the fractional error in the relative rate.
c
        if (abs(betar(nrc)) .gt. sscrew(4)) then
c
c         Test the absolute error in the relative rate.
c
          if (abs(alphar(nrc)) .gt. ss4) then
            if (abs(modr(nrc)) .le. 0.) then
              qbad = .true.
            else
c
c             Test the error in the moles destroyed for this step
c             against the total moles destroyed.
c
              rdelmd = (alphar(nrc)*delxi)/abs(modr(nrc))
c
              if (abs(rdelmd) .gt. sscrew(3)) qbad = .true.
            endif
          endif
        endif
c
        if (qbad) then
          qodeok = .false.
          if (iodb(2) .ge. 1) then
            j2 = ilnobl(ureac(nrc))
            write (noutpt,1010) ureac(nrc)(1:j2),rrelr1(nrc),
     $      rrelrp(nrc),betar(nrc),afrc1(nrc)
 1010       format(/3x,'Relative rate for ',a,':',
     $      /5x,'Calculated value= ',1pe12.5,
     $      /5x,' Predicted value= ',e12.5,
     $      /5x,'Fractional error= ',e12.5,
     $      /5x,'Affinity= ',e12.5,/)
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The error in the inverse rate is acceptable is either of two
c     tests is satisfied.
c
      qbad = .false.
c
c     Test the fractional error in the inverse rate.
c
      if (abs(betar(nrct1)) .gt. sscrew(4)) then
c
c       Test the absolute error in the time.
c
        atx = alphar(nrct1)*delxi
cXX
        if (abs(atx) .gt. 0.) then
          etx = time1 - tistrt
          if (etx.eq. 0.) then
            qbad = .true.
          else
c
            rdelti = atx/etx
c
c           Test the error in the time against the elapsed time.
c
            if (abs(rdelti) .gt. sscrew(3)) qbad = .true.
          endif
        endif
      endif
c
      if (qbad) then
        qodeok = .false.
        if (iodb(2) .ge. 1) then
          write (noutpt,1000) rirec1,rirecp,betar(nrct1),rdelti
 1000     format(/3x,'Inverse rate:',
     $    /5x,'Calculated value= ',1pe12.5,
     $    /5x,' Predicted value= ',e12.5,
     $    /5x,'Fractional error= ',e12.5,
     $    /5x,'Fractional error in time= ',e12.5,/)
        endif
      endif
c
      end
