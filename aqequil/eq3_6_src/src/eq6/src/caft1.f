      subroutine caft1(afrc1,aft1,nrct,nrctmx,rrelr1)
c
c     This subroutine computes the total affinity (aft1).
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
      integer nrctmx
c
      integer nrct
c
      real(8) aft1
      real(8) afrc1(nrctmx),rrelr1(nrctmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nrc
c
c-----------------------------------------------------------------------
c
      aft1 = 0.
      do nrc = 1,nrct
        if (rrelr1(nrc) .ne. 0.) then
          if (afrc1(nrc) .lt. 9999999.) then
            aft1 = aft1 + abs(afrc1(nrc)*rrelr1(nrc))
          else
c
c           Have a reactant with an "infinite" affinity. Set the
c           total affinity to "infinity" also.
c
            aft1 = 9999999.
            go to 100
          endif
        endif
      enddo
  100 continue
c
      end
