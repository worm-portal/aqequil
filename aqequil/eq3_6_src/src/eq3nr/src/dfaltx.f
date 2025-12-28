      subroutine dfaltx(itermx,rho,scamas,tdspkg,tdspl,tolbt,
     $ toldl,tolspf)
c
c     This subroutine sets the defaults for various run parameters. In
c     some cases, it forces the parameters to take on certain values,
c     or to fall in certain ranges.
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
      integer itermx
c
      real(8) rho,scamas,tdspkg,tdspl,tolbt,toldl,tolspf
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c       None
c
c-----------------------------------------------------------------------
c
      if (itermx .le. 0) itermx = 200
      if (rho .le. 0.) rho = 1.0
      if (tdspkg .lt. 0.) tdspkg = 0.
      if (tdspl .lt. 0.) tdspl = 0.
      if (scamas .le. 0.) scamas = 1.0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (tolbt .le. 0.) tolbt = 1.e-6
      if (tolbt .lt. 1.e-10) tolbt = 1.e-10
      if (tolbt .gt. 1.e-2) tolbt = 1.e-2
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (toldl .le. 0.) toldl = 1.e-6
      if (toldl .lt. 1.e-10) toldl = 1.e-10
      if (toldl .gt. 1.e-2) toldl = 1.e-2
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (tolspf .le. 0.) tolspf = 0.0005
      if (tolspf .lt. 0.01) tolspf = 0.01
      if (tolspf .gt. 0.00005) tolspf = 0.00005
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
