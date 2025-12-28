      module mod6xf
c
c     This module contains standard state thermodynamic data and
c     standard grid pressure data. This excludes interpolating
c     polynomial coefficient arrays belonging to the 'a' set (e.g.,
c     axhfsa, axlksa, axvfsa) which includes all such data read
c     from the data file. Included here are the "compressed"
c     equivalents.,
c
c     Pressure parameters.
c
      real(8) prehw,presg
c
      real(8), dimension(:,:), allocatable :: aprehw,apresg
c
c     "Eh reaction" parameters.
c
      real(8) xhfe,xlke,xvfe
c
      real(8), dimension(:), allocatable :: dhfe,dvfe
      real(8), dimension(:,:), allocatable :: axlke,axhfe,axvfe
      real(8), dimension(:,:,:), allocatable :: adhfe,advfe
c
c     Regular reaction parameters.
c
      real(8), dimension(:), allocatable :: xhfs,xhfsd,xlks,xlksd,
     $ xvfs,xvfsd
c
      real(8), dimension(:,:), allocatable :: dhfs,dhfsd,dvfs,dvfsd
c
      real(8), dimension(:,:,:), allocatable :: axhfs,axhfsd,axhfsx,
     $ axlks,axlksd,axlksx,axvfs,axvfsd,axvfsx
c
      real(8), dimension(:,:,:,:), allocatable :: adhfs,adhfsd,adhfsx,
     $ advfs,advfsd,advfsx
c
      end module mod6xf
