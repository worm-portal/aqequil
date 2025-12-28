      module mod6pt
c
c     This module contains data required to evaluate Pitzer's
c     equations. It excludes the variables and arrays that belong
c     to the 'a' set, which includes all the Pitzer data as read
c     from the data file. Included here are the "compressed"
c     equivalents.
c
c     Integer scalars.
c
      integer napt,nmut,nslt
c
cxx    integer ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,
cxx  $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta
c
c     Integer arrays. This are basically index arrays that allow
c     for efficient storage of the data.
c
      integer, dimension(:), allocatable :: nalpha
      integer, dimension(:,:), allocatable :: nmux,nmxi,nmxx,nslx,
     $ nsxi,nsxx
c
c     Coefficients for calculating empirical interaction coefficients
c     as functions of temperature and pressure.
c
      real(8), dimension(:,:), allocatable :: amu
      real(8), dimension(:,:,:), allocatable :: aslm
c
c     Empirical interaction coefficients and related functions at
c     temperature and pressure.
c
      real(8), dimension(:), allocatable :: pmu,pslm
      real(8), dimension(:,:), allocatable :: dpslm,gpit,palpha,pslamn
      real(8), dimension(:,:,:), allocatable :: dgpit
c
c     Theoretical higher-order electrostatic interaction coefficients
c     and related functions at temperature and pressure.
c
      real(8), dimension(:,:), allocatable :: elam,pelm
      real(8), dimension(:,:,:), allocatable :: delam,dpelm
c
      real(8), dimension(:), allocatable :: selm
      real(8), dimension(:,:), allocatable :: dselm
c
      end module mod6pt
