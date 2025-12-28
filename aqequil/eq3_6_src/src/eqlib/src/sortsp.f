      subroutine sortsp(iern1,iern2,istack,jcsort,jern1,jern2,
     $ jgext,jsitex,jetmax,jjsort,jssort,jstack,losp,lsort,ncmpr,
     $ nern1,nern2,netmax,noutpt,nphasx,npt,nptmax,nst,nstmax,nttyo)
c
c     This subroutine sorts the species in order of increasing mass.
c
c     This subroutine is called by:
c
c       EQLIB/ncmpex.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       jern1  = array marking the start of a site range for a
c                  generic ion exchange phase
c       jern2  = array marking the end of a site range for a
c                  generic ion exchange phase
c       iern1  = start of the phase range for generic ion exchange
c                  phases
c       iern2  = end of the phase range for generic ion exchange
c                  phases
c       jgext  = array giving the number of sites in generic
c                  generic ion exchange phases
c       jsitex = site index array, gives the site to which a
c                  species belongs
c       losp   = array of logarithms of numbers of moles of species
c       ncmpr  = species index array, gives the first and last
c                  species in a given phase
c       nern1  = start of the species range for generic ion exchange
c                  phases
c       nern2  = end of the species range for generic ion exchange
c                  phases
c       nphasx = phase index array, gives the phase to which a
c                  species belongs
c       nst    = number of species
c
c     Principal output:
c
c       jcsort = species index array, arranged in order of increasing
c                  mass, but arranged according to phase
c       jjsort = species index array, arranged in order of increasing
c                  mass, but arranged according to site and phase
c       jssort = species index array, arranged in order of increasing
c                  mass
c
c     Scratch arrays:
c
c       lsort, istack, and jstack
c
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations:
c
      integer jetmax,netmax,nptmax,nstmax
c
      integer noutpt,nttyo
c
      integer iern1,iern2,nern1,nern2
c
      integer istack(nstmax),jcsort(nstmax),jern1(jetmax,netmax),
     $ jern2(jetmax,netmax),jgext(netmax),jjsort(nstmax),
     $ jsitex(nstmax),jssort(nstmax),jstack(nstmax),ncmpr(2,nptmax),
     $ nphasx(nstmax)
c
      integer npt,nst
c
      real*8 losp(nstmax),lsort(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer je,ne,np,nrr1,nrr2,ns,nss,nsi
c
c-----------------------------------------------------------------------
c
c     Sort all species according to numbers of moles. There indices
c     in sorted order in the jssort array. Species belonging to
c     different phases are mixed together in this sort.
c
c     Calling sequence substitutions:
c       lsort for asort
c       losp for aval
c       jssort for jsort
c       nstmax for nmax
c       nst for nval
c
c     Caution: the jssort array from the last call is recycled as a
c     good starting point. Set jssort(1) to 0 to make a sort starting
c     from scratch.
c
      call qsortw(lsort,losp,istack,jssort,jstack,nstmax,noutpt,
     $ nttyo,nst)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the jcsort array using the jssort array already computed.
c     In this array, the species belonging to a given phase are grouped
c     together within the specified phase limits. The scratch array
c     istack is used here to represent the next available species
c     position within the range for any given phase. Note that istack
c     is dimensioned 1:nstmax. Here it must be of length at least
c     nptmax, which is necessarily no greater than nstmax.
c
      do np = 1,npt
        istack(np) = ncmpr(1,np)
      enddo
c
      do nss = 1,nst
        ns = jssort(nss)
        np = nphasx(ns)
        nsi = istack(np)
        jcsort(nsi) = ns
        istack(np) = istack(np) + 1
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute the jjsort array using the jcsort array already computed.
c     In this array, the species belonging to a given phase are grouped
c     together within the specified site limits. Presently only generic
c     ion exchange phases have sites in the present context. In the
c     future, other types of phases with distinct sites may be added
c     to the list. Phases without distinct sites are treated as having
c     "one" site.
c
c     The scratch array istack is used here to represent the next
c     available species position within the range for any given site
c     of the phase currently being processed. any given phase. Note
c     that istack is dimensioned 1:nstmax. Here it must have a length
c     of at least the maximum number of species in any possible phase,
c     which is necessarily no greater than nstmax.
c
      do ns = 1,nst
        jjsort(ns) = jcsort(ns)
      enddo
c
      do ns = nern1,nern2
        jjsort(ns) = 0
      enddo
c
      do np = iern1,iern2
        ne = np - iern1 + 1
        do je = 1,jgext(ne)
          istack(je) = jern1(je,ne)
        enddo
c
        nrr1 = ncmpr(1,np)
        nrr2 = ncmpr(2,np)
        do nss = nrr1,nrr2
          ns = jcsort(nss)
          je = jsitex(ns)
          nsi = istack(je)
          jjsort(nsi) = ns
          istack(je) = istack(je) + 1
        enddo
      enddo
c
      end
