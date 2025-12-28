      subroutine prtntt(nat,nata,natmax,nbt,nbta,nbtmax,nct,ncta,
     $ nctmax,ngt,ngta,ngtmax,nlt,nlta,nltmax,nmt,nmta,nmtmax,noutpt,
     $ npt,npta,nptmax,nst,nsta,nstmax,nxt,nxta,nxtmax)
c
c     This subroutine prints a table of statistics for species, phases,
c     and groups thereof showing for each entity the number on the data
c     file, the number the software is dimensioned for, and the number
c     appearing in the current problem.
c
c     This subroutine is called by:
c
c       EQ3NR/echox.f
c       EQ6/echoz.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
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
      integer natmax,nbtmax,nctmax,ngtmax,nltmax,nmtmax,nptmax,nstmax,
     $ nxtmax
c
      integer noutpt
c
      integer nat,nata,nbt,nbta,nct,ncta,ngt,ngta,nlt,nlta,nmt,nmta,
     $ npt,npta,nst,nsta,nxt,nxta
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      character*24 ux
c
c-----------------------------------------------------------------------
c
      write (noutpt,1000)
 1000 format(/7x,'--- Numbers of Phases, Species, and Groups Thereof',
     $ '---',//7x,'Entity',15x,'Date Base',4x,'Dimension',3x,
     $ 'Current Problem',/)
c
      ux = 'Chemical Elements       '
      write (noutpt,1010) ux,ncta,nctmax,nct
 1010 format(3x,a24,3x,i5,7x,i5,7x,i5)
c
      ux = 'Basis Species'
      write (noutpt,1010) ux,nbta,nbtmax,nbt
c
      ux = 'Phases'
      write (noutpt,1010) ux,npta,nptmax,npt
c
      ux = 'Species'
      write (noutpt,1010) ux,nsta,nstmax,nst
c
      ux = 'Aqueous Species'
      write (noutpt,1010) ux,nata,natmax,nat
c
      ux = 'Pure Minerals'
      write (noutpt,1010) ux,nmta,nmtmax,nmt
c
      ux = 'Pure Liquids'
      write (noutpt,1010) ux,nlta,nltmax,nlt
c
      ux = 'Gas Species'
      write (noutpt,1010) ux,ngta,ngtmax,ngt
c
      ux = 'Solid Soutions'
      write (noutpt,1020) ux,nxta,nxtmax,nxt
 1020 format(3x,a24,3x,i5,7x,i5,7x,i5,/)
c
      end
