      subroutine cwrpjt(noutpt)
c
c     This subroutine calculates and writes tables of the higher-order
c     electrostatic term functions J(x) and J'(x). These functions may
c     be computed from two different approximate formulations. The
c     first is one given by Pitzer (1975, eq. 47). The second is given
c     by Harvie (1981). Almost all modern application of the Pitzer
c     equations for aqueous electrolyte thermodynamics uses the more
c     recent Harvie formulation. These formulations are all
c     approximate in nature.
c
c     It is important to note that the tabulation in Table II of
c     Pitzer (1975) is not for the formulation represented by
c     eq. 47. Therefore, the results given here for the eq. 47
c     approximation will not precisely match the results given in
c     that table. Pitzer (1975) does not provide such a table for
c     the eq. 47 formulation.
c
c                            References
c
c     Harvie, C. E. 1981. Theoretical Investigations in Geochemistry
c       and Atom Surface Scattering. Ph.D. dissertation, University
c       of California, San Diego (Available as #8203026 from
c       University Microfilms International, Ann Arbor, Michigan).
c
c     Pitzer, K.S. 1975. Thermodynamics of electrolytes. V. Effects
c       of higher-order electrostatic terms. Journal of Solution
c       Chemistry, v. 4, p. 249-265.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
c       EQ3NR/eq3nr.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       noutpt = unit number of the output file
c
c     Principal output:
c
c       None returned
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ktable,n
c
      logical qpit75
c
      real(8) dhj0,d2hj0,hj0,hj1,hj2,x
c
c-----------------------------------------------------------------------
c
      write (noutpt,1000)
 1000 format(/1x,"Tables of J(x) and J'(x) (used with Pitzer's",
     $ ' equations)')
c
c     Tables using the Pitzer (1975) table format.
c
      write (noutpt,1010)
 1010 format(//3x,'Tables in the format of Pitzer (1975, Table II)',
     $ /3x,'(Note: that table was for a different approximation than',
     $ /3x,'than Pitzer, 1975, eq. 47, so results for eq. 47 will',
     $ /3x,'not precisely match those in Table II)')
c
      ktable = 1
      qpit75 = .true.
  100 continue
c
      if (qpit75) then
        write (noutpt,1020)
 1020   format(//5x,'Pitzer (1975, eq. 47) formulation')
      else
        write (noutpt,1030)
 1030   format(//5x,'Harvie (1981, Appendix B) formulation')
      endif
c
      write (noutpt,1040)
 1040 format(/7x,"x           J(x)        J'(x)",/)
c
      x = 0.
      do n = 1,10
        x = x + 0.01
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1050) x,hj0,dhj0
 1050   format(1x,f9.3,3x,f13.8,3x,f7.4)
      enddo
c
c     x should now be 0.1
c
      do n = 1,5
        x = x + 0.02
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1050) x,hj0,dhj0
      enddo
c
c     x should now be 0.2
c
      do n = 1,10
        x = x + 0.04
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1050) x,hj0,dhj0
      enddo
c
c     x should now be 0.6
c
      do n = 1,7
        x = x + 0.2
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1050) x,hj0,dhj0
      enddo
c
c     x should now be 2.0
c
      do n = 1,8
        x = x + 1.0
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1050) x,hj0,dhj0
      enddo
c
c     x should now be 10.0
c
      x = x + 2.0
      if (.not.qpit75) then
        call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
      else
        call gpj0(dhj0,d2hj0,hj0,x)
      endif
      write (noutpt,1050) x,hj0,dhj0
c
c     x should now be 12.0
c
      do n = 1,7
        x = x + 4.0
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1050) x,hj0,dhj0
      enddo
c
c     x should now be 40.0
c
      do n = 1,6
        x = x + 10.0
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1050) x,hj0,dhj0
      enddo
c
c     x should now be 100.0
c
      x = x + 100.0
      if (.not.qpit75) then
        call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
      else
        call gpj0(dhj0,d2hj0,hj0,x)
      endif
      write (noutpt,1050) x,hj0,dhj0
c
c     x should now be 200.0
c
      do n = 1,4
        x = x + 200.0
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1050) x,hj0,dhj0
      enddo
c
c     x should now be 1000.0
c
      x = x + 1000.0
      if (.not.qpit75) then
        call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
      else
        call gpj0(dhj0,d2hj0,hj0,x)
      endif
      write (noutpt,1050) x,hj0,dhj0
c
c     x should now be 2000.0
c
      do n = 1,4
        x = x + 2000.0
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1050) x,hj0,dhj0
      enddo
c
c     x should now be 10000.0
c
      if (ktable .le. 1) then
        qpit75 = .false.
        ktable = ktable + 1
        go to 100
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Tables using the Harvie (1981) table format (extended).
c
      write (noutpt,1060)
 1060 format(//3x,'Tables in the format of Harvie (1981, Appendix B)',
     $ /3x,'(extended format)')
c
      ktable = 1
      qpit75 = .true.
  120 continue
c
      if (qpit75) then
        write (noutpt,1020)
      else
        write (noutpt,1030)
      endif
c
      write (noutpt,1040)
c
      x = 0.
      do n = 1,10
        x = x + 0.001
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1070) x,hj0,dhj0
 1070   format(1x,f9.3,3x,f12.7,4x,f10.7)
      enddo
c
c     x should now be 0.01
c
      do n = 1,9
        x = x + 0.01
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1070) x,hj0,dhj0
      enddo
c
c     x should now be 0.1
c
      do n = 1,9
        x = x + 0.1
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1070) x,hj0,dhj0
      enddo
c
c     x should now be 1.0
c
      do n = 1,9
        x = x + 1.0
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1070) x,hj0,dhj0
      enddo
c
c     x should now be 10.0
c
      do n = 1,9
        x = x + 10.0
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1070) x,hj0,dhj0
      enddo
c
c     x should now be 100.0
c
      do n = 1,9
        x = x + 100.0
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1070) x,hj0,dhj0
      enddo
c
c     x should now be 1000.0
c
      do n = 1,9
        x = x + 1000.0
        if (.not.qpit75) then
          call ghj0(dhj0,d2hj0,hj0,hj1,hj2,x)
        else
          call gpj0(dhj0,d2hj0,hj0,x)
        endif
        write (noutpt,1070) x,hj0,dhj0
      enddo
c
c     x should now be 10000.0
c
      if (ktable .le. 1) then
        qpit75 = .false.
        ktable = ktable + 1
        go to 120
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,2000)
 2000 format(/' - - - - - - - - - - - - - - - - - - - - - - - - ',
     $ '- - - - - - ',/)
      end
