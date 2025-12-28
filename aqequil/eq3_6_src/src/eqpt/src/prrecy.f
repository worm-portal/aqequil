      subroutine prrecy(cdrs,nbtmx1,nbtmx2,nbt,ns,nfile,uspec)
c
c     This subroutine writes the reaction associated with the ns-th
c     species on the file whose unit number is nfile. This subroutine
c     is virtually identical in function to EQLIB/prreac.f, but
c     deals with data that is structured somewhat differently.
c
c     This subroutine is called by:
c
c       EQPT/rxnchk.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       cdrs   = array of reaction coefficients
c       ns     = species whose reaction is to be printed
c       nfile  = unit number of the file to write on
c       uspec  = array of species names
c
c     Principal input:
c
c       None
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbtmx1,nbtmx2
c
      integer nbt,nfile,ns
c
      character*24 uspec(nbtmx1)
c
      real*8 cdrs(nbtmx2,nbtmx1)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,nj,nb,nbt1
c
      integer ilnobl
c
      real*8 cx
c
c-----------------------------------------------------------------------
c
      nbt1 = nbt + 1
      cx = -cdrs(nbt1,ns)
      j2 = ilnobl(uspec(ns))
      write (nfile,1000) cx,uspec(ns)(1:j2)
 1000 format(/21x,f10.5,2x,a)
c
      do nb = 1,nbt
        cx = -cdrs(nb,ns)
        if (cx .gt. 0.) then
          j2 = ilnobl(uspec(nb))
          write (nfile,1010) cx,uspec(nb)(1:j2)
 1010     format(18x,' + ',f10.5,2x,a)
        endif
      enddo
c
      write (nfile,1020)
 1020 format(/31x,'='/)
c
      nj = 0
      do nb = 1,nbt
        cx = cdrs(nb,ns)
        if (cx .gt. 0.) then
          j2 = ilnobl(uspec(nb))
          if (nj .le. 0) then
            write (nfile,1030) cx,uspec(nb)(1:j2)
 1030       format(21x,f10.5,2x,a)
            nj = 1
          else
            write (nfile,1010) cx,uspec(nb)(1:j2)
          endif
        endif
      enddo
c
      end
