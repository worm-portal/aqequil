      subroutine echgex(axlks,cdrs,cgexj,iern1,iern2,jern1,jern2,
     $ jetmax,jgext,jpflag,jsflag,narxmx,narxt,ndrs,ndrsmx,ndrsr,
     $ netmax,noutpt,nptmax,ntprmx,ntprt,nstmax,press,tempc,ugexj,
     $ ugexmo,uphase,uspec,xlks)
c
c     This subroutine echoes a table for the generic ion exchangers,
c     describing the setup of species, reactions, and corresponding
c     thermodynamic data.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       noutpt = unit number of the output file
c
c     Principal output:
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
      integer jetmax,narxmx,ndrsmx,netmax,nptmax,ntprmx,nstmax
c
      integer noutpt
c
      integer jern1(jetmax,netmax),jern2(jetmax,netmax),jgext(netmax),
     $ jpflag(nptmax),jsflag(nstmax),narxt(ntprmx),ndrs(ndrsmx),
     $ ndrsr(2,nstmax)
c
      integer iern1,iern2,ntprt
c
      integer ilnobl
c
      character*48 uspec(nstmax)
      character*24 ugexmo(netmax),uphase(nptmax)
      character*8 ugexj(jetmax,netmax)
c
      real*8 axlks(narxmx,ntprmx,nstmax),cdrs(ndrsmx),
     $ cgexj(jetmax,netmax),xlks(nstmax)
c
      real*8 press,tempc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j,je,j3,j4,j5,ne,np,nr1,nr2,ns,nt
c
c-----------------------------------------------------------------------
c
      if (iern2 .lt. iern1) go to 999
c
      write (noutpt,1000)
 1000 format(/11x,' --- Generic Ion Exchange Phase Setup ---',/)
      write (noutpt,1010) tempc,press
 1010 format(6x,'Log K data, etc., are for ',f10.3,'C and ',g12.5,
     $ ' bars',/)
c
      do np = iern1,iern2
        ne = np - iern1 + 1
        j5 = ilnobl(ugexmo(ne))
        j4 = ilnobl(uphase(np))
        if (jpflag(np) .le. 0) then
          write (noutpt,1020) uphase(np)(1:j4)
 1020     format(/16x,'--- ',a,' ---')
        else
          write (noutpt,1030) uphase(np)(1:j4)
 1030     format(/16x,'--- ',a,' (Suppressed) ---')
        endif
        write (noutpt,1040) ugexmo(ne)(1:j5)
 1040   format(/11x,'Model type= ',a)
        do je = 1,jgext(ne)
          j3 = ilnobl(ugexj(je,ne))
          write (noutpt,1050) ugexj(je,ne)(1:j3)
 1050     format(/11x,'--- ',a,' ---',/)
          write (noutpt,1060) cgexj(je,ne)
 1060     format(6x,'Site stoichiometric factor N= ',1pg12.5,/)
          do ns = jern1(je,ne),jern2(je,ne)
c
c           Calling sequence substitutions:
c             noutpt for nf
c
            call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns,nstmax,uspec)
c
            nr1 = ndrsr(1,ns)
            nr2 = ndrsr(2,ns)
            nt = nr2 - nr1 + 1
            if (nt .gt. 1) then
              if (jsflag(ns) .le. 0) then
                write (noutpt,1070) xlks(ns)
 1070           format(/5x,'Log K= ',f12.4)
                do j = 1,ntprt
                  write (noutpt,1100) j
 1100             format(/7x,'Coefficients for temperature range ',i2,
     $            ':')
                  write (noutpt,1110) (axlks(i,j,ns), i = 1,narxt(j))
 1110             format( (4x,5(2x,g13.6)) )
                enddo
                write (noutpt,1150)
 1150           format(1x)
              else
                write (noutpt,1170)
 1170           format(2x,'Suppressed',/)
              endif
            endif
          enddo
        enddo
      enddo
c
  999 continue
      end
