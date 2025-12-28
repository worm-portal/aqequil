      subroutine reacts(cbsr,csts,delxi,drer0,iern1,ietmax,iktmax,
     $ iodb,jcode,jetmax,jgext,jreac,modr,modr0,morr,morr0,mrgers,
     $ mtb,mtb0,nbaspd,nbt,nbtmax,nbt1mx,ncmpr,nern1,nern2,nertmx,
     $ netmax,ngext,nodbmx,nord,noutpt,nptmax,nrct,nrctmx,nrd1mx,
     $ nrndex,nsrtmx,nstmax,nsts,nstsmx,nstsr,nttyo,nxridx,nxrtmx,
     $ rrelr0,rxbar,ureac,xirct,xirct0)
c
c     This subroutine computes the destroyed and current masses of the
c     reactants and the current mass balance totals for the equilibrium
c     system.
c
c        morr   = moles of reactant remaining
c        modr   = moles of irreversibly destroyed reactant
c        jreac  = reactant control switch
c           = -1  The reactant has saturated but continues to
c                 be available for irreversible reaction
c           =  0  The reactant is currently reacting
c           =  1  The reactant has been exhausted
c           =  2  The reactant has saturated and any remaining mass
c                 has been transferred to the equilibrium system
c                 (This occurs only if iopt(1) = 0)
c
c     This subroutine is called by:
c
c       EQ6/eqshel.f
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
c     Calling sequence variable declarations.
c
      integer ietmax,iktmax,jetmax,nbtmax,nbt1mx,nertmx,netmax,nodbmx,
     $ nptmax,nrctmx,nrd1mx,nsrtmx,nstmax,nstsmx,nxrtmx
c
      integer iodb(nodbmx),jcode(nrctmx),jgext(netmax),jreac(nrctmx),
     $ nbaspd(nbtmax),ncmpr(2,nptmax),ngext(jetmax,netmax),
     $ nrndex(nrctmx),nsts(nstsmx),nstsr(2,nstmax),nxridx(nrctmx)
c
      integer iern1,nbt,nern1,nern2,nord,noutpt,nrct,nttyo
c
      character*24 ureac(nrctmx)
c
      real*8 cbsr(nbt1mx,nsrtmx),csts(nstsmx),drer0(nrd1mx,nrctmx),
     $ modr(nrctmx),modr0(nrctmx),morr(nrctmx),morr0(nrctmx),
     $ mrgers(ietmax,jetmax,nertmx),mtb(nbtmax),mtb0(nbtmax),
     $ rrelr0(nrctmx),rxbar(iktmax,nxrtmx),xirct(nrctmx),
     $ xirct0(nrctmx)
c
      real*8 delxi
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ie,ik,je,j2,n,nb,ne,ner,np,nsr,nrc,nrn1,nrn2,nr1,nr2,
     $ ns,nxr
c
      integer ilnobl
c
      real*8 dlxrct,dmodr,dx
c
c-----------------------------------------------------------------------
c
c     Reset morr, modr, xirct, and mtb to values at the previous step.
c
      do nb = 1,nbt
        mtb(nb) = mtb0(nb)
      enddo
c
      do nrc = 1,nrct
        morr(nrc) = morr0(nrc)
        modr(nrc) = modr0(nrc)
        xirct(nrc) = xirct0(nrc)
      enddo
c
      if (delxi .le. 0.) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Increment morr, modr, xirct, and mtb to correspond to the
c     advance in reaction progress.
c
      if (iodb(1) .ge. 3) write (noutpt,1000)
 1000 format(' --- Advancement of Irreversible Reactions ---',
     $ /7x,'Reaction',14x,'dlxrct',7x,'dmodr')
c
      do nrc = 1,nrct
        if (jreac(nrc).eq.0 .or. jreac(nrc).eq.-1) then
c
          call integr(delxi,dlxrct,drer0,nord,nrc,nrctmx,
     $    nrd1mx,rrelr0)
c
cXX       It is assumed that each reactant has a reaction coefficient
cXX       of -1. Below it is implied that dlxrct is multiplied by -(-1).
c
          dmodr = dlxrct
          if (dmodr .gt. morr(nrc)) then
c
c           The amount destroyed in the current step slightly exceeds
c           the amount remaining.
            dmodr = morr(nrc)
c
cXX         It is assumed that each reactant has a reaction coefficient
cXX         of -1. Below it is implied that dlxrct is multiplied
cXX         by -(-1).
c
            dlxrct = dmodr
          endif
c
          if (iodb(1) .ge. 3) write (noutpt,1010) ureac(nrc),dlxrct,
     $    dmodr
 1010     format(3x,a24,2(3x,1pe12.5))
          if (dlxrct .ne. 0.) then
            xirct(nrc) = xirct(nrc) + dlxrct
            modr(nrc) = modr(nrc) + dmodr
            morr(nrc) = morr(nrc) - dmodr
c
            if (jcode(nrc) .eq. 0) then
c
c             Pure mineral.
c
              np = nrndex(nrc)
              ns = ncmpr(1,np)
              nr1 = nstsr(1,ns)
              nr2 = nstsr(2,ns)
              do n = nr1,nr2
                nb = nsts(n)
                mtb(nb) = mtb(nb) + dmodr*csts(n)
              enddo
c
            elseif (jcode(nrc) .eq. 1) then
c
c             Solid solution.
c
              nxr = nxridx(nrc)
              np = nrndex(nrc)
              nrn1 = ncmpr(1,np)
              nrn2 = ncmpr(2,np)
              ik = 0
              do ns = nrn1,nrn2
                ik = ik + 1
                nr1 = nstsr(1,ns)
                nr2 = nstsr(2,ns)
                dx = dmodr*rxbar(ik,nxr)
                do n = nr1,nr2
                  nb = nsts(n)
                  mtb(nb) = mtb(nb) + dx*csts(n)
                enddo
              enddo
c
            elseif (jcode(nrc) .eq. 2) then
c
c             Special reactant.
c
              nsr = nrndex(nrc)
              do nb = 1,nbt
                mtb(nb) = mtb(nb) + dmodr*cbsr(nb,nsr)
              enddo
c
            elseif (jcode(nrc).eq.3 .or. jcode(nrc).eq.4) then
c
c             Aqueous species or gas.
c
              ns = nrndex(nrc)
              nr1 = nstsr(1,ns)
              nr2 = nstsr(2,ns)
              do n = nr1,nr2
                nb = nsts(n)
                mtb(nb) = mtb(nb) + dmodr*csts(n)
              enddo
c
            elseif (jcode(nrc) .eq. 5) then
c
c             Generic ion exchanger.
c
              ner = nxridx(nrc)
              np = nrndex(nrc)
              ne = np - iern1 + 1
              ns = ncmpr(1,np) - 1
              do je = 1,jgext(ne)
                do ie = 1,ngext(je,ne)
                  dx = dmodr*mrgers(ie,je,ner)
                  ns = ns + 1
                  nr1 = nstsr(1,ns)
                  nr2 = nstsr(2,ns)
                  if (je .le. 1) then
c
c                   Have the first site. Increment the mass balance
c                   totals in a straightforward manner.
c
                    do n = nr1,nr2
                      nb = nsts(n)
                      mtb(nb) = mtb(nb) + dx*csts(n)
                    enddo
                  else
c
c                   Have a site beyond the first. Increment the mass
c                   balance totals in the usual manner, except do not
c                   increment here the total for the exchanger
c                   substrate. The complete increment for the substrate
c                   is obtained by considering only one site.
c
                    do n = nr1,nr2
                      nb = nsts(n)
                      ns = nbaspd(nb)
                      if (ns.lt.nern1 .or. ns.gt.nern2) then
                        mtb(nb) = mtb(nb) + dx*csts(n)
                      endif
                    enddo
                  endif
                enddo
              enddo
            else
              j2 = ilnobl(ureac(nrc))
              write (noutpt,64) jcode(nrc),ureac(nrc)(1:j2)
              write (nttyo,64) jcode(nrc),ureac(nrc)(1:j2)
   64         format(/' * Error - (EQ6/reacts) Programming error',
     $        ' trap: Have unknown',/7x,'reactant type code (jcode)',
     $        ' value of ',i5,' for reactant ',/7x,a,'.')
              stop
            endif
c
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
