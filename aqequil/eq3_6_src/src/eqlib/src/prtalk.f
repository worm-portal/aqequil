      subroutine prtalk(alki,alk1,alk2,mrmlra,noutpt,ntf1t,ntf2t,
     $ qrho,rho,tempc,wfh2o)
c
c     This subroutine writes a table of computed alkalinity parameters.
c
c     This subroutine is called by:
c
c       EQ3NR/scripx.f
c       EQ6/scripz.f
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
      integer noutpt
c
      integer ntf1t,ntf2t
c
      real(8) alki,alk1,alk2,mrmlra,rho,tempc,wfh2o
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2
c
      integer ilnobl
c
      logical qrho
c
      character(len=16) ux16
c
      real(8) alkicv,alkicw,alkihv,alkihw,alkipl,alk1cv,alk1cw,alk1hv,
     $ alk1hw,alk1pl,alk2cv,alk2cw,alk2hv,alk2hw,alk2pl
c
c-----------------------------------------------------------------------
c
c     Check the temperature. Alkalinity of any kind is only defined
c     in the temperature range 0-50 C.
c
      if (tempc.lt.0. .or. tempc.gt.50.) then
        write (ux16,'(f12.3)') tempc
        call lejust(ux16)
        j2 = ilnobl(ux16)
        write (noutpt,1000) ux16(1:j2)
 1000   format(//6x,'Alkalinity is not defined at ',a,' C.',
     $  /6x,'It is only defined in the temperature range 0-50 C.',/)
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write a table of the HCO3-CO3-OH total alkalinity.
c
      write (noutpt,1010)
 1010 format(/11x,'--- HCO3-CO3-OH Total Alkalinity ---',/)
c
c     Test for the required coefficients.
c
      if (ntf1t .le. 0) then
        write (noutpt,1020)
 1020   format(6x,'Coefficients to compute this alkalinity',
     $  ' are not available.',/)
        go to 110
      endif
c
      alk1cw = 50000*wfh2o*alk1
      alk1hw = 1.2192*alk1cw
c
      if (qrho) then
c
c       Have density data, include volumetric measures.
c
        alk1pl = alk1*mrmlra
        alk1cv = alk1cw*rho
        alk1hv = 1.2192*alk1cv
c
        write (noutpt,1030) alk1,alk1pl,alk1cw,alk1hw,alk1cv,alk1hv
 1030   format(16x,1pg12.5,' eq/kg.H2O',/16x,1pg12.5,' eq/L',
     $  /16x,1pg12.5,' mg/kg.sol CaCO3',/16x,1pg12.5,' mg/kg.sol HCO3-',
     $  /16x,1pg12.5,' mg/L CaCO3',/16x,1pg12.5,' mg/L HCO3-',/)
      else
c
c       Have no density data, do not include volumetric measures.
c
        write (noutpt,1040) alk1,alk1cw,alk1hw
 1040   format(16x,1pg12.5,' eq/kg.H2O',/16x,1pg12.5,' mg/kg.sol CaCO3',
     $  /16x,1pg12.5,' mg/kg.sol HCO3-',/)
      endif
c
  110 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write a table of the extended total alkalinity.
c
      write (noutpt,1050)
 1050 format(/11x,'--- Extended Total Alkalinity ---',/)
c
c     Test for the required coefficients.
c
      if (ntf2t .le. 0) then
        write (noutpt,1020)
        go to 120
      endif
c
      alk2cw = 50000*wfh2o*alk2
      alk2hw = 1.2192*alk2cw
c
      if (qrho) then
c
c       Have density data, include volumetric measures.
c
        alk2pl = alk2*mrmlra
        alk2cv = alk2cw*rho
        alk2hv = 1.2192*alk2cv
c
        write (noutpt,1030) alk2,alk2pl,alk2cw,alk2hw,alk2cv,alk2hv
      else
c
c       Have no density data, do not include volumetric measures.
c
        write (noutpt,1040) alk2,alk2cw,alk2hw
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  120 if (alki .gt. 0.) then
c
c       Write a table of the input alkalinity. This is relevant in
c       EQ3NR, but not in EQ6.
c
        write (noutpt,1060)
 1060   format(/11x,'--- Input Total Alkalinity ---',/)
c
        alkicw = 50000*wfh2o*alki
        alkihw = 1.2192*alkicw
c
        if (qrho) then
c
c         Have density data, include volumetric measures.
c
          alkipl = alki*mrmlra
          alkicv = alkicw*rho
          alkihv = 1.2192*alkicv
          write (noutpt,1030) alki,alkipl,alkicw,alkihw,alkicv,alkihv
        else
c
c         Have no density data, do not include volumetric measures.
c
          write (noutpt,1040) alki,alkicw,alkihw
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
