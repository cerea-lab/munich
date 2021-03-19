C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Kathleen Fahey
C     
C     This file is part of the Variable Size Resolved Model (VSRM),
C     based on the VSRM model of Carnegie Melon University.  It is a
C     component of the air quality modeling system Polyphemus.
C    
C     Polyphemus is developed in the INRIA - ENPC joint project-team
C     CLIME and in the ENPC - EDF R&D joint laboratory CEREA.
C    
C     Polyphemus is free software; you can redistribute it and/or modify
C     it under the terms of the GNU General Public License as published
C     by the Free Software Foundation; either version 2 of the License,
C     or (at your option) any later version.
C     
C     Polyphemus is distributed in the hope that it will be useful, but
C     WITHOUT ANY WARRANTY; without even the implied warranty of
C     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
C     General Public License for more details.
C     
C     For more information, visit the Polyphemus web site:
C     http://cerea.enpc.fr/polyphemus/
C-----------------------------------------------------------------------

      subroutine electro(x,con,spres,cmet,akeq,akhen,wv,temp,f,frel)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This routine computes the electroneutrality balance for the 
C     aqueous-phase model.
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     X    : H+ concentration.
C     CON  : concentration vector (28 species).
C     SPRES: gas-phase concentrations.
C     CMET : metal concentrations.
C     AKEQ : equilibrium rate constants.
C     AKHEN: Henry's rates.
C     WV   : water vapour.
C     TEMP : temperature.
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     -- OUTPUT VARIABLES
C     
C     F : evaluation of the electroneutrality relation (anion-cation).
C     FL: relative value (ABS(F)/(anion+cation).
C------------------------------------------------------------------------
C     
C     -- REMARKS
C     
C     The aqueous-phase composition may be found in values.f.
C     
C------------------------------------------------------------------------
C     
C     -- MODIFICATIONS
C     
C     1) Add relative value FL.
C     2) Change order or arguments (F).
C     3) Remove UU.
C     4) Optimize (coefloc and chlorine).
C------------------------------------------------------------------------
C     
C     -- AUTHOR(S)
C     
C     Kathleen Fahey, CEREA, , on the basis of the VSRM model 
C     (Carneggie Mellon University).
C     2005/10/3, cleaning and update, Bruno Sportisse, CEREA.
C     
C------------------------------------------------------------------------

      IMPLICIT NONE

      double precision con(28),spres(21),cmet(4),akeq(17),akhen(21)
      double precision temp,x,wv
      double precision f,frel

      double precision cc(46)
      double precision bparam,cparam,diak,cl,hcl
      double precision coefloc,dfac,dform
      double precision fanion,fcation
      double precision formald,hno2


C     1)Compute ions in aqueous-phase
C     
C     HNO2(g), CO2(g) and HCOOH(g) are 
C     computed from equilibrium
C     -----------------------------------

      cc(2)=(akeq(1)*con(1)*x)/(x*x+akeq(1)*x+akeq(1)*akeq(2)) !HSO3-
      cc(3)=(akeq(1)*akeq(2)*con(1))/
     &     (x*x+akeq(1)*x+akeq(1)*akeq(2)) !SO3--
      cc(5)=(akeq(3)*con(2)*x)/(x*x+akeq(3)*x+akeq(3)*akeq(4)) !HSO4-
      cc(6)=(akeq(3)*akeq(4)*con(2))/
     &     (x*x+akeq(3)*x+akeq(3)*akeq(4)) !SO4--

      coefloc = 8.314d-2*temp*wv
      
      dfac  = coefloc*akhen(3)*(1.d0+akeq(7)/x)
      hno2  = spres(3)/(1.d0+dfac) !New HNO2(g) in ppm
      cc(8) = akhen(3)*1.d-6*(akeq(7)/x)*hno2 !NO2-
      cc(10)= (akeq(6)*con(4))/(x+akeq(6)) !NO3-

      cc(12)= akeq(8)*akhen(5)*spres(5)*1.d-6/x !HCO3-
      cc(13)= akeq(9)*cc(12)/x  !CO3--
      
      cc(15)=(akeq(5)*con(6))/(x+akeq(5)) !HO2-

      dform = coefloc*akhen(8)*(1.d0+akeq(13)/x)
      formald=spres(8)/(1.d0+dform) !New HCOOH(g)
      cc(19)=akhen(8)*1.d-6*(akeq(13)/x)*formald !HCOO-

      cc(30)=(akeq(15)*con(17))/(x+akeq(15)) !O2-
      cc(38)=con(23)            !ClOH-
      cc(39)=con(24)            !SO4-
      cc(40)=con(25)            !SO5-
      cc(41)=con(26)            !HSO5-
      cc(42)=(x*con(27))/(x+akeq(17)) !HOCH2SO3-
      cc(43)=(akeq(17)*con(27))/(x+akeq(17)) !-OCH2SO3-
      cc(44)=con(28)            !CO3-
      cc(45)=akeq(11)/x         !OH-

C     Chlorine

      bparam = akeq(16)+con(15)-con(22)
      cparam = akeq(16)*con(22)
      diak   = dmax1(bparam*bparam+4.d0*cparam,1.0d-20)

      cl    = (-bparam+diak**0.5d0)/2.d0
      hcl   = (x*(con(15)-con(22)+cl))/(x+akeq(14))

      cc(27)=(akeq(14)*hcl)/x   !Cl-
      cc(36)=cl*cc(27)/akeq(16) !Cl2-

C     
      cc(33)=(akeq(10)*x*con(19))/
     &     (akeq(11)+akeq(10)*x) !NH4+
      cc(46)=x                  !H+

C     2) Compute electroneutrality balance
C     -----------------------------------

      fanion  = cc(2)+2.d0*cc(3)+cc(5)+2.d0*cc(6)+cc(8)+cc(10)
     &     + cc(12)+2.d0*cc(13)+cc(15)+cc(19)+cc(27)+cc(30)
     &     + cc(36)+cc(38)+cc(39)+cc(40)+cc(41)+cc(42)
     &     + 2.d0*cc(43)+cc(44)+cc(45)
      fcation = cc(33)+cc(46)
     &     + 3.d0*cmet(1)+2.d0*cmet(2)+cmet(3)+2.d0*cmet(4)

      f= fanion -fcation

      if (fanion+fcation.gt.0.d0) then
         frel=abs(f)/(fanion+fcation)
      else
         frel=-1.d0
      endif

      return
      end




