C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Kathleen Fahey and Bruno Sportisse
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

      subroutine react(c,cmet,con,akre,rr,arytm)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This routine computes the reaction rates for the 
C     aqueous-phase model.
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     C     : concentration vector for aqueous model ([...]).
C     CMET  : metal concentrations    ([...]).
C     CON   : gas + ions    ([...]).
C     AKRE  : kinetic rates.
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     -- OUTPUT VARIABLES
C     
C     RR    : reaction rates.
C     ARYTM : coefficient for the S(IV)+N(III) reaction (depends on pH).
C     
C------------------------------------------------------------------------
C     
C     -- REMARKS
C     
C------------------------------------------------------------------------
C     
C     -- MODIFICATIONS
C     
C     1) Clean the IF statements with ELSEIF
C     2) Remove GOTO and replace by IF statements.
C     3) Optimize RR(109).
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

      include 'aerpar.inc'
      include 'droppar.inc'

      double precision c(46),cmet(4),con(28)
      double precision akre(120),rr(120)
      double precision arytm,ph,r1,r2
      double precision r3,r4,r5,sn
      
      double precision hmsa

c     Due to the fact that a transfer of the HMSA produced in the 
c     aqueous phase to sulfate can skew the sulfate to an excessively 
c     high level, one should not turn on reaction 106 for the production
c     of HMSA unless HMSA is treated as a regular species in the model.

      hmsa = 0.d0

      rr(1)=akre(1)*c(14)*photo
      rr(2)=akre(2)*c(22)*photo
      rr(3)=akre(3)*c(28)*c(29)
      rr(4)=akre(4)*c(28)*c(30)
      rr(5)=akre(5)*c(28)*c(14)
      rr(6)=akre(6)*c(29)*c(29)
      rr(7)=akre(7)*c(29)*c(30)
      rr(8)=akre(8)*c(30)*c(30)
      rr(9)=akre(9)*c(29)*c(14)
      rr(10)=akre(10)*c(30)*c(14)
      rr(11)=akre(11)*c(28)*c(22)
      rr(12)=akre(12)*c(29)*c(22)
      rr(13)=akre(13)*c(30)*c(22)
      rr(14)=akre(14)*c(45)*c(22)
      rr(15)=akre(15)*c(15)*c(22)
      if (c(22) .le. 0.0d0) c(22)=1.0d-30
      rr(16)=akre(16)*c(14)*(c(22)**0.5d0)

      rr(17)=akre(17)*c(12)*c(28)
      rr(18)=akre(18)*c(12)*c(30)
      rr(19)=akre(19)*c(44)*c(30)
      rr(20)=akre(20)*c(44)*c(14)
      rr(21)=akre(21)*c(27)*c(28)*chlorine
      rr(22)=akre(22)*c(38)*chlorine
      rr(23)=akre(23)*c(46)*c(38)*chlorine
      rr(24)=akre(24)*c(37)*chlorine
      rr(25)=akre(25)*c(29)*c(36)*chlorine
      rr(26)=akre(26)*c(30)*c(36)*chlorine
      rr(27)=akre(27)*c(29)*c(37)*chlorine
      rr(28)=akre(28)*c(14)*c(36)*chlorine
      rr(29)=akre(29)*c(37)*c(14)*chlorine
      rr(30)=akre(30)*c(45)*c(36)*chlorine

      rr(31)=akre(31)*c(20)*c(21)
      rr(32)=akre(32)*c(21)*c(21)
      rr(33)=akre(33)*c(20)*c(28)
      rr(34)=akre(34)*c(21)*c(28)
      rr(35)=akre(35)*c(7)*photo
      rr(36)=akre(36)*c(8)*photo
      rr(37)=akre(37)*c(7)*c(28)
      rr(38)=akre(38)*c(8)*c(28)
      rr(39)=akre(39)*c(46)*c(14)*c(7)
      rr(40)=akre(40)*c(8)*c(22)
      rr(41)=akre(41)*c(8)*c(44)
      rr(42)=akre(42)*c(8)*c(36)*chlorine
      rr(43)=akre(43)*c(8)*c(31)
      rr(44)=akre(44)*c(10)*photo
      rr(45)=akre(45)*c(31)*photo
      rr(46)=akre(46)*c(31)*c(29)
      rr(47)=akre(47)*c(31)*c(30)
      rr(48)=akre(48)*c(31)*c(14)
      rr(49)=akre(49)*c(31)*c(27)*chlorine
      rr(50)=akre(50)*c(17)*c(28)
      rr(51)=akre(51)*c(17)*c(22)
      rr(52)=akre(52)*c(18)*c(28)
      rr(53)=akre(53)*c(18)*c(14)
      rr(54)=akre(54)*c(18)*c(31)
      rr(55)=akre(55)*c(18)*c(22)
      rr(56)=akre(56)*c(18)*c(36)*chlorine
      rr(57)=akre(57)*c(19)*c(28)
      rr(58)=akre(58)*c(19)*c(22)
      rr(59)=akre(59)*c(19)*c(31)
      rr(60)=akre(60)*c(19)*c(44)
      rr(61)=akre(61)*c(19)*c(36)*chlorine
      rr(62)=akre(62)*c(23)
      rr(63)=akre(63)*c(34)*c(29)
      rr(64)=akre(64)*c(34)*c(30)
      rr(65)=akre(65)*c(25)*photo
      rr(66)=akre(66)*c(25)*c(28)
      rr(67)=akre(67)*c(35)*c(28)
      rr(68)=akre(68)*c(35)*c(44)
      rr(69)=akre(69)*c(35)*c(36)*chlorine
      rr(70)=akre(70)*c(25)*c(28)
      rr(71)=akre(71)*c(35)*c(31)
      rr(72)=(akre(72)*c(1)+akre(73)*c(2)+akre(74)*
     &     c(3))*c(22)
      rr(73)=(akre(75)*c(14)*c(1))/(1.0d0+16.0d0*c(46))


c     RATE EXPRESSIONS FOR THE METAL CATALYSED OXIDATION OF S(IV)

      ph=-dLOG10(c(46))

      if (kiron .eq. 1) then
c     ** PHENOMENOLOGICAL EXPRESSION BY MARTIN et al. (1991) **
         
         if (ph .le. 3.d0) then
            rr(74)=6.d0*cmet(1)*con(1)/c(46)
         elseif (ph .gt. 3.d0 .AND. ph .le. 4.5d0) then
            rr(74) = 1.d9*con(1)*cmet(1)*cmet(1)
         elseif (ph .gt. 4.5d0 .AND. ph .le. 6.5d0) then
            rr(74) = 1.0d-3*con(1)
         elseif (ph .gt. 6.5d0) then
            rr(74)=1.0d-4*con(1)
         endif

c     ** EXPRESSION BY MARTIN (1984) **
      elseif (kiron .eq. 2) then
         if ((c(46) .GE. 1.0d-4).AND.(con(1) .GE. 1.0d-5)) then
            r1=(akre(76)*cmet(2)*cmet(2))/c(46)
            r2=(akre(77)*cmet(1)*con(1)/c(46))
            if (cmet(2) .le. 0.0d0) cmet(2)=1.0d-30
            r3=r2*(1.0d0+(1.7d3*cmet(2)**1.5d0)/(6.3d-6+cmet(1)))
            rr(74)=r1+r3
         else
            if (cmet(1)*cmet(2) .LT. 1.0d-15) then
               sn=1.0d0
            else
               sn=3.0d0
            end if

            if ((c(46).GE.1.0d-4).AND.(con(1).LT.1.0d-5)) then
               rr(74)=sn*(akre(78)*cmet(2)*c(2)+akre(77)*cmet(1)*
     &              con(1)/c(46))

            elseif ((c(46).LT.1.0d-4).AND.(con(1).GE.1.0d-5)) then
               r4=akre(76)*cmet(2)*cmet(2)/c(46)
               r5=akre(79)*cmet(1)*con(1)*con(1)
               rr(74)=r4+r5
            else 
               rr(74)=akre(78)*cmet(2)*c(2)
            endif
         endif
      endif   
      rr(75)=akre(80)*c(3)*c(28)
      rr(76)=akre(81)*c(2)*c(28)
      rr(77)=akre(82)*c(40)*c(2)+akre(117)*c(40)*c(3)
      rr(78)=akre(83)*c(40)*c(30)
      rr(79)=akre(84)*c(40)*c(18)
      rr(80)=akre(85)*c(40)*c(19)
      rr(81)=akre(86)*c(40)*c(40)
      rr(82)=akre(87)*c(41)*c(2)*c(46)
      rr(83)=akre(88)*c(41)*c(28)
      rr(84)=akre(89)*c(41)*c(39)
      rr(85)=akre(90)*c(41)*c(8)
      rr(86)=akre(91)*c(41)*c(27)*chlorine
      rr(87)=akre(92)*c(39)*c(2)
      rr(88)=akre(93)*c(39)*c(3)
      rr(89)=akre(94)*c(39)*c(29)
      rr(90)=akre(95)*c(39)*c(30)
      rr(91)=akre(96)*c(39)*c(45)
      rr(92)=akre(97)*c(39)*c(14)
      rr(93)=akre(98)*c(39)*c(8)
      rr(94)=akre(99)*c(39)*c(12)
      rr(95)=akre(100)*c(39)*c(19)
      rr(96)=akre(101)*c(39)*c(27)*chlorine
      rr(97)=akre(102)*c(39)*c(18)
      rr(98)=akre(103)*c(23)*c(2)/c(46)
      rr(99)=akre(104)*c(2)*c(25)*c(46)
      rr(100)=(akre(105)*c(46)+akre(106))*c(2)*c(24)
      rr(101)=akre(107)*c(29)*c(3)+akre(118)*c(3)*c(30)
      rr(102)=akre(108)*c(39)*c(35)
      rr(103)=akre(109)*c(2)*c(31)
      rr(104)=akre(110)*con(1)*c(21)*nitrogen ! kmf (08/15/02)

      if (c(46) .GE. 1.0d-3) then
         rr(105)=akre(111)*con(3)*con(1)*c(46)**0.5d0
         arytm=1.0d0
      else
         rr(105)=akre(112)*c(8)*c(2)*c(46)
         arytm=0.0d0
      end if

      rr(106)=(akre(113)*c(16)*c(2)+akre(119)*c(16)*c(3))*hmsa 
      rr(107)=akre(114)*c(42)*c(45)
      rr(108)=akre(115)*c(42)*c(28)
      rr(109)=akre(116)*c(36)*chlorine*(c(2)+c(3))
      
      return
      end
