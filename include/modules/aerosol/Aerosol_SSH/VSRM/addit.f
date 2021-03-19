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


      subroutine addit(rr,arytm,rp,rl)
C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This routine computes the source terms for the 28 species of the
C     aqueous-phase model.
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     RR    : reaction rates. 
C     ARYTM : coefficient for SO2+HNO2 
C     (depends on pH in subroutine react)
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     -- OUTPUT VARIABLES
C     
C     RP : production flux.
C     RL : destruction flux.
C------------------------------------------------------------------------
C     
C     -- REMARKS
C     
C     See subroutine VALUES.f for the description of species (vector CON)
C     
C------------------------------------------------------------------------
C     
C     -- MODIFICATIONS
C     
C------------------------------------------------------------------------
C     
C     -- AUTHOR(S)
C     
C     2005/10/3, Bruno Sportisse, CEREA.
C     
C------------------------------------------------------------------------
c     
      IMPLICIT NONE

      double precision rr(120),rp(28),rl(28)
      double precision arytm
      
c     
c     ** S(IV) **
c     
      rp(1)=rr(107)
      rl(1)=+rr(72)+rr(73)+rr(74)+rr(98)+rr(101)+rr(105)*arytm
     &     +rr(76)+rr(77)+rr(82)+rr(87)+rr(99)+rr(100)
     &     +2.0d0*rr(103)
     &     +rr(104)+2.0d0*rr(105)*(1.0d0-arytm)+rr(106)+rr(109)
     &     +rr(75)+rr(88)
c     
c     ** S(VI) **
c     
      rp(2)= rr(72)+rr(73)+rr(74)+rr(98)+rr(101)+rr(105)*arytm
     &     +rr(85)
     &     +2.d0*rr(82)+rr(84)+rr(86)+rr(87)+rr(88)+rr(89)+rr(90)
     &     +rr(91)+rr(92)+rr(93)+rr(94)+rr(95)+rr(96)+rr(97)+rr(99)
     &     +rr(100)+rr(102)+rr(103)+rr(104)
      rl(2)=0.0d0
c     
c     ** N(III) **
c     
      rp(3)=2.0d0*rr(31)+rr(32)+rr(33)+2.0d0*rr(104)
      rl(3)=rr(35)+rr(37)+rr(39)+rr(105)*arytm+rr(36)+rr(38)+rr(40)
     &     +rr(41)+rr(42)+rr(43)+rr(85)+rr(93)+rr(105)*(1.0d0-arytm)
c     
c     ** N(V) **
c     
      rp(4)=rr(32)+rr(34)+rr(39)+rr(40)+rr(43)+rr(46)+rr(47)+rr(48)+
     &     rr(49)+rr(54)+rr(59)+rr(62)+rr(71)+rr(85)+rr(103)
      rl(4)=rr(44)
c     
c     ** CO2 **
c     
      rp(5)=rr(52)+rr(54)+rr(55)+rr(56)+rr(57)+rr(58)+rr(59)+
     &     rr(60)+rr(61)+rr(79)+rr(80)+rr(95)+rr(97)+rr(19)+
     &     rr(20)+rr(60)+rr(68)+rr(41)
      rl(5)=rr(17)+rr(18)+rr(94)
c     
c     ** H2O2 **
c     
      rp(6)=rr(2)+rr(6)+rr(7)+rr(8)+rr(18)
      rl(6)=rr(1)+rr(5)+rr(9)+rr(10)+rr(16)+rr(20)+rr(29)+
     &     rr(39)+rr(48)+rr(53)+rr(73)+rr(92)+rr(15)+
     &     rr(28)
c     
c     ** HCHO **
c     
      rp(7)= rr(65)+rr(67)+rr(68)+rr(69)+rr(70)+rr(71)+rr(102)
     &     +rr(107)+rr(108)
      rl(7)= rr(106)+rr(50)+rr(51)
c     
c     ** HCOOH **
c     
      rp(8)=rr(50)
      rl(8)=rr(52)+rr(53)+rr(54)+rr(55)+rr(56)+rr(79)+rr(97)+
     &     rr(57)+rr(58)+rr(59)+rr(60)+rr(61)+rr(80)+rr(95)
c     
c     ** NO **
c     
      rp(9)=rr(35)+rr(36)+rr(45)
      rl(9)=rr(31)+rr(33)
c     
c     ** NO2 **
c     
      rp(10)=rr(37)+rr(38)+rr(41)+rr(42)+rr(43)+rr(44)+rr(93)
      rl(10)=rr(31)+2.0d0*rr(32)+rr(34)+2.0d0*rr(104)
c     
c     ** O3 **
c     
      rp(11)=0.0d0
      rl(11)=rr(2)+rr(11)+rr(12)+rr(13)+rr(14)+rr(15)+rr(16)+
     &     rr(40)+rr(51)+rr(55)+rr(58)+rr(72)
c     
c     ** PAN **
c     
      rp(12)=0.0d0
      rl(12)=rr(62)+rr(98)
c     
c     ** CH3COOOH **
c     
      rp(13)=0.0d0
      rl(13)=rr(100)
c     
c     ** CH3OOH **
c     
      rp(14)=rr(63)+rr(64)
      rl(14)=rr(65)+rr(66)+rr(70)+rr(99)
c     
c     ** HCl **
c     
      rp(15)=rr(22)+2.d0*rr(25)+2.d0*rr(26)+rr(27)+rr(29)+
     &     2.d0*rr(30)+2.d0*rr(42)+2.d0*rr(56)+2.d0*rr(61)+
     &     2.d0*rr(69)+2.d0*rr(109)+2.d0*rr(28)
      rl(15)=rr(21)+rr(49)+rr(86)+rr(96)
c     
c     ** OH **
c     
      rp(16)=2.0d0*rr(1)+rr(9)+rr(10)+rr(12)+
     &     rr(13)+rr(15)+rr(22)+rr(30)+rr(35)+rr(36)+rr(44)+
     &     rr(55)+rr(58)+rr(65)+rr(91)+rr(101)
      rl(16)=rr(3)+rr(4)+rr(5)+rr(11)+rr(17)+rr(21)+rr(33)+
     &     rr(34)+rr(37)+rr(38)+rr(50)+rr(52)+rr(57)+rr(66)+
     &     rr(67)+rr(75)+rr(76)+rr(83)+rr(108)
c     
c     ** HO2 **
c     
      rp(17)=rr(5)+rr(11)+rr(20)+rr(28)+rr(29)+rr(48)+rr(50)+
     &     rr(52)+rr(54)+rr(55)+rr(56)+rr(57)+rr(59)+rr(60)+
     &     rr(61)+rr(65)+rr(67)+rr(68)+rr(69)+rr(71)+rr(79)+
     &     rr(92)+rr(95)+rr(97)+rr(102)+rr(14)+rr(14)+rr(15)+
     &     rr(58)+rr(80)
      rl(17)=rr(3)+2.0d0*rr(6)+rr(7)+rr(9)+rr(12)+rr(25)+
     &     rr(27)+rr(46)+rr(63)+rr(89)+rr(101)+rr(4)+rr(7)+
     &     2.0d0*rr(8)+rr(10)+rr(13)+rr(18)+rr(19)+rr(26)+
     &     rr(47)+rr(64)+rr(78)+rr(90)
c     
c     ** NO3 **
c     
      rp(18)=0.0d0
      rl(18)= rr(43)+rr(45)+rr(46)+rr(47)+rr(48)+rr(49)+rr(54)+
     &     rr(59)+rr(71)+rr(103)
c     
c     ** NH3 **
c     
      rp(19)=0.0d0
      rl(19)=0.0d0
c     
c     ** CH3O2 **
c     
      rp(20)=rr(66)
      rl(20)=rr(63)+rr(64)
c     
c     ** CH3OH **
c     
      rp(21)=0.0d0
      rl(21)=rr(67)+rr(68)+rr(69)+rr(71)+rr(102)
c     
c     ** Cl2-, Cl **
c     
      rp(22)=rr(49)+rr(96)+rr(23)
      rl(22)=rr(25)+rr(26)+rr(28)+rr(30)+rr(42)+rr(56)+rr(61)+
     &     rr(69)+rr(109)+rr(24)+rr(27)+rr(29)
c     
c     ** ClOH- **
c     
      rp(23)=rr(21)+rr(24)
      rl(23)=rr(22)+rr(23)
c     
c     ** SO4- **
c     
      rp(24)=2.d0*rr(81)+rr(103)
      rl(24)= rr(84)+rr(87)+rr(88)+rr(89)+rr(90)+rr(91)+
     &     rr(92)+rr(93)+rr(94)+rr(95)+rr(96)+rr(97)+rr(102)
c     
c     ** SO5- **
c     
      rp(25)=rr(75)+rr(76)+rr(83)+rr(84)+rr(87)+rr(88)+rr(108)+
     &     rr(109)
      rl(25)=rr(78)+rr(79)+rr(80)+2.0d0*rr(81)
c     
c     ** HSO5- **
c     
      rp(26)=rr(77)+rr(78)+rr(79)+rr(80)
      rl(26)=+rr(82)+rr(83)+rr(84)+rr(85)+rr(86)
c     
c     ** HOCH2SO3- **
c     
      rp(27)=rr(106)
      rl(27)=rr(107)+rr(108)
c     
c     ** CO3- **
c     
      rp(28)=rr(17)+rr(18)+rr(94)
      rl(28)=rr(19)+rr(20)+rr(41)+rr(60)+rr(68)

      return
      end
