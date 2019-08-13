C------------------------------------------------------------------------
C     Copyright (C) 2001-2008, ENPC - INRIA - EDF R&D
C
C     This file is part of the air quality modeling system Polyphemus.
C
C     Polyphemus is developed in the INRIA - ENPC joint project-team
C     CLIME and in the ENPC - EDF R&D joint laboratory CEREA.
C
C     Polyphemus is free software; you can redistribute i and/or modify
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
C------------------------------------------------------------------------
 
      subroutine dratedc_RACM(ns,nr,rk,y,dw)
 
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes the derivative of reaction  rates.
C     This routine is automatically generated by SPACK.
C     Mechanism: RACM                
C     Species: _ciRACM             
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     Ns: chemical species number.
C     NR: reaction number.
C     RK: kinetic rates.
C     Y: chemical concentrations.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     DW: derivative of reaction rates wrt Y.
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     SPACK.
C
C------------------------------------------------------------------------
 
      implicit none
 
      integer nr,ns
      double precision rk(nr),y(ns),dw(nr,ns)
 
 
 
      dw(  1, 72) =  rk(  1)
      dw(  2, 50) =  rk(  2)
      dw(  3, 50) =  rk(  3)
      dw(  4, 23) =  rk(  4)
      dw(  5, 55) =  rk(  5)
      dw(  6, 15) =  rk(  6)
      dw(  7, 68) =  rk(  7)
      dw(  8, 68) =  rk(  8)
      dw(  9, 27) =  rk(  9)
      dw( 10, 66) =  rk( 10)
      dw( 11, 66) =  rk( 11)
      dw( 12, 61) =  rk( 12)
      dw( 13, 16) =  rk( 13)
      dw( 14, 53) =  rk( 14)
      dw( 15, 28) =  rk( 15)
      dw( 16, 59) =  rk( 16)
      dw( 17, 51) =  rk( 17)
      dw( 18, 51) =  rk( 18)
      dw( 19, 54) =  rk( 19)
      dw( 20, 48) =  rk( 20)
      dw( 21, 70) =  rk( 21)
      dw( 22, 49) =  rk( 22)
      dw( 23, 24) =  rk( 23)
      dw( 24, 35) =  rk( 24)
      dw( 25, 35) =  rk( 25) * Y( 50)
      dw( 25, 50) =  rk( 25) * Y( 35)
      dw( 26,  7) =  rk( 26)
      dw( 27,  7) =  rk( 27)
      dw( 28,  7) =  rk( 28)
      dw( 29, 50) =  rk( 29) * Y( 62)
      dw( 29, 62) =  rk( 29) * Y( 50)
      dw( 30, 50) =  rk( 30) * Y( 71)
      dw( 30, 71) =  rk( 30) * Y( 50)
      dw( 31, 62) =  rk( 31) * Y( 71)
      dw( 31, 71) =  rk( 31) * Y( 62)
      dw( 32, 27) =  rk( 32) * Y( 62)
      dw( 32, 62) =  rk( 32) * Y( 27)
      dw( 33, 71) =  rk( 33) * Y( 71)
      dw( 33, 71) =  rk( 33) * Y( 71)
      dw( 34, 71) =  rk( 34) * Y( 71)
      dw( 34, 71) =  rk( 34) * Y( 71)
      dw( 35, 35) =  rk( 35) * Y( 67)
      dw( 35, 67) =  rk( 35) * Y( 35)
      dw( 36, 35) =  rk( 36) * Y( 72)
      dw( 36, 72) =  rk( 36) * Y( 35)
      dw( 37, 35) =  rk( 37) * Y( 72)
      dw( 37, 72) =  rk( 37) * Y( 35)
      dw( 38, 62) =  rk( 38) * Y( 67)
      dw( 38, 67) =  rk( 38) * Y( 62)
      dw( 39, 62) =  rk( 39) * Y( 72)
      dw( 39, 72) =  rk( 39) * Y( 62)
      dw( 40, 62) =  rk( 40) * Y( 68)
      dw( 40, 68) =  rk( 40) * Y( 62)
      dw( 41, 71) =  rk( 41) * Y( 67)
      dw( 41, 67) =  rk( 41) * Y( 71)
      dw( 42, 71) =  rk( 42) * Y( 72)
      dw( 42, 72) =  rk( 42) * Y( 71)
      dw( 43, 15) =  rk( 43)
      dw( 44, 71) =  rk( 44) * Y( 68)
      dw( 44, 68) =  rk( 44) * Y( 71)
      dw( 45, 62) =  rk( 45) * Y( 23)
      dw( 45, 23) =  rk( 45) * Y( 62)
      dw( 46, 62) =  rk( 46) * Y( 55)
      dw( 46, 55) =  rk( 46) * Y( 62)
      dw( 47, 62) =  rk( 47) * Y( 15)
      dw( 47, 15) =  rk( 47) * Y( 62)
      dw( 48, 50) =  rk( 48) * Y( 67)
      dw( 48, 67) =  rk( 48) * Y( 50)
      dw( 49, 50) =  rk( 49) * Y( 72)
      dw( 49, 72) =  rk( 49) * Y( 50)
      dw( 50, 67) =  rk( 50) * Y( 67)
      dw( 50, 67) =  rk( 50) * Y( 67)
      dw( 51, 68) =  rk( 51) * Y( 67)
      dw( 51, 67) =  rk( 51) * Y( 68)
      dw( 52, 68) =  rk( 52) * Y( 72)
      dw( 52, 72) =  rk( 52) * Y( 68)
      dw( 53, 68) =  rk( 53) * Y( 72)
      dw( 53, 72) =  rk( 53) * Y( 68)
      dw( 54,  9) =  rk( 54)
      dw( 55, 68) =  rk( 55) * Y( 68)
      dw( 55, 68) =  rk( 55) * Y( 68)
      dw( 56, 62) =  rk( 56)
      dw( 57, 62) =  rk( 57) * Y(  1)
      dw( 57,  1) =  rk( 57) * Y( 62)
      dw( 58, 52) =  rk( 58) * Y( 62)
      dw( 58, 62) =  rk( 58) * Y( 52)
      dw( 59, 19) =  rk( 59) * Y( 35)
      dw( 59, 35) =  rk( 59) * Y( 19)
      dw( 60, 49) =  rk( 60) * Y( 35)
      dw( 60, 35) =  rk( 60) * Y( 49)
      dw( 61, 20) =  rk( 61) * Y( 62)
      dw( 61, 62) =  rk( 61) * Y( 20)
      dw( 62, 21) =  rk( 62) * Y( 62)
      dw( 62, 62) =  rk( 62) * Y( 21)
      dw( 63,  2) =  rk( 63) * Y( 62)
      dw( 63, 62) =  rk( 63) * Y(  2)
      dw( 64,  3) =  rk( 64) * Y( 62)
      dw( 64, 62) =  rk( 64) * Y(  3)
      dw( 65,  4) =  rk( 65) * Y( 62)
      dw( 65, 62) =  rk( 65) * Y(  4)
      dw( 66, 11) =  rk( 66) * Y( 62)
      dw( 66, 62) =  rk( 66) * Y( 11)
      dw( 67, 40) =  rk( 67) * Y( 62)
      dw( 67, 62) =  rk( 67) * Y( 40)
      dw( 68, 44) =  rk( 68) * Y( 62)
      dw( 68, 62) =  rk( 68) * Y( 44)
      dw( 69, 12) =  rk( 69) * Y( 62)
      dw( 69, 62) =  rk( 69) * Y( 12)
      dw( 70, 19) =  rk( 70) * Y( 62)
      dw( 70, 62) =  rk( 70) * Y( 19)
      dw( 71, 13) =  rk( 71) * Y( 62)
      dw( 71, 62) =  rk( 71) * Y( 13)
      dw( 72, 14) =  rk( 72) * Y( 62)
      dw( 72, 62) =  rk( 72) * Y( 14)
      dw( 73,  5) =  rk( 73) * Y( 62)
      dw( 73, 62) =  rk( 73) * Y(  5)
      dw( 74,  6) =  rk( 74) * Y( 62)
      dw( 74, 62) =  rk( 74) * Y(  6)
      dw( 75, 30) =  rk( 75) * Y( 62)
      dw( 75, 62) =  rk( 75) * Y( 30)
      dw( 76, 66) =  rk( 76) * Y( 62)
      dw( 76, 62) =  rk( 76) * Y( 66)
      dw( 77, 61) =  rk( 77) * Y( 62)
      dw( 77, 62) =  rk( 77) * Y( 61)
      dw( 78, 59) =  rk( 78) * Y( 62)
      dw( 78, 62) =  rk( 78) * Y( 59)
      dw( 79, 24) =  rk( 79) * Y( 62)
      dw( 79, 62) =  rk( 79) * Y( 24)
      dw( 80, 51) =  rk( 80) * Y( 62)
      dw( 80, 62) =  rk( 80) * Y( 51)
      dw( 81, 54) =  rk( 81) * Y( 62)
      dw( 81, 62) =  rk( 81) * Y( 54)
      dw( 82, 49) =  rk( 82) * Y( 62)
      dw( 82, 62) =  rk( 82) * Y( 49)
      dw( 83, 62) =  rk( 83) * Y( 48)
      dw( 83, 48) =  rk( 83) * Y( 62)
      dw( 84, 62) =  rk( 84) * Y( 10)
      dw( 84, 10) =  rk( 84) * Y( 62)
      dw( 85, 16) =  rk( 85) * Y( 62)
      dw( 85, 62) =  rk( 85) * Y( 16)
      dw( 86, 53) =  rk( 86) * Y( 62)
      dw( 86, 62) =  rk( 86) * Y( 53)
      dw( 87, 28) =  rk( 87) * Y( 62)
      dw( 87, 62) =  rk( 87) * Y( 28)
      dw( 88, 26) =  rk( 88) * Y( 62)
      dw( 88, 62) =  rk( 88) * Y( 26)
      dw( 89, 25) =  rk( 89) * Y( 62)
      dw( 89, 62) =  rk( 89) * Y( 25)
      dw( 90, 70) =  rk( 90) * Y( 62)
      dw( 90, 62) =  rk( 90) * Y( 70)
      dw( 91, 66) =  rk( 91) * Y( 68)
      dw( 91, 68) =  rk( 91) * Y( 66)
      dw( 92, 61) =  rk( 92) * Y( 68)
      dw( 92, 68) =  rk( 92) * Y( 61)
      dw( 93, 51) =  rk( 93) * Y( 68)
      dw( 93, 68) =  rk( 93) * Y( 51)
      dw( 94, 54) =  rk( 94) * Y( 68)
      dw( 94, 68) =  rk( 94) * Y( 54)
      dw( 95, 49) =  rk( 95) * Y( 68)
      dw( 95, 68) =  rk( 95) * Y( 49)
      dw( 96, 48) =  rk( 96) * Y( 68)
      dw( 96, 68) =  rk( 96) * Y( 48)
      dw( 97, 30) =  rk( 97) * Y( 68)
      dw( 97, 68) =  rk( 97) * Y( 30)
      dw( 98, 11) =  rk( 98) * Y( 68)
      dw( 98, 68) =  rk( 98) * Y( 11)
      dw( 99, 40) =  rk( 99) * Y( 68)
      dw( 99, 68) =  rk( 99) * Y( 40)
      dw(100, 44) =  rk(100) * Y( 68)
      dw(100, 68) =  rk(100) * Y( 44)
      dw(101, 12) =  rk(101) * Y( 68)
      dw(101, 68) =  rk(101) * Y( 12)
      dw(102, 19) =  rk(102) * Y( 68)
      dw(102, 68) =  rk(102) * Y( 19)
      dw(103, 13) =  rk(103) * Y( 68)
      dw(103, 68) =  rk(103) * Y( 13)
      dw(104, 14) =  rk(104) * Y( 68)
      dw(104, 68) =  rk(104) * Y( 14)
      dw(105, 25) =  rk(105) * Y( 68)
      dw(105, 68) =  rk(105) * Y( 25)
      dw(106, 11) =  rk(106) * Y( 50)
      dw(106, 50) =  rk(106) * Y( 11)
      dw(107, 40) =  rk(107) * Y( 50)
      dw(107, 50) =  rk(107) * Y( 40)
      dw(108, 44) =  rk(108) * Y( 50)
      dw(108, 50) =  rk(108) * Y( 44)
      dw(109, 12) =  rk(109) * Y( 50)
      dw(109, 50) =  rk(109) * Y( 12)
      dw(110, 19) =  rk(110) * Y( 50)
      dw(110, 50) =  rk(110) * Y( 19)
      dw(111, 13) =  rk(111) * Y( 50)
      dw(111, 50) =  rk(111) * Y( 13)
      dw(112, 14) =  rk(112) * Y( 50)
      dw(112, 50) =  rk(112) * Y( 14)
      dw(113, 49) =  rk(113) * Y( 50)
      dw(113, 50) =  rk(113) * Y( 49)
      dw(114, 48) =  rk(114) * Y( 50)
      dw(114, 50) =  rk(114) * Y( 48)
      dw(115, 25) =  rk(115) * Y( 50)
      dw(115, 50) =  rk(115) * Y( 25)
      dw(116, 29) =  rk(116) * Y( 72)
      dw(116, 72) =  rk(116) * Y( 29)
      dw(117, 29) =  rk(117) * Y( 71)
      dw(117, 71) =  rk(117) * Y( 29)
      dw(118, 17) =  rk(118) * Y( 72)
      dw(118, 72) =  rk(118) * Y( 17)
      dw(119, 17) =  rk(119)
      dw(120, 17) =  rk(120) * Y( 50)
      dw(120, 50) =  rk(120) * Y( 17)
      dw(121, 18) =  rk(121) * Y( 72)
      dw(121, 72) =  rk(121) * Y( 18)
      dw(122, 18) =  rk(122)
      dw(123, 18) =  rk(123) * Y( 50)
      dw(123, 50) =  rk(123) * Y( 18)
      dw(124, 22) =  rk(124) * Y( 72)
      dw(124, 72) =  rk(124) * Y( 22)
      dw(125, 22) =  rk(125)
      dw(126, 22) =  rk(126) * Y( 50)
      dw(126, 50) =  rk(126) * Y( 22)
      dw(127, 69) =  rk(127) * Y( 72)
      dw(127, 72) =  rk(127) * Y( 69)
      dw(128, 26) =  rk(128)
      dw(129, 56) =  rk(129) * Y( 72)
      dw(129, 72) =  rk(129) * Y( 56)
      dw(130, 25) =  rk(130)
      dw(131, 63) =  rk(131) * Y( 67)
      dw(131, 67) =  rk(131) * Y( 63)
      dw(132, 60) =  rk(132) * Y( 67)
      dw(132, 67) =  rk(132) * Y( 60)
      dw(133, 64) =  rk(133) * Y( 67)
      dw(133, 67) =  rk(133) * Y( 64)
      dw(134, 36) =  rk(134) * Y( 67)
      dw(134, 67) =  rk(134) * Y( 36)
      dw(135, 37) =  rk(135) * Y( 67)
      dw(135, 67) =  rk(135) * Y( 37)
      dw(136, 38) =  rk(136) * Y( 67)
      dw(136, 67) =  rk(136) * Y( 38)
      dw(137, 41) =  rk(137) * Y( 67)
      dw(137, 67) =  rk(137) * Y( 41)
      dw(138, 45) =  rk(138) * Y( 67)
      dw(138, 67) =  rk(138) * Y( 45)
      dw(139, 39) =  rk(139) * Y( 67)
      dw(139, 67) =  rk(139) * Y( 39)
      dw(140, 42) =  rk(140) * Y( 67)
      dw(140, 67) =  rk(140) * Y( 42)
      dw(141, 43) =  rk(141) * Y( 67)
      dw(141, 67) =  rk(141) * Y( 43)
      dw(142, 31) =  rk(142) * Y( 67)
      dw(142, 67) =  rk(142) * Y( 31)
      dw(143, 32) =  rk(143) * Y( 67)
      dw(143, 67) =  rk(143) * Y( 32)
      dw(144, 33) =  rk(144) * Y( 67)
      dw(144, 67) =  rk(144) * Y( 33)
      dw(145, 69) =  rk(145) * Y( 67)
      dw(145, 67) =  rk(145) * Y( 69)
      dw(146, 56) =  rk(146) * Y( 67)
      dw(146, 67) =  rk(146) * Y( 56)
      dw(147, 58) =  rk(147) * Y( 67)
      dw(147, 67) =  rk(147) * Y( 58)
      dw(148, 47) =  rk(148) * Y( 67)
      dw(148, 67) =  rk(148) * Y( 47)
      dw(149, 46) =  rk(149) * Y( 67)
      dw(149, 67) =  rk(149) * Y( 46)
      dw(150, 63) =  rk(150) * Y( 71)
      dw(150, 71) =  rk(150) * Y( 63)
      dw(151, 60) =  rk(151) * Y( 71)
      dw(151, 71) =  rk(151) * Y( 60)
      dw(152, 64) =  rk(152) * Y( 71)
      dw(152, 71) =  rk(152) * Y( 64)
      dw(153, 36) =  rk(153) * Y( 71)
      dw(153, 71) =  rk(153) * Y( 36)
      dw(154, 37) =  rk(154) * Y( 71)
      dw(154, 71) =  rk(154) * Y( 37)
      dw(155, 38) =  rk(155) * Y( 71)
      dw(155, 71) =  rk(155) * Y( 38)
      dw(156, 41) =  rk(156) * Y( 71)
      dw(156, 71) =  rk(156) * Y( 41)
      dw(157, 45) =  rk(157) * Y( 71)
      dw(157, 71) =  rk(157) * Y( 45)
      dw(158, 39) =  rk(158) * Y( 71)
      dw(158, 71) =  rk(158) * Y( 39)
      dw(159, 42) =  rk(159) * Y( 71)
      dw(159, 71) =  rk(159) * Y( 42)
      dw(160, 43) =  rk(160) * Y( 71)
      dw(160, 71) =  rk(160) * Y( 43)
      dw(161, 31) =  rk(161) * Y( 71)
      dw(161, 71) =  rk(161) * Y( 31)
      dw(162, 32) =  rk(162) * Y( 71)
      dw(162, 71) =  rk(162) * Y( 32)
      dw(163, 33) =  rk(163) * Y( 71)
      dw(163, 71) =  rk(163) * Y( 33)
      dw(164, 69) =  rk(164) * Y( 71)
      dw(164, 71) =  rk(164) * Y( 69)
      dw(165, 69) =  rk(165) * Y( 71)
      dw(165, 71) =  rk(165) * Y( 69)
      dw(166, 56) =  rk(166) * Y( 71)
      dw(166, 71) =  rk(166) * Y( 56)
      dw(167, 56) =  rk(167) * Y( 71)
      dw(167, 71) =  rk(167) * Y( 56)
      dw(168, 58) =  rk(168) * Y( 71)
      dw(168, 71) =  rk(168) * Y( 58)
      dw(169, 47) =  rk(169) * Y( 71)
      dw(169, 71) =  rk(169) * Y( 47)
      dw(170, 46) =  rk(170) * Y( 71)
      dw(170, 71) =  rk(170) * Y( 46)
      dw(171, 63) =  rk(171) * Y( 63)
      dw(171, 63) =  rk(171) * Y( 63)
      dw(172, 60) =  rk(172) * Y( 63)
      dw(172, 63) =  rk(172) * Y( 60)
      dw(173, 64) =  rk(173) * Y( 63)
      dw(173, 63) =  rk(173) * Y( 64)
      dw(174, 36) =  rk(174) * Y( 63)
      dw(174, 63) =  rk(174) * Y( 36)
      dw(175, 37) =  rk(175) * Y( 63)
      dw(175, 63) =  rk(175) * Y( 37)
      dw(176, 38) =  rk(176) * Y( 63)
      dw(176, 63) =  rk(176) * Y( 38)
      dw(177, 41) =  rk(177) * Y( 63)
      dw(177, 63) =  rk(177) * Y( 41)
      dw(178, 45) =  rk(178) * Y( 63)
      dw(178, 63) =  rk(178) * Y( 45)
      dw(179, 39) =  rk(179) * Y( 63)
      dw(179, 63) =  rk(179) * Y( 39)
      dw(180, 42) =  rk(180) * Y( 63)
      dw(180, 63) =  rk(180) * Y( 42)
      dw(181, 43) =  rk(181) * Y( 63)
      dw(181, 63) =  rk(181) * Y( 43)
      dw(182, 31) =  rk(182) * Y( 63)
      dw(182, 63) =  rk(182) * Y( 31)
      dw(183, 32) =  rk(183) * Y( 63)
      dw(183, 63) =  rk(183) * Y( 32)
      dw(184, 33) =  rk(184) * Y( 63)
      dw(184, 63) =  rk(184) * Y( 33)
      dw(185, 69) =  rk(185) * Y( 63)
      dw(185, 63) =  rk(185) * Y( 69)
      dw(186, 69) =  rk(186) * Y( 63)
      dw(186, 63) =  rk(186) * Y( 69)
      dw(187, 56) =  rk(187) * Y( 63)
      dw(187, 63) =  rk(187) * Y( 56)
      dw(188, 56) =  rk(188) * Y( 63)
      dw(188, 63) =  rk(188) * Y( 56)
      dw(189, 58) =  rk(189) * Y( 63)
      dw(189, 63) =  rk(189) * Y( 58)
      dw(190, 47) =  rk(190) * Y( 63)
      dw(190, 63) =  rk(190) * Y( 47)
      dw(191, 46) =  rk(191) * Y( 63)
      dw(191, 63) =  rk(191) * Y( 46)
      dw(192, 60) =  rk(192) * Y( 69)
      dw(192, 69) =  rk(192) * Y( 60)
      dw(193, 64) =  rk(193) * Y( 69)
      dw(193, 69) =  rk(193) * Y( 64)
      dw(194, 36) =  rk(194) * Y( 69)
      dw(194, 69) =  rk(194) * Y( 36)
      dw(195, 37) =  rk(195) * Y( 69)
      dw(195, 69) =  rk(195) * Y( 37)
      dw(196, 38) =  rk(196) * Y( 69)
      dw(196, 69) =  rk(196) * Y( 38)
      dw(197, 41) =  rk(197) * Y( 69)
      dw(197, 69) =  rk(197) * Y( 41)
      dw(198, 45) =  rk(198) * Y( 69)
      dw(198, 69) =  rk(198) * Y( 45)
      dw(199, 39) =  rk(199) * Y( 69)
      dw(199, 69) =  rk(199) * Y( 39)
      dw(200, 42) =  rk(200) * Y( 69)
      dw(200, 69) =  rk(200) * Y( 42)
      dw(201, 43) =  rk(201) * Y( 69)
      dw(201, 69) =  rk(201) * Y( 43)
      dw(202, 31) =  rk(202) * Y( 69)
      dw(202, 69) =  rk(202) * Y( 31)
      dw(203, 32) =  rk(203) * Y( 69)
      dw(203, 69) =  rk(203) * Y( 32)
      dw(204, 33) =  rk(204) * Y( 69)
      dw(204, 69) =  rk(204) * Y( 33)
      dw(205, 69) =  rk(205) * Y( 69)
      dw(205, 69) =  rk(205) * Y( 69)
      dw(206, 56) =  rk(206) * Y( 69)
      dw(206, 69) =  rk(206) * Y( 56)
      dw(207, 58) =  rk(207) * Y( 69)
      dw(207, 69) =  rk(207) * Y( 58)
      dw(208, 47) =  rk(208) * Y( 69)
      dw(208, 69) =  rk(208) * Y( 47)
      dw(209, 46) =  rk(209) * Y( 69)
      dw(209, 69) =  rk(209) * Y( 46)
      dw(210, 47) =  rk(210) * Y( 47)
      dw(210, 47) =  rk(210) * Y( 47)
      dw(211, 47) =  rk(211) * Y( 46)
      dw(211, 46) =  rk(211) * Y( 47)
      dw(212, 46) =  rk(212) * Y( 46)
      dw(212, 46) =  rk(212) * Y( 46)
      dw(213, 63) =  rk(213) * Y( 68)
      dw(213, 68) =  rk(213) * Y( 63)
      dw(214, 60) =  rk(214) * Y( 68)
      dw(214, 68) =  rk(214) * Y( 60)
      dw(215, 64) =  rk(215) * Y( 68)
      dw(215, 68) =  rk(215) * Y( 64)
      dw(216, 36) =  rk(216) * Y( 68)
      dw(216, 68) =  rk(216) * Y( 36)
      dw(217, 37) =  rk(217) * Y( 68)
      dw(217, 68) =  rk(217) * Y( 37)
      dw(218, 38) =  rk(218) * Y( 68)
      dw(218, 68) =  rk(218) * Y( 38)
      dw(219, 41) =  rk(219) * Y( 68)
      dw(219, 68) =  rk(219) * Y( 41)
      dw(220, 45) =  rk(220) * Y( 68)
      dw(220, 68) =  rk(220) * Y( 45)
      dw(221, 39) =  rk(221) * Y( 68)
      dw(221, 68) =  rk(221) * Y( 39)
      dw(222, 42) =  rk(222) * Y( 68)
      dw(222, 68) =  rk(222) * Y( 42)
      dw(223, 43) =  rk(223) * Y( 68)
      dw(223, 68) =  rk(223) * Y( 43)
      dw(224, 31) =  rk(224) * Y( 68)
      dw(224, 68) =  rk(224) * Y( 31)
      dw(225, 32) =  rk(225) * Y( 68)
      dw(225, 68) =  rk(225) * Y( 32)
      dw(226, 33) =  rk(226) * Y( 68)
      dw(226, 68) =  rk(226) * Y( 33)
      dw(227, 69) =  rk(227) * Y( 68)
      dw(227, 68) =  rk(227) * Y( 69)
      dw(228, 56) =  rk(228) * Y( 68)
      dw(228, 68) =  rk(228) * Y( 56)
      dw(229, 58) =  rk(229) * Y( 68)
      dw(229, 68) =  rk(229) * Y( 58)
      dw(230, 47) =  rk(230) * Y( 68)
      dw(230, 68) =  rk(230) * Y( 47)
      dw(231, 46) =  rk(231) * Y( 68)
      dw(231, 68) =  rk(231) * Y( 46)
      dw(232, 65) =  rk(232) * Y( 71)
      dw(232, 71) =  rk(232) * Y( 65)
      dw(233, 65) =  rk(233) * Y( 63)
      dw(233, 63) =  rk(233) * Y( 65)
      dw(234, 65) =  rk(234) * Y( 69)
      dw(234, 69) =  rk(234) * Y( 65)
      dw(235, 65) =  rk(235) * Y( 65)
      dw(235, 65) =  rk(235) * Y( 65)
      dw(236, 65) =  rk(236) * Y( 67)
      dw(236, 67) =  rk(236) * Y( 65)
      dw(237, 65) =  rk(237) * Y( 68)
      dw(237, 68) =  rk(237) * Y( 65)
 
      RETURN
      END
 