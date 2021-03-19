C-----------------------------------------------------------------------
C     Copyright (C) 2007, ENPC - INRIA - EDF R&D
C     Author(s): Maryline Tombette
C
C     This file is part of the Simple Aqueous model (SIMPLE_AQUEOUS), a
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

      subroutine electro_simple(x, conc_in, henrys, cst_dissoc, temp,
     $     lwc, f, frel, conc_out)

C-----------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes the electroneutrality balance for the
C     aqueous-phase model.
C
C-----------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     X    : H+ concentration. [mol.L-1].
C     CONC_IN  : concentration vector:
C                conc_in(1): S(IV) ([atm]).
C                conc_in(2): total sulfate ([atm]).
C                conc_in(3): total ammonium ([atm]).
C                conc_in(4): total nitrate ([atm]).
C                conc_in(5): total CO2 ([atm]).
C     cst_dissoc : equilibrium rate constants ([mol.L-1]).
C     henrys : Henry's rates ([mol.L-1.K-1]).
C     temp : temperature ([K]).
C     lwc : liquid water content ([vol water / vol air]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     F : evaluation of the electroneutrality relation (anion-cation).
C     FL: relative value (ABS(F)/(anion+cation).
C     conc_out: output concentrations of ions ([mol/L_water]), 
C               partial pressures ([atm]),
C               and aqueous species (mol/L_water]):
C               conc_out(1): [H+].
C               conc_out(2): [SO2(aq)].
C               conc_out(3): [HSO3-].
C               conc_out(4): [SO3^2-].
C               conc_out(5): [HNO3(aq)].
C               conc_out(6): [NO3-].
C               conc_out(7): [NH3(aq)].
C               conc_out(8): [NH4+].
C               conc_out(9): partial pressure HNO3(g).
C               conc_out(10): partial pressure NH3(g).
C
C-----------------------------------------------------------------------
C
C     -- REMARKS
C
C-----------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C     Added comments to ease the reading/checking (Yelva Roustan, 2015).
C
C-----------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Marilyne Tombette, CEREA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

      include 'droppar_simple.inc'
      include 'data_mass.inc'

C     Input variables.
      double precision conc_in(5), henrys(ngas_aq)
      double precision x, temp, lwc
      double precision cst_dissoc(nreact_dissoc)
      double precision cst_oxydation(nreact_oxydation)

C     Output variables.
      double precision f, frel
      double precision conc_out(10)

C     Local variables.
      double precision cOHm
      double precision ppHNO3
      double precision cHNO3, cNO3m
      double precision ppSO2
      double precision cSO2, cHSO3m, cSO32m
      double precision cH2SO4, cHSO4m, cSO42m
      double precision ppNH3
      double precision cNH3, cNH4p
      double precision cCO2, cHCO3m, cCO32m
      double precision fanion, fcation
      double precision cst_RT

      cst_RT = 8.206d-2 * temp  ! in atm.L_air.mol-1

C     1)Compute ions in aqueous-phase.
C     --------------------------------

C     * Dissociation H2O -> H+ + OH-
C     -  [OH-] = cst_dissoc(H2O) / [H+]
      cOHm = cst_dissoc(1) / x  ! [OH-] in mol.L-1_water

C     [HNO3]_tot = [HNO3]_g + [HNO3]_a + [NO3-]_a (all in atm)
C     with:
C     *  [X]_a the concentration in the liquid in the volume of air (in atm).
C     -  [X]_a = [X]_aq * cst_RT * lwc
C     -  [X]_aq the concentration in the liquid (in mol.L-1_water).
C     * Henry equilibrium
C     -  [HNO3]_g = [HNO3]_aq / henrys(HNO3)
C     * Dissociation HNO3 -> NO3- + H+
C     -  [HNO3]_aq = ([NO3-]_aq * [H+]) / cst_dissoc(HNO3)
      ppHNO3 = conc_in(4) /
     $     (1.d0 + cst_RT * lwc * henrys(ighno3)
     $     * (1.d0 + cst_dissoc(7) / x)) ! [HNO3]_g in atm
      cHNO3 = ppHNO3 * henrys(ighno3) ! [HNO3]_aq in mol.L-1_water
      cNO3m = cHNO3 * cst_dissoc(7) / x ! [NO3-]_aq in mol.L-1_water

C     [NH3]_tot = [NH3]_g + [NH3]_a + [NH4+]_a (all in atm)
C     with:
C     *  [X]_a the concentration in the liquid in the volume of air (in atm).
C     -  [X]_a = [X]_aq * cst_RT * lwc
C     -  [X]_aq the concentration in the liquid (in mol.L-1_water).
C     * Henry equilibrium
C     -  [NH3]_g = [NH3]_aq / henrys(NH3)
C     * Dissociation NH3.H2O -> NH4+ + OH-
C     -  [NH3]_aq = ([NH4+]_aq * [OH-]) / cst_dissoc(NH3)
      ppNH3 = conc_in(3) /
     $     (1.d0 + cst_RT * lwc * henrys(ignh3)
     $     * (1.d0 + cst_dissoc(6) * x / cst_dissoc(1))) ! [NH3]_g in atm
      cNH3 = henrys(ignh3) * ppNH3 ! [NH3]_aq in mol.L-1_water
      cNH4p = cst_dissoc(6) * x / cst_dissoc(1) * cNH3 ! [NH4+]_aq in mol.L-1_water

C     [SIV]_tot = [SO2]_g + [SO2]_a + [HSO3-]_a + [SO3--]_a (all in atm)
C     with:
C     * Henry equilibrium
C     -  [SO2]_g = [SO2]_aq / henrys(SO2)
C     * Dissociation SO2.H2O -> HSO3- + H+
C     -  [SO2]_aq = ([HSO3-]_aq * [H+]) / cst_dissoc(SO2)
C     * Dissociation HSO3- -> SO3-- + H+
C     -  [HSO3-]_aq = ([SO3--]_aq * [H+]) / cst_dissoc(HSO3-)
      cSO2 = conc_in(1) / (1.d0 / henrys(igso2)
     $     + cst_RT * lwc * (1.d0 + cst_dissoc(2) / x
     $     + cst_dissoc(2) * cst_dissoc(3) / (x * x))) ! [SO2]_aq in mol.L-1_water
      cHSO3m = cSO2 * cst_dissoc(2) / x ! [HSO3]_aq in mol.L-1_water
      cSO32m = cHSO3m * cst_dissoc(3) / x ! [SO3--]_aq in mol.L-1_water

C     [SVI]_tot = [H2SO4]_a + [HSO4-]_a + [SO4--]_a (all in atm)
C     with:
C     * Dissociation H2SO4.H2O -> HSO4- + H+
C     -  [H2SO4]_aq = ([HSO4-]_aq * [H+]) / cst_dissoc(H2SO4)
C     * Dissociation HSO4- -> SO4-- + H+
C     -  [HSO4-]_aq = ([SO4--]_aq * [H+]) / cst_dissoc(HSO4-)
      cH2SO4 = conc_in(2) / (cst_RT * lwc) / (1.d0 + cst_dissoc(4) / x
     $     + cst_dissoc(4) * cst_dissoc(5) / (x * x)) ! [H2SO4]_aq in mol.L-1_water
      cHSO4m = cH2SO4 * cst_dissoc(4) / x ! [HSO4-]_aq in mol.L-1_water
      cSO42m = cHSO4m * cst_dissoc(5) / x ! [SO4--]_aq in mol.L-1_water

C     [CO2]_tot = [CO2]_g + [CO2]_a + [HCO3-]_a + [CO3--]_a (all in atm)
C     with:
C     *  [X]_a the concentration in the liquid in the volume of air (in atm).
C     -  [X]_a = [X]_aq * cst_RT * lwc
C     -  [X]_aq the concentration in the liquid (in mol.L-1_water).
C     * Henry equilibrium
C     -  [CO2]_g = [CO2]_aq / henrys(CO2)
C     * Dissociation CO2.H2O -> HCO3- + H+
C     -  [CO2]_aq = ([HCO3-]_aq * [H+]) / cst_dissoc(CO2)
C     * Dissociation HCO3- -> CO3-- + H+
C     -  [HCO3-]_aq = ([CO3--]_aq * [H+]) / cst_dissoc(HCO3-)
      cCO2 = conc_in(5) * henrys(igco2)
     $     / (1.d0 + lwc * cst_RT * henrys(igco2)
     $     * (1.d0 + cst_dissoc(8) / x
     $     + cst_dissoc(8) * cst_dissoc(9) / (x * x))) ! [CO2]_aq in mol.L-1_water
      cHCO3m = cCO2 * cst_dissoc(8) / x ! [HCO3-]_aq in mol.L-1_water
      cCO32m = cHCO3m * cst_dissoc(9) / x ! [CO3--]_aq in mol.L-1_water

      conc_out(1) = x
      conc_out(2) = cSO2
      conc_out(3) = cHSO3m
      conc_out(4) = cSO32m
      conc_out(5) = cHNO3
      conc_out(6) = cNO3m
      conc_out(7) = cNH3
      conc_out(8) = cNH4p
      conc_out(9) = ppHNO3
      conc_out(10) = ppNH3

C     2) Compute electroneutrality balance.
C     -------------------------------------

      fanion = cOHm
     $     + cNO3m
     $     + cHSO3m + 2.0d0 * cSO32m
     $     + cHSO4m + 2.0d0 * cSO42m
     $     + cHCO3m + 2.0d0 * cCO32m

      fcation = x + cNH4p

      f = fanion - fcation

      if ((fanion + fcation).gt.0.d0) then
         frel = abs(f) / (fanion + fcation)
      else
         frel = -1.d0
      endif

      return
      end
