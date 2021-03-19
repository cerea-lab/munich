C-----------------------------------------------------------------------
C     Copyright (C) 2007, ENPC - INRIA - EDF R&D
C     Author(s): Kathleen Fahey, Maryline Tombette
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

      subroutine compute_ph(conc_in, henrys, cst_dissoc, temp, lwc, pH,
     $     conc_out)

C-----------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes the H+ concentration such that
C     electroneutrality is met.
C     If no convergence occurs, a default value is used.
C
C-----------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     CONC_IN   : concentration vector.
C                 conc_in(1): S(IV) ([atm]).
C                 conc_in(2): total sulfate ([atm]).
C                 conc_in(3): total ammonium ([atm]).
C                 conc_in(4): total nitrate ([atm]).
C                 conc_in(5): total CO2 ([atm]).
C     henrys : Henry's constant associated with CONC_IN,
C              ([mol.L-1.K-1]).
C     cst_dissociation : dissociation constants ([mol.L-1]),
C                       (see droppar_simple.inc).
C     temp : Temperature ([K]).
C     lwc : Liquid Water Content ([vol water / vol air]).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     pH : pH.
C     conc_out: output concentrations of ions ([mol/L_water]),
C               partial pressures([atm]),
C               and aqueous species(mol/L_water]):
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
C     Adapted for the simple aqueous module (Marilyne Tombette, 2007).
C     Modified to improve the numerical behaviour (Yelva Roustan, 2015).
C
C-----------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Kathleen Fahey & Maryline Tombette, CEREA.
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

      include 'num_aq.inc'
      include 'droppar_simple.inc'

      double precision conc_in(5), conc_out(10)
      double precision akeq(17)
      double precision xsol, pH, temp,lwc
      double precision henrys(ngas_aq)
      double precision cst_dissoc(nreact_dissoc)
      double precision cst_oxydation(nreact_oxydation)

      double precision aa, bb, error
      double precision f, fa, fb, frel, farel, fbrel
      double precision x, xm

      integer i, ind_ok

C     1) Numerical setup and default value.
C     -------------------------------------
      xsol = 10.d0**(-phdef)   ! in mol L-1_water

C     2) Initial interval [aa,bb] for the bisection method.
C     -----------------------------------------------------
      aa = 10.d0**(-14)   ! in mol L-1_water

      call electro_simple(aa, conc_in, henrys, cst_dissoc,
     $     temp, lwc, fa, farel, conc_out)

      ind_ok = 0

      do i = -13, 1
         if (ind_ok.eq.0) then
            x = 10.0d0**i

            call electro_simple(x, conc_in, henrys, cst_dissoc,
     $           temp, lwc, f, frel, conc_out)

            if ((f * fa).ge.0.d0) then
               aa = x
               fa = f
               farel = frel
            else
               bb = x
               fb = f
               fbrel = frel
               ind_ok = 1
            endif
         endif
      enddo

C     3) Bisection method on [aa,bb] if bb has been found.
C     Default value otherwise.
C     ----------------------------------------------------

      if (ind_ok.eq.1) then

         ind_ok = 0

         do  i = 1, NIT_PH
            if (ind_ok.eq.0) then
               error = (fbrel + farel) / 2.d0

               if (error.le.rtolph) then
                  xsol = (aa + bb) / 2.d0
                  ind_ok = 1
               else
                  x = (aa + bb) / 2.d0

                  call electro_simple(x, conc_in, henrys, cst_dissoc,
     $                 temp, lwc, f, frel, conc_out)

                  if ((fa * f).gt.0.d0) then
                     aa = x
                     fa = f
                     farel = frel
                  else
                     bb = x
                     fb = f
                     fbrel = frel
                  endif
               endif

            endif
         enddo

      endif

      conc_out(1) = xsol
      ph = -dlog10(xsol)

      return
      end
