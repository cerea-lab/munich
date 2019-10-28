C-----------------------------------------------------------------------
C     Copyright (C) 2001-2007, ENPC - INRIA - EDF R&D
C
C     This file is part of the air quality modeling system Polyphemus.
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



      SUBROUTINE roschem (ns, nr, DLconc,ZCsourc,ZCsourcf,
     $     convers_factor, convers_factor_jac,
     s     ts,tf,DLRki,DLRkf,DLconc_old,DLk1,DLk2,
     $     option_chemistry)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes one timestep for gas-phase chemistry RACM
C     in one grid cell. The solver is the second-order Rosenbrock
C     method (see the user's guid). The linear systems to be
C     solved are optimized.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     TS: initial time (GMT, computed from January 1st, [s]).
C     TF: final time (GMT, computed from January 1st, [s]).
C     ZCSOURC: array of chemical volumic emissions at initial time.
C     ZCSOURCF: array of chemical volumic emissions at final time.
C     DLRKI: kinetic rates at initial time.
C     DLRKF: kinetic rates at final time.
C     JTESTI/J/K: index I/J/K of grid cell (for clipping report).
C
C     -- INPUT/OUTPUT VARIABLES
C
C     DLCONC: array of chemical concentrations ([\mu.g/m^3]).
C     # Before entry, it is given at initial time of the timestep.
C     # On exit, it is computed at final time of the timestep.
C
C     -- OUTPUT VARIABLES
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
C     Denis Quélo, CEREA, June 2001.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      integer ns, nr

      INTEGER Ji, Jj
      DOUBLE PRECISION DLconcbis(ns), DLconc(ns),DLr(ns)
      DOUBLE PRECISION DLstep,tf,ts
      DOUBLE PRECISION Igamma
      DOUBLE PRECISION DLk1(ns), DLk2(ns)
      DOUBLE PRECISION DLb1(ns), DLb2(ns)
      DOUBLE PRECISION DLmat(ns,ns)
      DOUBLE PRECISION DLmatlu(ns,ns)
      DOUBLE PRECISION DLdrdc(ns,ns)
      DOUBLE PRECISION ZCsourc(ns),ZCsourcf(ns)
      DOUBLE PRECISION DLRki(Nr),DLRkf(Nr)

      double precision convers_factor(ns),
     $     convers_factor_jac(ns,ns)

      DOUBLE PRECISION DLconc_new(ns)
      DOUBLE PRECISION Threshold

      DOUBLE PRECISION DLconc_old(ns)
      INTEGER option_chemistry

      Threshold = 0.D0

C     Numerical initialization.

      DLstep=tf-ts
      Igamma = 1.D0 + 1.D0/DSQRT(2.D0)

      Do Ji=1,ns
        DLconc_old(Ji) = DLconc(Ji)
      Enddo

C------------------------------------------------------------------------
C     1 - First step

C     Compute chemical production terms at initial time (DLb1).

      IF (option_chemistry .eq. 1) then
         CALL fexchem_racm (ns,nr,DLconc,DLRki,
     s        ZCsourc,convers_factor,DLb1)
      ELSEIF (option_chemistry .eq. 2) then
         CALL fexchem_racm2 (ns,nr,DLconc,DLRki,
     s        ZCsourc,convers_factor,DLb1)
      ELSEIF (option_chemistry .eq. 3) then
         CALL fexchem_cb05 (ns,nr,DLconc,DLRki,
     s        ZCsourc,convers_factor,DLb1)
      ELSEIF (option_chemistry .eq. 4) then
         CALL fexchem_leighton (ns,nr,DLconc,DLRki,
     s        ZCsourc,convers_factor,DLb1)
      ENDIF   

C     Compute the Jacobian at initial time (DLRDC).

      IF (option_chemistry .eq. 1) then
         CALL jacdchemdc_racm (ns, nr, DLconc,convers_factor,
     $     convers_factor_jac,DLRKi,DLdrdc)
      ELSEIF (option_chemistry .eq. 2) then
         CALL jacdchemdc_racm2 (ns, nr, DLconc,convers_factor,
     $     convers_factor_jac,DLRKi,DLdrdc)
      ELSEIF (option_chemistry .eq. 3) then
         CALL jacdchemdc_cb05 (ns, nr, DLconc,convers_factor,
     $     convers_factor_jac,DLRKi,DLdrdc)
      ELSEIF (option_chemistry .eq. 4) then
         CALL jacdchemdc_leighton (ns, nr, DLconc,convers_factor,
     $     convers_factor_jac,DLRKi,DLdrdc)
      ENDIF

C     Compute K1 by solving (1-Igamma*dt*DLRDC) DLK1 = DLR

      DO Jj=1,ns
         DO Ji=1,ns
            DLmat(Ji,Jj) = - Igamma*DLstep*DLdrdc(Ji,Jj)
         ENDDO
         DLmat(Jj,Jj) = 1.D0 + DLmat(Jj,Jj)
      ENDDO

      CALL solvlin(ns,0,DLmat,DLmatlu,DLk1,DLb1,
     $     option_chemistry)

C------------------------------------------------------------------------
C     2 - Second step

C     Compute first-order approximation (DLCONCBIS).

      DO Ji=1,ns
         DLconcbis(Ji) = DLconc(Ji) + DLstep * DLk1(Ji)
         IF (DLconcbis(Ji) .LE. Threshold) THEN
            DLconcbis(Ji) = Threshold
            DLk1(Ji) = (DLconcbis(Ji) - DLconc(Ji)) / DLstep
         ENDIF
      ENDDO

C     Compute chemical production terms (DLR) at final time with
C     the first-order approximation.

      IF (option_chemistry .eq. 1) then
         CALL fexchem_racm (ns, nr, DLconcbis,DLRkf,
     $     ZCsourcf,convers_factor, DLr)
      ELSEIF (option_chemistry .eq. 2) then
         CALL fexchem_racm2 (ns, nr, DLconcbis,DLRkf,
     $     ZCsourcf,convers_factor, DLr)
      ELSEIF (option_chemistry .eq. 3) then
         CALL fexchem_cb05 (ns, nr, DLconcbis,DLRkf,
     $     ZCsourcf,convers_factor, DLr)
      ELSEIF (option_chemistry .eq. 4) then
         CALL fexchem_leighton (ns, nr, DLconcbis,DLRkf,
     $     ZCsourcf,convers_factor, DLr)
      ENDIF

C     Compute DLK2 by solving (optimized scheme)
C     (1-Igamma*dt*DLRDC) DLK2 = DLR-2 DLK1

      DO Ji=1,ns
         DLb2(Ji) =  DLr(Ji) - 2.D0*DLk1(Ji)
      ENDDO

      CALL solvlin(ns,1,DLmat,DLmatlu,DLk2,DLb2,
     $     option_chemistry)


C------------------------------------------------------------------------
C     3 - Compute concentrations at final time (optimized scheme):
C     DLCONC+(3*DLRK1+DLRK2)*dt/2

      DO Ji=1,ns
         DLconc_new(Ji) = DLconc(Ji) + 1.5D0 * DLstep * DLk1(Ji)
     &        + 0.5D0 * DLstep * DLk2(Ji)
         IF (DLconc_new(Ji) .LE. Threshold) THEN
            DLconc(Ji) = Threshold
         ELSE
            DLconc(Ji) = DLconc_new(Ji)
         ENDIF
      ENDDO

      END
