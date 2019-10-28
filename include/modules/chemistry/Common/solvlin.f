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



      SUBROUTINE solvlin (ns, Kindlu,DLa,DLalu,DLx,DLb,option_chemistry)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine solves DLA * DLX = DLB where DLA is an input matrix,
C     and DLB is an input vector.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     KINDLU: 0 if DLALU is not a LU factorization of DLA. If KINDLU is
C     # not zero, DLALU is assumed to be a LU factorization of DLA.
C     DLA: matrix (NS x NS).
C     DLB: right-hand-side vector (NS) of the equation to be solved.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     DLALU: if KINDLU is not zero, DLALU is an LU factorisation of DLA.
C     # Otherwise, on exit, DLALU is an LU factorization of DLA.
C     IPVT: pivot indices; for 1 <= i <= NS, row i of the
C     # matrix was interchanged with row IPVT(i).
C
C     -- OUTPUT VARIABLES
C
C     DLX: solution of DLA * DLX = DLB.
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

      integer ns

      INTEGER Kindlu
      INTEGER Ji, Jj
      DOUBLE PRECISION DLa(ns,ns)
      DOUBLE PRECISION DLalu(ns,ns)
      DOUBLE PRECISION DLx(ns), DLb(ns)
      INTEGER option_chemistry


      DO Ji=1,ns
         DLx(Ji)=DLb(Ji)
      ENDDO

C------------------------------------------------------------------------
C     1 - Solve DLa * Dlx = Dlb

      IF (Kindlu .EQ. 0) THEN   ! DLalu is not
                                ! an LU factorization of DLa.
         DO Jj=1,ns
            DO Ji=1,ns
               DLalu(Ji,Jj)=DLa(Ji,Jj)
            ENDDO
         ENDDO

         IF (option_chemistry .eq. 1) then
            CALL LU_decompose_racm(ns,DLalu)
            CALL LU_solve_racm(ns,DLalu,DLx)
         ELSEIF (option_chemistry .eq. 2) then
            CALL LU_decompose_racm2(ns,DLalu)
            CALL LU_solve_racm2(ns,DLalu,DLx)
         ELSEIF (option_chemistry .eq. 3) then
            CALL LU_decompose_cb05(ns,DLalu)
            CALL LU_solve_cb05(ns,DLalu,DLx)
         ELSEIF (option_chemistry .eq. 4) then
            CALL LU_decompose_leighton(ns,DLalu)
            CALL LU_solve_leighton(ns,DLalu,DLx)
         ENDIF

      ELSE                      ! DLalu is an LU factorization of DLa.

         IF (option_chemistry .eq. 1) then
            CALL LU_solve_racm(ns,DLalu,DLx)
         ELSEIF (option_chemistry .eq. 2) then
            CALL LU_solve_racm2(ns,DLalu,DLx)
         ELSEIF (option_chemistry .eq. 3) then
            CALL LU_solve_cb05(ns,DLalu,DLx)
         ELSEIF (option_chemistry .eq. 4) then
            CALL LU_solve_leighton(ns,DLalu,DLx)
         ENDIF

      ENDIF

      END
