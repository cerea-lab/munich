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

      subroutine fullequil(con,spres,cmet,akeq,akhen,wv,temp,xsol)
      
C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This routine computes the H+ concentration such that 
C     electroneutrality is met.
C     If no convergence occurs, a default value is used.
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     CON   : concentration vector ([...]).
C     SPRES : gas-phase concentration ([...]).
C     CMET  : metal concentrations    ([...]).
C     AKEQ  : kinetic rates for equilibrium.
C     AKHENR: Henry's rate.
C     WV    : water vapor ([...]).
C     TEMP  : temperature ([K]).
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     -- OUTPUT VARIABLES
C     
C     XSOL : H+ concentration.
C     
C------------------------------------------------------------------------
C     
C     -- REMARKS
C     
C------------------------------------------------------------------------
C     
C     -- MODIFICATIONS
C     
C     1) Delete all single/double conversions.
C     2) Remove the GOTO.
C     3) Define maximum number of iterations (NIT_PH) and default value
C     no convergence occurs (PHDEF).
C     4) Remove UU for the call to ELECTRO.
C     5) Remove include (for_aq.inc).
C     6) Change list of arguments for ELECTRO (add relative function).
C     7) Fix bug for BB for the initial interval.
C     8) Prevent for no convergence.
C     9) Loop from -13 (instead of -14).
C     10)Numerical setup in num_aq.inc
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

      include 'num_aq.inc'

      double precision con(28), spres(21), cmet(4)
      double precision akeq(17), akhen(21)
      double precision wv,temp,xsol

      double precision aa,bb,error
      double precision f,fa,fm,frel,farel,fmrel
      double precision x,xm

      integer i,ind_ok

      
C     1) Numerical setup
C     and default value
C     ------------------
C     defined in num_aq.inc

      xsol    = 10.0d0**(-phdef)

C     2) INITIAL INTERVAL [aa,bb] FOR THE BISECTION METHOD
C     ----------------------------------------------------

      x=10.0d0**(-14)
      call electro(x,con,spres,cmet,akeq,akhen,wv,temp,fa,farel)
      aa=x

      ind_ok=0

      DO i=-13,1
         IF (ind_ok.eq.0) THEN

            x=10.0d0**i
            CALL electro(x,con,spres,cmet,akeq,akhen,wv,temp,f,frel)

            IF (f*fa.ge.0.0d0) THEN
               aa=x
               fa=f
            ELSE
               bb=x
               ind_ok=1
            ENDIF

         ENDIF
      ENDDO

C     3) BISECTION METHOD on [aa,bb]
C     if bb has been found.
C     Default value otherwise.
C     ------------------------------

      IF (ind_ok.eq.1) THEN

         ind_ok=0

         DO  i=1,NIT_PH
            IF (ind_ok.eq.0) THEN

               error= dabs(bb-aa)/aa
               IF (error .le. rtolph) THEN
                  xsol   =(aa+bb)/2.0d0
                  ind_ok = 1
               ENDIF
               xm=(aa+bb)/2.0d0
               call electro(xm,con,spres,cmet,akeq,akhen,wv,temp,fm,
     &              fmrel)

               IF (fa*fm .gt.  0.0d0) THEN
                  aa=xm
                  fa=fm
               ELSE
                  bb=xm
               ENDIF

            ENDIF
         ENDDO

      ENDIF

      RETURN
      END


