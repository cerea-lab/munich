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

      subroutine initactiv(NS,dsf_aero,ifirstact)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This routine computes the label for the first activated bin.
C     
C     The activated bins range from IFIRSTACT to NS.
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     -- INPUT/OUTPUT VARIABLES
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
C     2005/11/14, Bruno Sportisse, CEREA.
C     
C------------------------------------------------------------------------

      IMPLICIT NONE

      include 'aerpar.inc'
      include 'droppar.inc'

      integer i,indok
      integer ifirstact

      INTEGER NS
      DOUBLE PRECISION dsf_aero(NS)

      indok = 0
      ifirstact = NS + 1

      DO i = 1,NS
         IF ((dsf_aero(i).GE.dactiv).AND.(indok.EQ.0)) THEN
            ifirstact = i
            indok    = 1
         ENDIF
      ENDDO

      return
      end
