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
 
      subroutine dratedc_leighton  (ns,nr,rk,y,dw)
 
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes the derivative of reaction  rates.
C     This routine is automatically generated by SPACK.
C     Mechanism: Leighton            
C     Species: _ciCB05             
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
 
 
 
      dw(  1, 52) =  rk(  1)
      dw(  2, 43) =  rk(  2)
      dw(  3, 47) =  rk(  3) * Y( 48)
      dw(  3, 48) =  rk(  3) * Y( 47)
 
      RETURN
      END
 