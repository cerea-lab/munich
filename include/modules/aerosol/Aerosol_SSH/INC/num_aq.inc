C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Kathleen Fahey
C
C     This file is part of the Variable Size Resolved Model (VSRM),
C     based on the VSRM model of Carnegie Melon University. It is a
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

C     Numerical setuo for aqueous-phase model

      double precision TINYAQ,TINYAQ2,RTOLSULF
      double precision RTOLPH,PHDEF
	  double precision TINYNAQ	  

      integer NITSUBAQ,NITVSRM,NIT_PH

C     For time integration
      
      parameter (NITSUBAQ = 10) !subsycling timesteps
      parameter (NITVSRM  = 10)  !number of possible restarts if
                                !Sulfate conservation not met.

      parameter (TINYAQ   = 0.D0) !1.D-20) !tiny parameters for outputs
      parameter (TINYAQ2  = 0.D0) !1.D-12) !tiny parameters for outputs 
                                ! of AQOPERATOR.f
      parameter (TINYNAQ = 0.D0) !1.D-10) !tiny parameters for number

      parameter (RTOLSULF = 1.D-2) ! Tolerance error for sulfate
                                ! balance conservation


C     For PH computation

      parameter (RTOLPH  = 1.d-5) !RTOLPH : RELATIVE TOLERANCE 
                                  !for convergence
      parameter (PHDEF   = 4.5d0) !PHDEF  : default value 
                                  !if no convergence
      parameter (NIT_PH  = 1000)  !maximal number of iterations 
                                    !for bisection method

