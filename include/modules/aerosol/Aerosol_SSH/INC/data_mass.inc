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
      double precision mmNH3, mmHNO3, mmHCl 
      double precision mmSO2,mmH2SO4
      double precision mmSO4,mmHSO5,mmHMSA
      double precision mmH2O2,mmHCHO, mmHNO2,mmHCOOH
      double precision mmO3,mmOH,mmHO2,mmNO3,mmNO
      double precision mmNO2,mmPAN,mmN2O5,mmNH4
      double precision mmN,mmS

      parameter (mmNH3  = 17.031d0)
      parameter (mmHNO3 = 63.013d0)
      parameter (mmHCl  = 36.461d0)
      parameter (mmSO2  = 64.064d0)
      parameter (mmH2SO4= 98.0d0)
      parameter (mmH2O2 = 34.015d0)
      parameter (mmHCHO = 30.026d0)
      parameter (mmHNO2 = 46.026d0)
      parameter (mmO3   = 47.998d0)
      parameter (mmOH   = 17.007d0)
      parameter (mmHO2  = 33.007d0)
      parameter (mmNO3  = 62.005d0)
      parameter (mmNO   = 30.006d0)
      parameter (mmNO2  = 46.006d0)	
      parameter (mmPAN  = 121.049d0)
      parameter (mmN2O5 = 108.0104d0)
      parameter (mmNH4  = 18.0383d0)
      parameter (mmSO4  = 96.063d0)
      parameter (mmHSO5 = 113.0704d0)
      parameter (mmHMSA = 111.098d0)
      parameter (mmN    = 14.0067d0)
      parameter (mmS    = 32.066d0)
      parameter (mmHCOOH= 46.03d0)
C     
C     parameter (mm   = )

