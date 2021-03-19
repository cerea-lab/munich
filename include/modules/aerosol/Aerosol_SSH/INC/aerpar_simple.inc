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

c******************************************************************
c     AEROSOL PARAMETERS 
c******************************************************************
c     
c     USEFUL CONSTANTS
c     
      INTEGER nas,naa
      INTEGER nan,nac,na4,naw,nae
      INTEGER nao,nar
      INTEGER ngtotal,naers
      
      DOUBLE PRECISION  dpmax
      DOUBLE PRECISION  dpmin,rho

c     number and distribution of sections
c     section distribution:
c     
c     |       |       |       | ...  |           |          |
c     dp(1)   dp(2)   dp(3)   ...   dp(NS-1)        dp(NS)
c     =dpmin                                       =dpmax
c     
c     dpmin and dpmax are the geometric average diameter of particles in the
c     smallest and largest bins.  Cut points are given by sqrt(dp(i)*dp(i+1)).
c     
c     USE 10 NM AS LOWER BOUND TO AVOID EVAPORATION PROBLEMS IN AQUEOUS MODULE      
c     
      parameter (dpmax = 10.0d-6) ! 10.0 um max. diameter [m] 
      parameter (dpmin = 0.05d-6) ! 0.1 um min. diam. [m]  NEW VALUE 
      parameter (rho = 1.4d12)  ! particle density [ug/m^3]
c     
c     AEROSOL COMPONENTS IN THE AEROSOL CONCENTRATION VECTOR
c     
      parameter (nas =  1)      ! sodium
      parameter (naa =  2)      ! ammonium
      parameter (nan =  3)      ! nitrate
      parameter (nac =  4)      ! chloride
      parameter (na4 =  5)      ! sulfate
      parameter (naw =  6)      ! water
      parameter (nae =  7)      ! elemental carbon
      parameter (nao =  8)      ! organics
      parameter (nar =  9)      ! crustal

c     total number of gas phase species so we know where the aerosol starts
c     
      parameter (ngtotal = 50)
      parameter (naers = 9)

      DOUBLE PRECISION liquid_density(naers)
c *********
c liquid_density(1) = 0.97D-06  ! sodium
c liquid_density(2) = 0.91D-06  ! ammonium
c liquid_density(3) = 1.50D-06  ! nitrate
c liquid_density(4) = 1.15D-06  ! chloride
c liquid_density(5) = 1.84D-06  ! sulfate
c liquid_density(6) = 1.00D-06  ! water
c liquid_density(7) = 2.25D-06  ! elemental carbon
c liquid_density(8) = 1.30D-06  ! organics
c liquid_density(9) = 2.94D-06  ! crustal ? a peu pres ?
      parameter (liquid_density = (/0.97D-06, 0.91D-06,
     & 1.50D-06, 1.15D-06, 1.84D-06, 1.00D-06, 2.25D-06, 
     & 1.30D-06, 2.94D-06/))
