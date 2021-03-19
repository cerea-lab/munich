C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Edouard Debry
C
C     This file is part of the Size Resolved Aerosol Model (SIREAM),
C     which is a component of the air quality modeling system 
C     Polyphemus.
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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       
C	RUNNING POINTERS                           
C       
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        INTEGER E1,E2,ICG2,ICE2,INU2,
     &	KDSLV2,IKV2,ICUT2,EH2O,EMD,EBC,ENa,ECl,ESO4,
     &	ENH3,ENO3
	INTEGER nesp_isorropia,nesp_aec,nesp_pankow
	INTEGER nesp_pom
	INTEGER isorropia_species(NEXT)
	INTEGER aec_species(NEXT)
	INTEGER pankow_species(NEXT)
	INTEGER poa_species(NEXT)
	INTEGER aerosol_species_interact(NEXT)
	 
	COMMON /varpnt/
     &  isorropia_species,aec_species,pankow_species,
     &  poa_species,aerosol_species_interact,
     &  E1,E2,ICG2,ICE2,INU2,KDSLV2,IKV2,ICUT2,
     &  EH2O,EMD,EBC,ENa,ECl,ESO4,ENH3,ENO3,
     &	nesp_isorropia,nesp_aec,nesp_pankow,nesp_pom
!$OMP   THREADPRIVATE(/varpnt/)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
