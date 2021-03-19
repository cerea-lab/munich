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


      SUBROUTINE compute_scavenging_coefficient_aer(Nsection,
     $     Temperature, Pressure, Rain, Diameter, Density,
     $     Zref, CloudBaseHeight,ScavengingCoefficient)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Compute scavenging coefficients for aerosol species
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     Nsection: Number of size section.
C     Temperature: Temperature in the cell. ([K])
C     Pressure: Pressure in the cell. ([Pa])
C     Rain: Rain intensity in the cell. ([mm/h])
C     Diameter: Representative diameter of aerosol particles. ([m])
C     Density: Density of aerosol particles. ([Kg/m^3])
C     RelativeHumidity: Relative humidity in the cell. ([0,1])
C     Zref: Height of the vertical node. ([m])
C     CloudBaseHeight: Height of cloud basis in the cell. ([m])
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     ScavengingCoefficient: Scavenging coefficient. ([1./s]).
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C     Parameterization used in this subroutine is presented and
C     justified in " Numerical and theoretical investigation of a
C     simplified model for the parameterization of below-cloud
C     scavenging by falling raindrops " (B. Sportisse and L. Du Bois,
C     Atmospheric Environment, Vol. 36, pp 5719-5727, 2002 )
C
C------------------------------------------------------------------------

      IMPLICIT NONE

C     -- Input.
      INTEGER Nsection
      DOUBLE PRECISION Temperature
      DOUBLE PRECISION Pressure
      DOUBLE PRECISION Rain
      DOUBLE PRECISION Diameter(Nsection)
      DOUBLE PRECISION Density(Nsection)
      DOUBLE PRECISION Zref
      DOUBLE PRECISION CloudBaseHeight

C     -- Local variables.
C     Pi.
      DOUBLE PRECISION PI
C     Air mean free path. ([m])
      DOUBLE PRECISION AirMeanFreePath
C     Kinetic viscosity.
      DOUBLE PRECISION Nuair
C     Air Density.
      DOUBLE PRECISION AirDensity
      DOUBLE PRECISION Muair
      DOUBLE PRECISION DropDiameter
      DOUBLE PRECISION DropVelocity
      DOUBLE PRECISION SettlingVelocity
      DOUBLE PRECISION DLrey
      DOUBLE PRECISION CC,DLdb,DLreypart
      DOUBLE PRECISION DLtau,DLSt,DLSh
      DOUBLE PRECISION DLSstar,DLphi,DLome,DLE
C     Water density and dynamic viscosity.
      DOUBLE PRECISION WaterDensity,Muwater

      INTEGER Is

C     -- Output.
      DOUBLE PRECISION ScavengingCoefficient(Nsection)

C------------------------------------------------------------------------
C     0. Setup

      PI=3.14159265358979323846D0

      DO Is=1,Nsection
         ScavengingCoefficient(Is)=0.D0
      ENDDO

c     [J/Kg/K] = [atm.L/mol/K] * [Pa/atm] / [Kg/mol] / [L/m^3]
      WaterDensity = 1.d3
      Muwater = 1.1d-3

      IF (Rain.GT.0.D0) THEN
         IF (CloudBaseHeight.GT.Zref) THEN
C------------------------------------------------------------------------
C     1. Computation of dynamic viscosity through Sutherland's law
C     and of air mean free path (in meter).

            call compute_AIR_FREE_MEAN_PATH(Temperature,
     &           Pressure,AirMeanFreePath,Muair)
            AirMeanFreePath = AirMeanFreePath*1.D-6

C     Mair/R USI = 28.97d-3/8.314 .

            AirDensity = Pressure * 3.48d-3 / Temperature

            Nuair = Muair / AirDensity

C------------------------------------------------------------------------
C     2. Computation of rain drop diameter, velocity and Reynolds number.

c     DropDiameter in cm
C     DropDiameter = 0.1d0 * 0.976d0 * Rain**0.21d0
            DropDiameter = 0.97d-1 * Rain**0.158d0

c     DropVelocity in m/s
C     DropVelocity = 9.58d0 * (1.d0 - DEXP(
C     &           (- (( DropDiameter / 0.171d0 )**1.147d0) ) ) )
            DropVelocity = 485.4d0 * DropDiameter
     &           * DEXP(-19.5d0*DropDiameter)

            DLrey = 1.d-2 * DropVelocity * DropDiameter / Nuair

C------------------------------------------------------------------------
C     2. Computation of raindrop-aerosol collision efficiency.

            DO Is=1,Nsection

C     Cunningham factor.
               call compute_CC
     &              (AirMeanFreePath,Diameter(Is),CC)

C     Gravitationnal settling velocity.
               call compute_VSTOKES(Diameter(Is),Density(Is),
     &              CC, Muair, SettlingVelocity)

C     Aerosol diffusivity.

               DLdb = 1.381d-23 * Temperature * CC
     &              / (3.d0 * PI * Muair * Diameter(Is))

C     Reynolds Number.

               DLreypart = DLrey / 2.d0

C     Relaxation time and Stokes number.

               DLtau =  Diameter(Is)**2.D0 * (Density(Is)-AirDensity)
     &              * CC / Muair / 18.d0
               DLSt = 2.d2* DLtau*(DropVelocity-SettlingVelocity)
     &              / DropDiameter

C     Schmidt Number.

               DLSh = Nuair / DLdb
C
               DLSstar = (1.2d0+(1.d0/12.d0)
     &              * DLOG(1.d0+DLreypart))
     &              / ( 1.d0 + DLOG(1.d0+DLreypart))

               IF(DLSt.LT.DLSstar)THEN
                  DLSt=DLSstar
               ENDIF

               DLphi = 100.d0 * Diameter(Is)/ DropDiameter
               DLome = Muair / Muwater

C     Efficiency coefficient.

               DLE = 4.d0 / ( DLreypart * DLSh )
     &              * ( 1.d0 + 0.4d0 * DLreypart**0.5d0
     &              * DLSh**(1.d0/3.d0)
     &              + 0.16d0*DLreypart**0.5d0*
     &              DLSh**(0.5d0))

               DLE = DLE + 4.d0 * DLphi
     &              * (DLome+(1.d0+2.d0*DLreypart
     &              **0.5d0 )* DLphi )

               DLE = DLE + ( (DLSt - DLSstar)
     &              / (DLSt - DLSstar + 2.d0/3.d0))**(1.5d0)
     &              * (Density(Is)/WaterDensity)**(0.5d0)

C     Scavenging coefficient for aerosol.

               ScavengingCoefficient(Is) = 1.5d-1 * DLE
     &              * Rain / (3600.d0 * DropDiameter)

            ENDDO               ! Loops bins.

         ENDIF
      ENDIF

      END
