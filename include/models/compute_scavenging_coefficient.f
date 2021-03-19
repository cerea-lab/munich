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


      subroutine compute_scavenging_coefficient(Nx, Ny, Nz, Nscav, Nxh,
     $     Nyh, Nzh, level, gas_phase_diffusivity, henry_constant,
     $     temperature, pressure, rain, cloud_height, lambda_scavenging)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Compute below-cloud scavenging coefficients for gaseous species.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     NX                    : number of cells along X.
C     NY                    : number of cells along Y.
C     NZ                    : number of cells along Z.
C     NXH                   : number of cells along X for Henry constant.
C     NYH                   : number of cells along Y for Henry constant.
C     NZH                   : number of cells along Z for Henry constant.
C     NSCAV                 : number of scavenged species.
C     LEVEL                 : height of vertical level. ([m])
C     GAS_PHASE_DIFFUSIVITY : gas phase diffusivity. ([cm^2/s])
C     HENRY_CONSTANT        : Henry constant. ([mol/L/atm])
C     TEMPERATURE           : temperature. ([K])
C     PRESSURE              : pressure. ([Pa])
C     RAIN                  : rain intensity. ([mm/h])
C     CLOUD_HEIGHT          : height of cloud basis. ([m])
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     LAMBDA_SCAVENGING: scavenging coefficient. ([1/s]).
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C     The parameterization used in this subroutine is presented and
C     justified in " Numerical and theoretical investigation of a
C     simplified model for the parameterization of below-cloud
C     scavenging by falling raindrops " (B. Sportisse and L. Du Bois,
C     Atmospheric Environment, Vol. 36, pp 5719-5727, 2002)
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C     Yelva Roustan, CEREA, december 2007.
C     Edouard Debry, CEREA, July 2007.
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Yelva Roustan, CEREA, February 2003.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      integer Nx, Ny, Nz, Nscav, Nxh, Nyh, Nzh
      integer i, j, k, ih, jh, kh, jsca

      double precision gas_phase_diffusivity(Nscav) !([cm^2/s])
      double precision henry_constant(Nxh, Nyh, Nzh, Nscav) !([mol/L/atm])
      double precision temperature(Nx, Ny, Nz)
      double precision pressure(Nx, Ny, Nz)
      double precision rain(Nx, Ny)
      double precision cloud_height(Nx, Ny)
      double precision level(Nz)

      double precision lambda_scavenging(Nx, Ny, Nz, Nscav)

C     air kinematic viscosity, volumic mass and dynamic viscosity.
      double precision eta_air, rho_air, nu_air
C     drop diameter, rain drop velocity.
      double precision d_drop, u_rain
C     Reynolds and Sherwood numbers.
      double precision reynolds, sherwood
C     mass transfer coefficient
      double precision k_transfer
C     perfect gas constant
      double precision Pr

      integer henry_3d

C------------------------------------------------------------------------
C     0. Setup

      if ((Nxh.eq.Nx).and.(Nyh.eq.Ny).and.(Nzh.eq.Nz)) then
         henry_3d = 1
      else if ((Nxh.eq.1).and.(Nyh.eq.1).and.(Nzh.eq.1)) then
         henry_3d = 0
      else
         write(6, *) "Error in compute_scavenging_coefficient: "
     $        // "incompatible dimensions of array for Henry constants."
         stop
      endif

      Pr = 8.20d-2              !Perfect gaz constant ([atm.L/mol/K]).

      do jsca = 1,Nscav
         do k = 1,Nz
            do j = 1,Ny
               do i = 1,Nx
                  lambda_scavenging(i, j, k, jsca) = 0.d0
               enddo
            enddo
         enddo
      enddo

      do j=1,Ny
         do i=1,Nx
            if (rain(i,j).gt.0.D0) then
               do k=1,Nz
                  if (cloud_height(i, j).gt.level(k)) then

C------------------------------------------------------------------------
C     1. Computation of kinematic viscosity with Sutherland's law.

                     eta_air = 1.83d-5 * ( 416.16d0 /
     $                    ( temperature(i, j, k) + 120.d0) ) *
     $                    ( temperature(i, j, k) / 296.16d0 )**1.5d0

C     Mair/R USI = 28.97d-3/8.314 .

                     rho_air = pressure(i, j, k) * 3.48d-3
     $                    / temperature(i, j, k)

                     nu_air = eta_air / rho_air

C------------------------------------------------------------------------
C     2. Computation of rain drop diameter, velocity and Reynolds number.

C     d_drop in cm
C     d_drop = 0.1d0 * 0.976d0 * rain(i,j)**0.21d0
                     d_drop = 0.97d-1 * rain(i,j)**0.158d0
C     u_rain in m/s
C     u_rain = 9.58d0 * (1.d0 - DEXP(
C     &                    (- ( d_drop / 0.171d0 )**1.147d0 ) ) )
                     u_rain = 485.4d0 * d_drop * DEXP(-19.5d0 * d_drop)

                     reynolds = 1.d-2 * u_rain * d_drop / nu_air

C------------------------------------------------------------------------
C     3. Computation of mass transfer coefficient.

                     do jsca=1,Nscav

C     Sherwood number.
                        sherwood = 2.d0 + 0.6d0 * reynolds**0.5d0
     $                       * ( 1.d4 * nu_air
     $                       / gas_phase_diffusivity(jsca))**(1.d0/3.d0)

C     Transfert coefficient ([cm/s]).
                        k_transfer = ( gas_phase_diffusivity(jsca)
     $                       / d_drop ) * sherwood

C------------------------------------------------------------------------
C     4. Computation of wet scavenging coefficient for gas.

C     1.67d-6 = 6 * 1.d-6 / 3.6d0.
                        ih = (i-1) * henry_3d + 1
                        jh = (j-1) * henry_3d + 1
                        kh = (k-1) * henry_3d + 1

                        lambda_scavenging(i, j, k, jsca) =
     $                       1.67d-6 * rain(i, j) /
     $                       ( u_rain * d_drop / k_transfer )
     $                       * dexp( -(cloud_height(i, j) - level(k))
     $                       / ( d_drop / 6.d0 / k_transfer * u_rain
     $                       * henry_constant(ih, jh, kh, jsca) * Pr
     $                       * temperature(i, j, k)) )

                     enddo
                  endif
               enddo
            endif
         enddo
      enddo

      end


C------------------------------------------------------------------------


      subroutine compute_scavenging_coefficient_pudykiewicz(
     $     Nx, Ny, Nz, Nscav, level, temperature, pressure,
     $     humidity, lambda_scavenging)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Compute in-cloud scaveging coefficients on the basis of
C     specific humidity for gaseous and particulate species.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     NX          : number of cells along X.
C     NY          : number of cells along Y.
C     NZ          : number of cells along Z.
C     NSCAV       : number of scavenged species.
C     LEVEL       : height of vertical level. ([m])
C     TEMPERATURE : temperature. ([K])
C     PRESSURE    : pressure. ([Pa])
C     HUMIDITY    : specific humidity. ([kg/kg])
C
C     -- INPUT/OUTPUT VARIABLES
C
C     -- OUTPUT VARIABLES
C
C     LAMBDA_SCAVENGING: scavenging coefficient. ([1/s])
C
C------------------------------------------------------------------------
C
C     -- REMARKS
C
C     The parameterization used in this subroutine is presented
C     in " Modelling transport and deposition of caesium and iodine from
C     the Chernobyl accident " (J. Brandt, J.H. Christensen and L.M.
C     Frohn, Atmos. Chem. Phys., Vol. 2, pp 397-417, 2002) and proposed
C     in " Simulation of the Chernobyl dispersion with a 3-D hemispheric
C     tracer model " (J. Pudykiewicz, Tellus, Vol. 41B, 391-412, 1989)
C
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C     Yelva Roustan, CEREA, 2007.
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Denis Quélo, CEREA, 2006.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      integer Nx, Ny, Nz, Nscav
      integer i, j, k, jsca

      double precision temperature(Nx, Ny, Nz)
      double precision pressure(Nx, Ny, Nz)
      double precision humidity(Nx, Ny, Nz)
      double precision level(Nz)

      double precision lambda_scavenging(Nx, Ny, Nz, Nscav)

      DOUBLE PRECISION relative_humidity(Nx, Ny, Nz)
      DOUBLE PRECISION A, esat
      DOUBLE PRECISION relative_humidity_min, relative_humidity_max
      DOUBLE PRECISION relative_humidity_threshold
      DOUBLE PRECISION tkelvin

C------------------------------------------------------------------------
C     0. Setup
      tkelvin = 273.15d0

      A = 3.5d-5

      relative_humidity_min = 0.d0
      relative_humidity_max = 1.d0
      relative_humidity_threshold = 0.8d0

C     Compute relative humidity.
      do k = 1,Nz
         do j = 1,Ny
            do i = 1,Nx
C     Saturation water pressure.
               esat = 610.78d0 * DEXP(17.2694d0
     $              * (Temperature(i, j, k) - tkelvin)
     $              / (Temperature(i, j, k) - 35.86d0))
               relative_humidity(i, j, k) = humidity(i, j, k)
     $              * Pressure(i, j, k) / esat
     $              / (0.622d0 + 0.378d0 * humidity(i, j, k))
               relative_humidity(i, j, k) = DMIN1(relative_humidity_max,
     $              relative_humidity(i, j, k))
               relative_humidity(i, j, k) = DMAX1(relative_humidity_min,
     $              relative_humidity(i, j, k))
            enddo
         enddo
      enddo

C     Compute in-cloud scavenging coefficient.
      do Jsca = 1,Nscav
         do k = 1,Nz
            do j = 1,Ny
               do i = 1,Nx
                  lambda_scavenging(i, j, k, Jsca) = A
     $                 * DMAX1((relative_humidity(i, j, k)
     $                 - relative_humidity_threshold), 0.d0)
     $                 /  (1.d0 - relative_humidity_threshold)
               enddo
            enddo
         enddo
      enddo

      end
