C-----------------------------------------------------------------------
C     Copyright (C) 2005-2007, ENPC - INRIA - EDF R&D
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

      subroutine constants(temp)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This routine computes the kinetic/equilibrium rates for the aqueous-phase
C     mechanism (17 dissociation reactions, 21 Henry's reactions and
C     120 forward reactions).
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     TEMP   : temperature             [K].
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
C     1) Optimize computation by defining 1/TEMP-1./298.
C     2) Remove all "include" (not used).
C------------------------------------------------------------------------
C     
C     -- AUTHOR(S)
C     
C     2005/10/3, Kathleen Fahey, CEREA.
C     
C------------------------------------------------------------------------
c     
      IMPLICIT NONE 

      include 'aerpar.inc'
      include 'droppar.inc'
      include 'aqrates.inc'
      double precision dheq(17),dhhen(21),dhre(120)
      double precision bkeq(17),bkhen(21),bkre(120)
      double precision temp

      double precision coefloc
      integer i

      data dheq/1960.d0,1500.d0,0.,2720.d0,-3730.d0,8700.d0,-1260.d0,
     &     -1000.d0,-1760.d0,
     &     -450.d0,-6710.d0,4020.d0,-20.d0,6900.d0,0.d0,0.d0,0.d0/
      data bkeq/1.23d-2,6.61d-8,1.d3,1.02d-2,2.2d-12,15.4d0,5.1d-4,
     &     4.46d-7,4.68d-11,1.75d-5,1.0d-14,1.82d3,1.78d-4,1.74d6,
     &     3.5d-5,
     &     5.26d-6,2.0d-12/
      data dhhen/3120.d0,0.d0,4780.d0,8700.d0,2420.d0,6620.d0,
     &     6460.d0,5740.d0,1480.d0,
     &     2500.d0,2300.d0,5910.d0,6170.d0,5610.d0,2020.d0,5280.d0,
     &     6640.d0,8700.d0,3400.d0,
     &     5600.d0,4900.d0/
      data bkhen/1.23d0,0.d0,49.d0,2.1d5,3.4d-2,7.45d4,6.3d3,
     &     3.5d3,1.9d-3,
     &     0.01d0,1.13d-2,2.9d0,473.d0,227.d0,727.d0,25.d0,
     &     2000.d0,2.1d5,
     &     75.d0,6.d0,220.d0/


      data dhre/0.0d0,0.0d0,-1500.d0,-1500.d0,-1700.d0,-2365.d0,
     &     -1500.d0,0.0d0,0.0d0,
     &     0.0d0,0.0d0,0.0d0,-1500.d0,0.0d0,-2520.d0,0.0d0,-1910.d0,
     &     0.0d0,-1500.d0,-2820.d0,
     &     -1500.d0,0.0d0,0.0d0,0.0d0,-1500.d0,-1500.d0,-1500.d0,
     &     -3370.d0,0.0d0,-2160.d0,
     &     -1500.d0,-1500.d0,-1500.d0,-1500.d0,0.0d0,0.0d0,-1500.d0,
     &     -1500.d0,-6693.d0,-6950.d0,
     &     0.0d0,-1500.d0,-1500.d0,0.0d0,0.0d0,-1500.d0,-1500.d0,
     &     -2800.d0,-1500.d0,-1500.d0,
     &     0.0d0,-1500.d0,-5180.d0,-3200.d0,0.0d0,-4300.d0,-1500.d0,
     &     0.0d0,-1500.d0,-3400.d0,
     &     -2600.d0,0.0d0,-3000.d0,-1600.d0,0.0d0,-1700.d0,-1500.d0,
     &     -4500.d0,-4400.d0,-1800.d0,
     &     -2800.d0,0.0d0,-5530.d0,-5280.d0,-4430.d0,-13700.d0,
     &     -11000.d0,-13700d0,-11000.d0,
     &     -1500.d0,-1500.d0,-3100.d0,-1500.d0,-5300.d0,-4000.d0,
     &     -1500.d0,-4755.d0,-1900.d0,
     &     0.0d0,-6650.d0,-7050.d0,-1500.d0,-1500.d0,-1500.d0,
     &     -1500.d0,-1500.d0,-2000.d0,
     &     -1500.d0,-2100.d0,-1500.d0,-1500.d0,-2700.d0,0.0d0,
     &     -3800.d0,-4000.d0,0.0d0,0.0d0,
     &     -1800.d0,0.0d0,0.0d0,0.0d0,-6100.d0,-4900.d0,
     &     -4500.d0,-1500.d0,-1500.d0,
     &     -2000.d0,0.d0,-1800.d0,120.d0/


      data bkre/2.5d-6,2.0d-5,7.0d9,1.0d10,2.7d7,
     &     8.6d5,1.0d8,0.3d0,0.5d0,
     &     0.13d0,2.0d9,1.0d4,1.5d9,70.d0,2.8d6,7.8d-3,1.5d7,
     &     1.5d6,4.0d8,8.0d5,
     &     4.3d9,6.1d9,2.1d10,1.3d3,4.5d9,1.0d9,3.1d9,
     &     1.4d5,4.5d7,7.3d6,
     &     2.0d8,1.0d8,2.0d10,1.3d9,3.7d-5,6.3d-6,1.0d9,
     &     1.0d10,6.3d3,5.0d5,
     &     4.0d5,2.5d8,1.2d9,1.0d-7,1.0d-5,4.5d9,1.0d9,
     &     1.0d6,1.0d8,2.0d9,
     &     0.1d0,1.6d8,4.6d-6,2.1d5,5.0d0,6.7d3,2.5d9,100.0d0,
     &     6.0d7,1.1d5,1.9d6,
     &     4.0d-4,7.6d5,5.0d7,5.4d-7,2.7d7,4.5d8,2.6d3,
     &     3.5d3,1.9d7,1.0d6,
     &     2.4d4,3.7d5,1.5d9,1.3d6,4.7d0,0.82d0,5.0d3,1.0d7,
     &     4.6d9,4.2d9,3.0d5,
     &     1.0d8,200.d0,1.4d4,2.d8,7.5d7,1.7d7,1.d5,0.31d0,
     &     1.8d-3,1.3d9,5.3d8,
     &     5.0d9,5.0d9,8.0d7,1.2d7,8.8d8,9.1d6,1.7d8,
     &     2.0d8,1.4d6,6.7d-3,
     &     1.9d7,5.0d7,6.0d2,1.0d6,2.5d7,1.0d8,2.0d6,
     &     1.42d2,4.77d3,2.94d2,
     &     3.6d3,1.4d9,3.4d8,
     &     2.5d4,1.0d5,2.5d7,120.d0/


      coefloc = 1.d0/temp-1.d0/298.d0

      do i=1,17
         akeq(i)=bkeq(i)*dexp(dheq(i)*coefloc)
      enddo

      do i=1,21
         akhen(i)=bkhen(i)*dexp(dhhen(i)*coefloc)
      enddo

      do i=1,120
         akre(i)=bkre(i)*dexp(dhre(i)*coefloc)
      enddo

      return
      end
