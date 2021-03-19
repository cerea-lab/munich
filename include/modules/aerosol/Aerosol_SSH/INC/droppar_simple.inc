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
c     UPDATED DROPPAR.INC FOR THE SIZE RESOLVED MODEL
c     LAST UPDATE : 22 November 2005
c*************************************************************************
c     DROPPAR.INC
c*************************************************************************
c     
c     AQUEOUS-PHASE PARAMETERS AND VARIABLES

      INTEGER igso2,igh2o2,ignh3
      INTEGER ighno3,igo3,igco2
C      INTEGER ngch3o2,ngch3o2h,ngch3oh
C      INTEGER ngch3co3h,meqn1,meqn2,kiron
      DOUBLE PRECISION temp_ref
  
      parameter(temp_ref=298.d0)

      INTEGER ngas_aq,ngas_henry

      parameter (ngas_aq = 6)

      DOUBLE PRECISION chenry(ngas_aq),dhhenry(ngas_aq)

      INTEGER nreact_dissoc, nreact_oxydation
      parameter(nreact_dissoc = 9)
      parameter(nreact_oxydation = 5)
  
  
      DOUBLE PRECISION ckdissoc(nreact_dissoc)
      DOUBLE PRECISION dhkdissoc(nreact_dissoc)
      DOUBLE PRECISION ckoxydation(nreact_oxydation)
      DOUBLE PRECISION dhkoxydation(nreact_oxydation)
c     
c     NSP : NUMBER OF TOTAL PARTICLE SECTIONS
c     
c     
c     NBOUNDS : NUMBER OF SECTIONAL BOUNDARIES
c     
c     AQUEOUS-PHASE COMPONENTS
c     
c     
c     GASES IN LOCAL ARRAY
c     
c     INCORRECT  ******* TO BE FIXED *************      


      parameter (igso2  = 1)
      parameter (ignh3  = 2)
      parameter (ighno3 = 3)
      parameter (igh2o2 = 4)
      parameter (igo3   = 5)
      parameter (igco2  = 6)
c     
c     
c     ACTIVATION DIAMETER (Dry)
c     
      double precision dactiv,dsep !,avdiam
      parameter (dactiv = 0.2D0) ! in \mu m

c     
c     SEPARATION DIAMETER FOR TWO SECTIONS (DRY)
c     
      parameter (dsep = 1.5D0)  ! in \mu m
     
c     WET DIAMETERS OF THE TWO SECTIONS

c     If bulk, diameter = avdiam
C      parameter (avdiam = 20)   ! in \mu m     
c     
c     
c     PARAMETER FOR SCALING OF PHOTOLYSIS RATES
c     
C      double precision photo
C      parameter (photo = 1.d0)  
c     
c     fraction of liquid water content that goes to each s.r. section
c     
C      double precision frac1,frac2
C      parameter (frac1 = 0.8d0) ! fraction of lwc in sect. 1
C      parameter (frac2 = 0.2d0) ! fraction of lwc in sect. 2
     
      double precision minso2
      parameter (minso2 = 1.5d-4)

C     LWCFOG: minimum threshold for fog
      
      double precision LWCFOG
      parameter (LWCFOG=0.05d0) ! in g.m-3 

C     Parameters for bimodal activated distribution: 
C     DM: mean diameter; SD: variance.

      double precision DM1_aq, DM2_aq, SD1_aq, SD2_aq
      
      parameter(DM1_aq = 0.4d0)
      parameter(DM2_aq = 2.5d0)
      parameter(SD1_aq = 1.8d0)
      parameter(SD2_aq = 2.15d0)

C      double precision DG_AQ, ACC_AQ
C      parameter(DG_AQ=1.d-5)
C      parameter(ACC_AQ=0.1d0)

C     Henry constants at 298K (in mol L-1 atm-1) 
C     and heat of dissolution (in kcal mol-1)
C     for resp. SO2, NH3, HNO3, H2O2, O3 and CO2
      data chenry /1.23d0, 62.d0, 2.1d5, 7.45d4, 1.13d-2, 3.4d-2/
      data dhhenry /-6.25d0, -8.17d0, -17.27d0, -14.5d0, -5.04d0, 
     &             -4.85d0 /

C     Kinetic constants at 298K (in mol L-1)and delta H/R (in kcal mol-1)
C     For dissociation reactions (resp.):
C     H2O -> H+ + OH-
C     SO2.H2O -> HSO3- + H+
C     HSO3- -> H+ + SO3--
C     H2SO4(aq) -> H+ + HSO4-
C     HSO4- -> H+ + SO4--
C     NH3.H2O -> NH4+ + OH-
C     HNO3 -> NO3- + H+
C     CO2.H2O -> HCO3- + H+
C     HCO3- -> H+ + CO3--
      data ckdissoc /1.d-14, 1.3d-2, 6.6d-8,  1.0d3, 1.02d-2,
     &     1.7d-5, 15.4d0, 4.3d-7, 4.7d-11/
      data dhkdissoc/-6710.d0, 1960.d0, 1500.d0, 0.d0, 2720.d0, 
     &                  -450.d0, 8700.d0, -1000.d0, -1760.d0/ 
C      data dhkdissoc /13.35d0,-4.16d0,-2.23d0,5.4d0,7.91d0,8.65d0,-17.3d0/

C     Kinetic constants at 298K and delta H
C     For sulfur oxydation reactions (resp.):
C     S(IV) + O3 -> S(IV) + O2 (3 constantes)
C     S(IV) + H2O2 -> S(VI) + H2O (2 constantes)
      data ckoxydation / 2.4d4, 3.7d5, 1.5d9, 7.5d7, 13.0d0/
      data dhkoxydation / 0.0d0, 10.99d0, 10.48d0, 8.79d0, 0.0d0/