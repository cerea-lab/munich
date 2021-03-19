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

      subroutine aqoperator(NS,tbeg,deltat,gas,aerosol,lwc,
     &     temp,press,
     &     fHCHO,fhcooh,fso2,fh2o2,fNH3,fhno3,fhcl,
     &     ifirstact,
     &     pH,
     &     aqso2,aqh2o2,
     &     fdist)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This routine computes evolution of gas and aerosol species for one 
C     timestep.
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     TBEG   : initial time            [min].
C     DELTAT : timestep                [min].
C     LWC    : liquid water content    [g.m^{-3}].
C     TEMP   : temperature             [K].
C     PRESS  : pressure                [atm].
C     F*     : conversion factor \mu.m^3 to ppm.
C     FDIST* : fraction of the activated bimodal distribution that go to
C     the initial distribution.
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     GAS    : gas-phase concentration [ppm].
C     AEROSOL: aerosol concentration   [\mug.m^-3].
C     PH     : pH.
C     
C     -- OUTPUT VARIABLES
C     
C     aqSO2  : aqueous part of SO2  [\mug.m^-3]
C     aqH2O2 : aqueous part of H2O2 [\mug.m^-3]
C     
C------------------------------------------------------------------------
C     
C     -- REMARKS
C     
C     The vector YAQ is defined in the following way (in ug/m3):
C     
C     YAQ(1) : total HCHO 
C     YAQ(2) : total HCOOH 
C     YAQ(3) : total SO2 
C     YAQ(4) : total H2O2
C     YAQ(5) : NH3(g)   
C     YAQ(6) : SO2 (aq)
C     YAQ(7) : H2O2 (aq)
C     YAQ(8) : NITRATE (aq)
C     YAQ(9) : CHLORIDE (aq)
C     YAQ(10): AMMONIUM (aq)
C     YAQ(11): SULFATE (aq) 
C     YAQ(12): HSO5-
C     YAQ(13): HMSA 
C     
C------------------------------------------------------------------------
C     
C     -- MODIFICATIONS
C     
C     1) Optimize conversion (with data_mass.inc).
C     2) Remove variables SBEF, SAFT, SAFT1 (not used).
C     3) Optimize differences to be added.
C     4) Define TINYAQ2 in num_aq.inc.
C     5) Rename P in PRESS. 
C     6) Rename NGCN abd NGCC to NGN and NGC.
C     7) Replace DD by DAER.
C     8) Remove include 'dropcom.inc'.
C     9) Remove the initialization (in VSRMCHEM.f now) and IAQ as
C     input.
C     10) Remove the N2O5 trick.
C     11) Conversion factors f* as input parameters (not computed).
C     12) Replace NSECT by NS.
C     13) Remove common /integration/
C     14) Replace DAER by DSF.
C     15) Remove include 'dropcom.inc'.
C     16) Use ifirstact.
C     17) 2005/11/25: Add aqSO2 and aqH2O2 in output (Edouard Debry). 
C------------------------------------------------------------------------
C     
C     -- AUTHOR(S)
C     
C     Kathleen Fahey, CEREA, , on the basis of the VSRM model
C     (Carneggie Mellon University).
C     2005/10/3, cleaning and update, Bruno Sportisse, CEREA.
C     
C------------------------------------------------------------------------
c     
      IMPLICIT NONE
c     
      include 'data_mass.inc'
      include 'aerpar.inc'
      include 'droppar.inc'
      include 'num_aq.inc'

      INTEGER NS

      double precision gas(ngas_aq), aerosol(NS,naers)
c      double precision sodium(NS)
      double precision deltat, lwc, temp
      double precision gascon(ngas_aq)
      double precision yaq(meqn1)
      double precision watcont, watvap
      double precision tbeg,press,pressure, therm

      double precision fHCHO,fhcooh,fso2,fh2o2,fNH3,fhno3,fhcl
      double precision hno3, hcl
      double precision tnitold,chlorold,ammonold
      double precision sulfateold,hso5old,hmsaold
      double precision crustal,salt,stime,stout

      double precision aqprod(13),aqdest(13),yaqprime(13),pH
      double precision aqso2,aqh2o2
      double precision fdist(NS)

      integer isp,i,isect
      integer ifirstact
      
      common / gascon / gascon
      common / aqrate / pressure, therm
      common/aqchem1/watcont,watvap,crustal,salt
!$OMP THREADPRIVATE(/gascon/, /aqrate/, /aqchem1/)
      
      DO i=1,13
         YAQ(i)=0.D0
      ENDDO

C     1) Compute kinetic rates and TRANSFER 
C     the input data VIA COMMONs TO AQRATE
C     ------------------------------------------
      pressure = press
      call constants(temp)

      watcont = lwc
      therm = temp
      do i=1, ngas_aq
         gascon(i) = gas(i)
      enddo

c     2) UNIT CONVERSION from ppm to ug/m3
c     LOADING OF THE CONCENTRATIONS IN THE YAQ VECTOR
c     -----------------------------------------------
      
c     THE TOTAL FORMALDEHYDE/FORMIC ACID ARE TRANSFERRED AS GAS-PHASE SPECIES
c     (THEY ARE NOT INCLUDED IN THE AQUEOUS OR AEROSOL MATRICES)

      yaq(1) = gas(nghcho) /fhcho !total HCHO (ug/m3)
      yaq(2) = gas(nghcooh)/fhcooh !total HCOOH (ug/m3)

c     GAS-PHASE SPECIES (in ug/m3)

      yaq(3) = gas(ngso2) /fso2 !total SO2 
      yaq(4) = gas(ngh2o2)/fh2o2 !total H2O2
      yaq(5) = gas(nga)   /fNH3 !NH3(g) 

c     AQUEOUS-PHASE SPECIES (in ug/m3)
c     THE DROPLETS ARE ASSUMED TO BE WITHOUT SO2 OR H2O2 
C     IN THE BEGINNING OF A TIMESTEP.

      yaq(6) = 1.d-10           !SO2 (aq) 
      yaq(7) = 1.d-10           !H2O2 (aq)

C     Projection for activated aerosols

      do isect=ifirstact, NS
         yaq(8) = yaq(8) + aerosol(isect,nan) !NITRATE (aq)
         yaq(9) = yaq(9) + aerosol(isect,nac) !CHLORIDE (aq)
         yaq(10) = yaq(10) + aerosol(isect,naa) !AMMONIUM (aq)
         yaq(11) = yaq(11) + aerosol(isect,na4) !SULFATE (aq) 
         yaq(12) = yaq(12) + aerosol(isect,nahso5) !HSO5-
         yaq(13) = yaq(13) + aerosol(isect,nahmsa) !HMSA 
      enddo

c     CALCULATION OF THE DISSOLUTION OF THE AVAILABLE HNO3 AND HCl
c     TO THE DROPLETS/PARTICLES. WE ASSUME THAT THE TIMESCALE OF DISSOLUTION
c     IS SMALL AND THAT ALL HNO3 AND HCl ARE DISSOLVED (ZERO VAPOR PRESSURE)

      hno3 = gas(ngn)/fhno3     ! HNO3(g)
      hcl =  gas(ngc)/fhcl      ! HCl(g)

c     3) SAVE THE OLD AQUEOUS-PHASE CONCENTRATIONS BEFORE THE INTEGRATION
c--------------------------------------------------------------------
      
      tnitold   = yaq(8)        !NITRATE (aq)
      chlorold  = yaq(9)        !CHLORIDE (aq)
      ammonold  = yaq(10)       !AMMONIUM (aq)
      sulfateold= yaq(11)       !SULFATE (aq) 
      hso5old   = yaq(12)       !HSO5-
      hmsaold   = yaq(13)       !HMSA 

c     4) CALCULATION OF SODIUM VECTOR (NEEDED FOR THE INTEGRATION ROUTINE)
c     (TRANSFERRED VIA COMMON STATEMENT AQRATESA)
c     --------------------------------------------------------------------
      
c      do isect=1, NS
c         sodium(isect) = aerosol(isect,nas)
c      enddo

c     5) CALCULATE CRUSTAL SPECIES CONCENTRATION INSIDE DROPLETS
c     IT IS USED IN THE AQUEOUS-PHASE CHEMISTRY CALCULATION TO
c     ESTIMATE Fe and Mn CONCENTRATIONS
c     ----------------------------------------------------------
      
      crustal = 0.d0
      salt =0.d0
      do isect=ifirstact,NS
         crustal=crustal+aerosol(isect,nar)
         salt=salt+aerosol(isect,nas)
      enddo   
      
C     6) ADD THE GAS-PHASE CONCENTRATIONS TO THE TOTAL FOR HCl AND HNO3
C     -----------------------------------------------------------------

      gas(ngn) = 0.d0           ! all HNO3 is dissolved
      gas(ngc) = 0.d0           ! all HCl is dissolved
      yaq(8) = yaq(8)+hno3      ! HNO3 increase
      yaq(9) = yaq(9)+hcl       ! HCl increase  

c     7) Perform integration
c     ----------------------

      stime = tbeg*60.d0        ! in seconds
      stout = (tbeg+deltat) * 60.d0 ! in seconds
      
      call aqintegr(yaq, stime, stout)

C     7bis) Compute PH
C     ----------------
      CALL aqrates(yaq, 1,aqprod, aqdest, yaqprime,ph)

C     8) ADD THE MASS CHANGE TO THE CORRESPONDING SECTION 
c     ---------------------------------------------------
      
      do isect=1, NS
         aerosol(isect,nan)    = aerosol(isect,nan)
     &        + fdist(isect) * (yaq(8) - tnitold) ! NITRATE (aq) 
         aerosol(isect,nac)    = aerosol(isect,nac)
     &        + fdist(isect) * (yaq(9) - chlorold) ! CHLORIDE (aq)
         aerosol(isect,naa)    = aerosol(isect,naa) 
     &        + fdist(isect) * (yaq(10) - ammonold) ! AMMONIUM (aq)
         aerosol(isect,na4)    = aerosol(isect,na4)
     &        + fdist(isect) * (yaq(11) - sulfateold) ! SULFATE (aq)
         aerosol(isect,nahso5) = aerosol(isect,nahso5)
     &        + fdist(isect) * (yaq(12) - hso5old) ! HSO5 (aq)
         aerosol(isect,nahmsa) = aerosol(isect,nahmsa)
     &        + fdist(isect) * (yaq(13) - hmsaold) ! HMSA (aq)
      enddo

      do isect=1, NS
         do isp=1, naers     
            aerosol(isect,isp)=dmax1(aerosol(isect,isp),tinyaq2)
         enddo 
      enddo 

c     9) Update gas and aerosols 
c     --------------------------  
      aqSO2 = yaq(6)
      aqH2O2= yaq(7)

c     ** GAS-PHASE SPECIES **     
c     IMPORTANT : THE DISSOLVED S(IV) and H2O2 ARE TRANSFERRED
c     BACK TO THE GAS PHASE
c     
C     0.79 = MSO2/MHSO3- = 64.06/81.07
      
      gas(nghcho) = yaq(1)         *fhcho !total HCHO (ppm)
      gas(nghcooh)= yaq(2)         *fhcooh !total HCOOH (ppm)
      gas(ngso2)  = (yaq(3)+yaq(6)*0.7901d0)
     &     *fso2                !SO2(g) (ppm)
      gas(ngh2o2) = (yaq(4)+yaq(7))*fh2o2 !H2O2(g) (ppm)              
      gas(nga)    = yaq(5)         *fNH3 !NH3(g) (ppm)

      return
      end
