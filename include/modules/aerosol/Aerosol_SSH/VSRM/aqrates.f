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

      subroutine aqrates(yaq, iph, aqprod, aqdest, yaqprime,ph)
      
C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This routine computes the time derivatives for the aqueous-phase
C     system in the bulk case (one cloud section), the production term 
C     and the destruction rate.
C     If IPH=1, we only compute pH.
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     YAQ: concentration of the aqueous-phase system ([...]).
C     IPH: 0 for all computations, 1 for pH only.
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     
C     -- OUTPUT VARIABLES
C     
C     AQPROD  : production term.
C     AQDEST  : destruction rate.
C     YAQPRIME: first-order derivative defined as 
C     AQPROD-AQDEST*YAQ  ([...]).
C     PH      : pH.
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
C     1) Remove call to QSATURATION.
C     2) Remove initialization to 0 (computed anyhow).
C     3) Remove many debbuging tests.
C     4) Introduce tiny value TINYAQ=1.D-20 from num_aq.inc.
C     5) Minimize calculations in many loops (through local
C     coefficients (COEFLOC).
C     6) Minimize calculations for final computation of production/
C     destruction terms.
C     7) If Sulfate is too small, the final computation is avoided.
C     8) Remove call to DIFFER (net change directly computed).
C     9) Molar mass in data_mass.inc.
C     10)Dimension is equal to 13 (remove NT from the input variables).
C     11)Omit double calculation of VALUES and define default values for
C     QSSA species before the call.
C     12)Remove the lines for computing QSSA species (IRADICAL=0 anyhow).
C     13)Remove many variables not used.
C     14/Add as local data GMOL and WMOL..
C     15)Optimize with coefloc2.
C     16)Remove commons /count/, /istep/, /tnow/, /integration/, /tprof/.
C     17)Remove tlwc,wlm,wl
C     18)Add pres in the arguments of MASS.F
C     19)gmol, wmol no more in dropcom.inc.
C     20)Replace form by formald (F77 function otherwise ??).
C     21)Add pH tricks.
C------------------------------------------------------------------------
C     
C     -- AUTHOR(S)
C     
C     Kathleen Fahey, CEREA, , on the basis of the VSRM model
C     (Carneggie Mellon University).
C     2005/10/3, cleaning and update, Bruno Sportisse, CEREA.
C     
C------------------------------------------------------------------------

      IMPLICIT NONE

      include 'data_mass.inc'
      include 'aerpar.inc'
      include 'droppar.inc'
      include 'aqrates.inc'
      include 'num_aq.inc'
c     
      double precision coefloc,coefloc2
      double precision yaq(13),yaqprime(13)
      double precision aqprod(13), aqdest(13)
      double precision con(28),cmet(4), c(46), gascon(ngas_aq)
      double precision spres(21), gcon(21)
      double precision fgl(21), flg(21), gfgl(21), gflg(21) 
      double precision rp(28), rl(28)
      double precision rr(120)
      double precision ph
      
      double precision hc1, hc2, vlwc, hyd
      double precision arytm, hno2, af
      double precision ah, chen, rad, wvol 
      double precision formald, crustal, salt, cph
      double precision pressure, watvap, watcont
      double precision temp,therm

      double precision wmol(29),gmol(22)

      integer i,iph
      
      common/aqrate/pressure,therm
      common/gascon/gascon      ! in ppm
      common/aqchem1/watcont,watvap,crustal,salt  
                                ! lwc,wat. vapor in g/m3
                                ! crustal species concentration INSIDE DROPLETS (ug/m3)
                                ! Na concentration INSIDE DROPLETS (ug/m3)
      common/sstate/gcon,con,cmet,rad, wvol, cph
!$OMP THREADPRIVATE(/aqrate/, /gascon/, /aqchem1/, /sstate/)

c     
c     MOLECULAR WEIGHTS
c     
      data wmol/81.0d0,  96.0d0,  47.0d0, 62.0d0, 62.0d0,
     &     34.0d0,  48.0d0,  46.0d0, 30.0d0, 46.0d0,
     &     48.0d0,  121.0d0, 76.0d0, 48.0d0, 35.5d0,
     &     17.0d0,  33.0d0,  62.0d0, 18.0d0, 47.0d0,
     &     32.0d0,  35.5d0,  52.5d0, 96.0d0, 112.0d0,
     &     113.0d0, 111.0d0, 60.0d0, 18.0d0/
c     
      data gmol/64.0d0,  98.08d0,  47.02d0, 63.02d0, 44.01d0,
     &     34.02d0, 30.03d0,  46.00d0, 30.01d0, 46.01d0,
     &     48.00d0, 121.05d0, 76.00d0, 48.00d0, 36.50d0,
     &     17.00d0, 33.01d0,  62.01d0, 17.00d0, 47.00d0,
     &     32.00d0, 18.00d0/

      
c     1) Meteorological/thermodynamics conditions  
c     -------------------------------------------     
      temp = therm
      vlwc = watcont*1.d-6      ! vol/vol
      wvol = vlwc

      coefloc2 = 1.d9*wvol
c     radius
      rad = 0.5d0*avdiam*1.d-6  ! from \mu m to m

c     2) Compute the gas-phase concentrations in ppm (SPRES)
c     and in mol/l (GCON)
c     (SOME ARE ASSUMED TO REMAIN CONSTANT DURING THE AQUEOUS-PHASE
c     CHEMICAL PROCESSES, THE REST ARE UPDATED)
c     ALL HCHO AND HCOOH ARE STILL IN THE GAS-PHASE
c     -------------------------------------------------------------    
      coefloc = 8.314d-5*temp/pressure

      spres(1)  = yaq(3)*coefloc/mmSO2 ! SO2 (g)
      spres(2)  = TINYAQ        ! H2SO4 (g)
      spres(3)  = gascon(nghno2) ! HNO2 (g)
      spres(4)  = TINYAQ        ! HNO3 (g) [ALREADY DISSOLVED]
      spres(5)  = 350.d0        ! CO2 (g)
      spres(6)  = yaq(4)*coefloc/mmH2O2 ! H2O2 (g)

      hc1=8.314d-2*temp*vlwc*akhen(7)
      hc2=coefloc/mmHCHO
      spres(7)  = yaq(1)*hc2/(hc1+1.d0) ! HCHO (g)
      spres(8)  = yaq(2)*coefloc/mmHCOOH ! HCOOH (g)
      spres(9)  = gascon(ngno)  ! NO (g)
      spres(10) = gascon(ngno2) ! NO2 (g)
      spres(11) = gascon(ngo3)  ! O3 (g)
      spres(12) = gascon(ngpan) ! PAN (g)
      spres(13) = TINYAQ        ! CH3COOOH (g)
      spres(14) = TINYAQ        ! CH3OOH (g)
      spres(15) = TINYAQ        ! HCl (g)  [ALREADY DISSOLVED]
      spres(16) = gascon(ngoh)  ! OH (g)
      spres(17) = gascon(ngho2) ! HO2 (g)
      spres(18) = gascon(ngno3) ! NO3 (g)
      spres(19) = yaq(5)*coefloc/mmNH3 ! NH3 (g)
      spres(20) = gascon(ngch3o2) ! CH3O2 (g)
      spres(21) = TINYAQ        ! CH3OH (g)
c     spres(22) = watvap*8.314d-2*temp/(1000.d0*18.d0*pressure) ! H2O (g)

      coefloc = 0.08206d6*temp
      do i=1,21
         gcon(i) = spres(i)/coefloc
      enddo

c     3) Compute metal species on the basis of the crustal aerosol mass
c     NA+ is an input.
c     -----------------------------------------------------------------
      coefloc = 1.d-3/watcont

      cmet(1) = firon*crustal*coefloc/55.85d0 ! Fe(3+) mol/l
      cmet(2) = fman*crustal*coefloc/54.9d0 ! Mn(2+) mol/l
      cmet(3) = salt*coefloc/23.d0 ! Na(+)  mol/l
      cmet(4) = caratio*crustal*coefloc/40.08d0 ! Ca(2+) mol/l
c     
c     DO NOT LET THE Fe(3+) AND Mn(2+) CONCENTRATIONS EXCEED A
c     CERTAIN LIMIT BECAUSE THEY CAUSE EXTREME STIFFNESS DUE TO
c     THE REACTION S(IV)->S(VI) 
c     
c     if (cmet(1) .gt. 1.0e-5) cmet(1)=1.0e-5
c     if (cmet(2) .gt. 1.0e-5) cmet(2)=1.0e-5
c     
c     4) Compute the aqueous concentrations ([M])
c     -------------------------------------------
      
      con(1) = yaq(6)*coefloc/wmol(1) ! S(IV)
      con(2) = yaq(11)*coefloc/wmol(2) ! S(VI)
      con(3) = 0.d0             ! N(III) (LATER)
      con(4) = yaq(8)*coefloc/wmol(4) ! N(V)
      con(5) = 0.d0             ! CO2 (LATER)
      con(6) = yaq(7)*coefloc/wmol(6) ! H2O2
      con(7) = akhen(7)*spres(7)*1.d-6 ! HCHO
      con(8) = 0.d0             ! HCOOH (LATER)
      con(9) = 1.0d-6*akhen(9)*spres(9) ! NO
      con(10) = 1.0d-6*akhen(10)*spres(10) ! NO2
      con(11) = 1.0d-6*akhen(11)*spres(11) ! O3
      con(12) = 1.0d-6*akhen(12)*spres(12) ! PAN
      con(13) = 1.0d-6*akhen(13)*spres(13) ! CH3COOOH
      con(14) = 1.0d-6*akhen(14)*spres(14) ! CH3OOH
      con(15) = yaq(9)*coefloc/wmol(15) ! HCl
      con(16) = 0.d0            ! OH (LATER)
      con(17) = 0.d0            ! HO2 (LATER)
      con(18) = 0.d0            ! NO3 (LATER)
      con(19) = yaq(10)*coefloc/wmol(19) ! NH3
      con(20) = 1.0d-6*akhen(20)*spres(20) ! CH3O2
      con(21) = 1.0d-6*akhen(21)*spres(21) ! CH3OH
      con(22) = 0.d0            ! Cl (LATER)
      con(23) = 0.d0            ! ClOH- (LATER)
      con(24) = 0.d0            ! SO4- (LATER)
      con(25) = 0.d0            ! SO5- (LATER)
      con(26) = yaq(12)*coefloc/wmol(26) ! HSO5-
      con(27) = yaq(13)*coefloc/wmol(27) ! HOCH2SO3-
      con(28) = 0.d0            ! CO3- (LATER)

      do i=1, 28
         con(i) = dmax1(con(i), TINYAQ)
      enddo

c     5) CALCULATION OF pH and VOLATILE CONCENTRATIONS (CO2, N(III), HCOOH)
c     ---------------------------------------------------------------------
      
      call fullequil(con,spres,cmet,akeq,akhen,vlwc,temp,hyd) 
      ph=-DLOG10(hyd)

C     Stop if only pH required
C     
      IF (IPH.EQ.0) THEN
c     N(III) (HNO2)     
         ah = 8.314d-2*temp*vlwc*akhen(3)*(1.d0+akeq(7)/hyd)/pressure
         hno2=spres(3)/(1.d0+ah)
         con(3)=akhen(3)*(1.d0+akeq(7)/hyd)*1.d-6*hno2

c     CO2 (Tot)aq ([M])    
         chen=akhen(5)*(1.d0+akeq(8)/hyd+akeq(8)*akeq(9)/hyd**2d0)
         con(5)=chen*spres(5)*1.d-6 

c     HCOOH     
         af=8.314d-2*temp*vlwc*akhen(8)*(1.d0+akeq(13)/hyd)/pressure
         formald=spres(8)/(1.d0+af) ! NEW HCOOH(g) ppm
         con(8)=akhen(8)*(1.d0+akeq(13)/hyd)*1.d-6*formald

c     6) Default values for radical species
c     -------------------------------------

         con(16)=1.d-5*TINYAQ   !OH(g)
         con(17)=con(16)        !HO2(g)
         con(18)=con(16)        !NO3(g)
         con(23)=con(16)        !ClOH- 
         con(24)=con(16)        !SO4-
         con(25)=con(16)        !SO5-
         con(28)=con(16)        !CO3-

c     7) Compute ionic species
c     ------------------------
         
         call values(hyd, con, akeq, c)

c     8) Compute reaction rates
c     -------------------------

         call react(c,cmet,con,akre,rr,arytm)

c     9) Compute production/consumption for aqueous-phase chemistry
c     -------------------------------------------------------------

         call addit(rr, arytm, rp, rl)

c     10) Compute mass transfer rates
c     -------------------------------
         
         call mass(wvol,rad,temp,pressure,
     &        gcon,con,c,akeq,akhen,fgl,flg,gfgl,gflg)

c     11) Compute first derivative, production/consumption term     
c     ---------------------------------------------------------

C     SULFATE SPECIES 
         if (gascon(ngso2) .le. minso2) then

            aqprod(3) = 0.d0
            aqdest(3) = 0.d0
            aqprod(6) = 0.d0
            aqdest(6) = 0.d0
            aqprod(11) = 0.d0
            aqdest(11) = 0.d0
            aqprod(12) = 0.d0
            aqdest(12) = 0.d0
            aqprod(13) = 0.d0
            aqdest(13) = 0.d0

         else

c     ** GAS-PHASE SPECIES **
c     (RATES IN UG/M3 AIR S)

C     SO2(g)
            coefloc = 1.d9*gmol(1)
            aqprod(3) = coefloc*gflg(1)
            aqdest(3) = coefloc*gfgl(1)
c     
c     ** AQUEOUS-PHASE SPECIES **
c     

C     S(IV)
            coefloc = coefloc2*wmol(1)
            aqprod(6) = coefloc*(rp(1)+fgl(1))
            aqdest(6) = coefloc*(rl(1)+flg(1))

C     S(VI)
            coefloc = coefloc2*wmol(2)
            aqprod(11) = coefloc*(rp(2)+fgl(2))
            aqdest(11) = coefloc*(rl(2)+flg(2))
            
C     HSO5-
            coefloc = coefloc2*wmol(26)
            aqprod(12) = coefloc*rp(26)
            aqdest(12) = coefloc*rl(26)

C     HMSA
            coefloc = coefloc2*wmol(27)
            aqprod(13) = coefloc*rp(27)
            aqdest(13) = coefloc*rl(27)

         endif

c     ** GAS-PHASE SPECIES **
c     (RATES IN UG/M3 AIR S)

C     HCHO (TOT)
         coefloc = coefloc2*wmol(7)
         aqprod(1) = coefloc*rp(7)
         aqdest(1) = coefloc*rl(7)

C     HCOOH (TOT)
         coefloc = coefloc2*wmol(8)
         aqprod(2) = coefloc*rp(8)
         aqdest(2) = coefloc*rl(8)

C     H2O2(g)
         coefloc = 1.d9*gmol(6)
         aqprod(4) = coefloc*gflg(6)
         aqdest(4) = coefloc*gfgl(6)

C     NH3(g)
         coefloc = 1.d9*gmol(19)
         aqprod(5) = coefloc*gflg(19)
         aqdest(5) = coefloc*gfgl(19)
         
c     
c     ** AQUEOUS-PHASE SPECIES **
c     
C     H2O2
         coefloc = coefloc2*wmol(6)
         aqprod(7) = coefloc*(rp(6)+fgl(6))
         aqdest(7) = coefloc*(rl(6)+flg(6))

C     N(V)
         coefloc = coefloc2*wmol(4)
         aqprod(8) = coefloc*rp(4)
         aqdest(8) = coefloc*rl(4)
         
C     Cl-
         coefloc = coefloc2*wmol(15)
         aqprod(9) = coefloc*(rp(15)+fgl(15))
         aqdest(9) = coefloc*(rl(15)+flg(15))

C     NH4+
         coefloc = coefloc2*wmol(19)
         aqprod(10) = coefloc*(rp(19)+fgl(19))
         aqdest(10) = coefloc*(rl(19)+flg(19))
c     
c     CALCULATION OF DERIVATIVES
C     AND OF THE APPROPRIATE DESTRUCTION RATE 
C     (in place of the destruction flux).
c     
         do i=1, 13
            yaqprime(i) =  aqprod(i)-aqdest(i) 
            if (yaq(i) .le.TINYAQ) then
               aqdest(i) = TINYAQ*100.d0
            else 
               aqdest(i) = aqdest(i)/yaq(i)
               aqdest(i) = dmax1(aqdest(i), TINYAQ*100.d0)
            endif
         enddo

      ENDIF

      end
