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

      subroutine vsrmchem(NGAS,NAER,NS,NB,NICD,RHO_AERO,
     &     fixed_rho_aero,DBF_AERO,dsf_aero,xbf_aero,C_GAS,
     &     C_AER,HUMID,press2,temp,lwc_c,t0,t1,rain_rate,pH,
     &     qscav_gas,qscav_aer,qscav_num,INUM,C_NUM,dqlimit,
     &     CLD_SP_INT,ID_SIZE,List_species, n_isorropia,n_aec, 
     &     section_pass,IREDIST,with_fixed_density)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This routine solves the aqueous-phase model.
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C
C     HUMID : specific humidity    ([kg.kg^{-1}]).
C     PRESS2: pressure             ([Pa]).
C     TEMP  : temperate            ([T]).
C     LWC_C : liquid water content ([g.m^{-3}]).
C     T0/T1 : initial/final time   ([s]).
C     rain_rate : rain precipitation rate ([mm/hr]).
C     CLD_SP_INT: CLD_SP_INT
C     ID_SIZE: relations between bin index and size sections
C     NS: number of bins
C     NGAS: number of gas species
C     NAER: number of aerosols species
C     NICD: number of cloud interact species
C     NB: number of size section
C
C     -- INPUT/OUTPUT VARIABLES
C
C     C_GAS: gas-phase concentration ([\mu.g/m^3]).
C     C_AER: aerosol concentration ([\mu.g/m^3]).
C     C_NUM: number concentration ([#/m^3]).

C     -- OUTPUT VARIABLES
C
C     PH : pH.
C     qscav_gas : scavenged gas quantity by in-cloud scavenging
C     qscav_aer : scavenged aerosol quantity by in-cloud scavenging
C     qscav_num : scavenged number quantity by in-cloud scavenging
C     for dissolved aerosols species ([\mu.g/m^3])
C     and dissolves gas      species ([ppm])
C     
C------------------------------------------------------------------------
C     
C     -- REMARKS
C     
C------------------------------------------------------------------------
C     
C     -- MODIFICATIONS
C     
C     1) Optimize conversion and define molar mass as parameters.
C     2) Introduce pointers for aerosols (E*) and for gas (ICTM*).
C     3) Clipping at the end (TINYAQ).
C     4) Rewrite mass balance for S and N.
C     5) Add inclusion of data_mass.inc and pointer_ctm.inc.
C     6) Remove call to VSRM (now included as such).
C     7) Define NITSUBAQ (number of subcycling timesteps = 10 here).
C     8) Remove include 'dropcom.inc'.
C     9) Move the initialization (before in AQOPERATOR.f).
C     10)Remove IAQ.
C     11)Clarify the computation of N2O5 and update the computation
C     of nitrate mass.
C     12)Remove the GOTO (with NITVSRM possible restart).
C     13)Replace NSECT by NS.
C     14)Remove tfin.
C     15)Conversion factors f* as input parameters (not computed) for
C     AQOPERATOR.F.
C     16)Remove call to DROPINIT.F (included as such).
C     17)Remove computations of DAER (replaced by DSF).
C     18)Update computation of total mass.
C     19)Add 3bis to avoid numerical difficulties for low SO2 (before:
C     in DECISIONS.F).
C     20)Comment call to DECISIONS.F
C     21)Uniform distribution of N2O5 in bins.
C     22)Remove include 'dropcom.inc'.
C     23) 2005/11/25: Implement in-cloud scavenging (Edouard Debry).
C     24) Shupeng ZHU: add composition support for SCRAM (Oct.2014)
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

      include 'param.inc'
      include 'CONST_A.INC'
      include 'data_mass.inc'
      include 'aerpar.inc'
      include 'droppar.inc'
      include 'varp.inc'
      include 'varp_cloud.inc'
      include 'num_aq.inc'
C       include 'vara.inc'
C       include 'dynaero.inc'

      INTEGER NGAS,NAER,NS,NB,NF,NICD
      INTEGER k,l,IREDIST
      INTEGER nesp,f,jesp,b
      INTEGER section_pass      
      INTEGER CLD_SP_INT(NICD)
      INTEGER ID_SIZE(NS,2)
      INTEGER INUM
      INTEGER List_species(NAER)
      INTEGER n_isorropia,n_aec

      INTEGER ifirstact!first activated bin
      INTEGER isect,isp,i,istep,nistep,ind_ok,j,jj,jb
      
      double precision d(NB)
      double precision mass_rdb(NB,naers)
      double precision totQ(NB)
      double precision Nub(NB)
      double precision dqlimit
      double precision dold(NS)
      double precision volume(NS)      
      
      double precision coefloc,tinit,coefloc2
      double precision dh2SO4,dso2
      double precision gas(ngas_aq), aerosol(NS,naers)   
      double precision gasav(ngas_aq), aerosav(NS,naers)
      double precision rh,temp,press,lwc_c,t0,t1,initso2
      double precision dnew(NS),deltat,press2
      double precision numberconc(NS),totmass(NS)
      double precision C_GAS(NGAS)
      double precision C_AER(NS,NAER)
      double precision C_NUM(NS)      
      double precision HUMID
      double precision nitbef,nitaf
      double precision sulfbef,sulfaf,sbal
      double precision pH
      double precision qscav_gas(NGAS)
      double precision qscav_aer(NS,NAER)
      double precision qscav_num(NS)
      double precision rain_rate
      double precision fdist(NS), fdist2(NS)
      double precision foa(n_aec, NS), ftot

      double precision collision_eff,drop_diam,scav_coef
      double precision scav_frac,aqSO2,aqH2O2
      double precision DBF_AERO(NB+1),RHO_AERO(NS)
      double precision dsf_aero(NS),xbf_aero(NB+1)
      double precision fixed_rho_aero,totmasstot
      double precision fixed_diameter(NS)

      double precision totmass_init(NS)
      double precision fNH3,fHNO3,fHCl,fSO2
      double precision fHCHO,fO3,fOH,fHCOOH
      double precision fNO,fNO3,fNO2,fPAN,fHO2,fH2O2,fHNO2      

      data collision_eff /0.9d0/

      double precision liquid_density(naers)
      data liquid_density /0.97D-06, 0.85D-06, 
     &     0.91D-06, 1.50D-06, 1.15D-06, 1.84D-06, 1.00D-06,
     &     2.25D-06, 1.30D-06, 2.33D-06, 1.84D-06, 1.835D-06/

      integer with_fixed_density ! YK

c     1) Initialisation
c     ------------------
c   Initialize system pointer
      ictmNH3=CLD_SP_INT(1)
      ictmHNO3=CLD_SP_INT(2)
      ictmHCl=CLD_SP_INT(3)
      ictmSO2=CLD_SP_INT(4)
      ictmH2O2=CLD_SP_INT(5)
      ictmHCHO=CLD_SP_INT(6)
      ictmHNO2=CLD_SP_INT(7)
      ictmO3=CLD_SP_INT(8)
      ictmOH=CLD_SP_INT(9)
      ictmHO2=CLD_SP_INT(10)
      ictmNO3=CLD_SP_INT(11)
      ictmNO=CLD_SP_INT(12)
      ictmNO2=CLD_SP_INT(13)
      ictmPAN=CLD_SP_INT(14)
      ictmH2SO4=CLD_SP_INT(15)

      nesp_isorropia=n_isorropia
      nesp_aec=n_aec

      EMD=List_species(1)
      EBC=List_species(2)
      ENa=List_species(3)
      ESO4=List_species(4)
      ENH3=List_species(5)
      ENO3=List_species(6)
      ECl=List_species(7)
      EH2O=NAER
      E1=1
      E2=NAER-1

      do jesp=1,n_isorropia
	isorropia_species(jesp) = List_species(2+jesp)
      enddo

      do jesp=1,n_aec
	aec_species(jesp) = List_species(2+n_isorropia+jesp)
      enddo

      call initactiv(NS,dsf_aero,ifirstact)

      nistep = NITSUBAQ

      do i=1,ngas_aq
         gas(i)=0.d0
      enddo
      
      do isect=1,NS
         do isp=1,naers
            aerosol(isect,isp)=0.d0
         enddo
      enddo

      do isp=1,NGAS!SZ
         qscav_gas(isp)=0.d0
      enddo

      do j=1,NS!SZ
        do isp=1,NAER!SZ
         qscav_aer(j,isp)=0.d0
        enddo
      enddo

      do isp=1,NS!SZ
         qscav_num(isp)=0.d0
      enddo      

      deltat = (t1-t0)/60.d0/nistep !timestep in min

      press  = press2/101325.d0 ! pressure in atm

C     2) Compute sectional diameters and bimodal distribution
C     -------------------------------------------------------

C      call newdist(NS,DBF_AERO,dsf_aero,fdist,fdist2,ifirstact)

C     3) Conversion from \mug/m^3 to ppm
C     -----------------------------------

      coefloc = 8.314d-5*temp/press

      fNH3 = coefloc/mmNH3      !NH3
      fHNO3= coefloc/mmHNO3     !HNO3
      fHCl = coefloc/mmHCl      !HCl
      fSO2 = coefloc/mmSO2      !SO2
      fH2O2= coefloc/mmH2O2     !H2O2
      fHCHO= coefloc/mmHCHO     !HCHO
      fHNO2= coefloc/mmHNO2     !HNO2
      fO3  = coefloc/mmO3       !O3
      fOH  = coefloc/mmOH       !OH
      fHO2 = coefloc/mmHO2      !HO2
      fNO3 = coefloc/mmNO3      !NO3
      fNO  = coefloc/mmNO       !NO
      fNO2 = coefloc/mmNO2      !NO2
      fPAN = coefloc/mmPAN      !PAN
      fHCOOH = coefloc/mmHCOOH  !HCOOH

      coefloc2 = 8.314d0*temp/press2

      gas(nga)    = C_GAS(ictmNH3)  * coefloc2 / mmNH3  !NH3
      gas(ngn)    = C_GAS(ictmHNO3) * coefloc2 / mmHNO3 !HNO3
      gas(ngc)    = C_GAS(ictmHCl)  * coefloc2 / mmHCl  !HCl
      gas(ngso2)  = C_GAS(ictmSO2)  * coefloc2 / mmSO2  !SO2
      gas(ngh2o2) = C_GAS(ictmH2O2) * coefloc2 / mmH2O2 !H2O2
      gas(nghcho) = C_GAS(ictmHCHO) * coefloc2 / mmHCHO !HCHO
      gas(nghno2) = C_GAS(ictmHNO2) * coefloc2 / mmHNO2 !HNO2
      gas(ngo3)   = C_GAS(ictmO3)   * coefloc2 / mmO3   !O3
      gas(ngoh)   = C_GAS(ictmOH)   * coefloc2 / mmOH   !OH
      gas(ngho2)  = C_GAS(ictmHO2)  * coefloc2 / mmHO2  !HO2
      gas(ngno3)  = C_GAS(ictmNO3)  * coefloc2 / mmNO3  !NO3
      gas(ngno)   = C_GAS(ictmNO)   * coefloc2 / mmNO   !NO
      gas(ngno2)  = C_GAS(ictmNO2)  * coefloc2 / mmNO2  !NO2
      gas(ngpan)  = C_GAS(ictmPAN)  * coefloc2 / mmPAN  !PAN

      gas(nghcooh)   = 0.1d0*gas(ngh2o2) !HCOOH
      gas(ngch3o2h)  = 0.2d0*gas(ngh2o2) !CH3OOH(g)  ppm = 0.2*H2O2
      gas(ngch3o2)   = 1.0d-6            !CH3O2(g)   ppm
      gas(ngch3oh)   = 1.0d-3            !CH3OH(g)   ppm = 1 ppb
      gas(ngch3co3h) = 0.05d0*gas(ngh2o2)!CH3C(O)OOH ppm = 0.05*H2O2

      initso2 = gas(ngso2)


      DO i = 1,NS
         aerosol(i,nas) = C_AER(i,ENa) ! Na
         aerosol(i,na4) = C_AER(i,ESO4) ! SO4
         aerosol(i,naa) = C_AER(i,ENH3) ! NH4
         aerosol(i,nan) = C_AER(i,ENO3) ! NO3
         aerosol(i,nac) = C_AER(i,ECl) ! Cl
         aerosol(i,nar) = C_AER(i,EMD) ! DUST
         aerosol(i,nae) = C_AER(i,EBC) ! EC
         aerosol(i,nao) = 0.0                ! POA
         ftot=0.0
         DO j = 1,nesp_aec
            jesp = aec_species(j)
            aerosol(i,nao) = aerosol(i,nao) + C_AER(i,jesp)
            foa(j,i) = C_AER(i,jesp)
            ftot = ftot + C_AER(i,jesp)
         ENDDO
         DO j=1, nesp_aec
            if (ftot .GT. 0.d0) then
            foa(j,i)=foa(j,i)/ftot
            else
               foa(j,i)=0.d0
            endif
         ENDDO
         aerosol(i,naw) = C_AER(i,EH2O) ! H2O
      ENDDO      


C     Calculation of aerosol number and diameter of each section before aqueous chemistry
C     --------------------------------

      DO i = 1, NB
         fixed_diameter(i) = sqrt(DBF_AERO(i)*DBF_AERO(i+1))
      ENDDO


      do i = 1,NS
         b=ID_SIZE(i,1)
         totmass_init(i)= aerosol(i,naa)+aerosol(i,na4)+
     &        aerosol(i,nan)+aerosol(i,nac)+aerosol(i,nas)+
     &        aerosol(i,nao)+aerosol(i,nae)+aerosol(i,nar)

         if(INUM.EQ.1) then     ! Case when the number concentration is followed
            ! Real value
            numberconc(i) = C_NUM(i)
            ! Calcul of diameter for redistribution euler
            if (numberconc(i) .GT. 0.d0) then
               dold(i) =(totmass_init(i)/
     &              (cst_pi6*RHO_AERO(i)*numberconc(i)))
     &              **cst_FRAC3
            else
               dold(i)=sqrt(DBF_AERO(b)*DBF_AERO(b+1))
            endif

         else
            ! Calculated value
            numberconc(i)=totmass_init(i)/(dsf_aero(i)**3.d0)
     &           /cst_pi6/fixed_rho_aero
            dold(i)=sqrt(DBF_AERO(b)*DBF_AERO(b+1))

         endif

      enddo


C     3bis) At low SO2 concentration, transfer all SO2 to
C     SO4-- in order to avoid numerical difficulties
C     ---------------------------------------------------
      IF (gas(ngso2).LE.minso2) THEN
         IF (gas(ngh2o2).GE.gas(ngso2)) THEN
            gas(ngh2o2) = gas(ngh2o2) - gas(ngso2)
         ELSE
            gas(ngh2o2) = 0.d0
         ENDIF
         dso2 = (C_GAS(ictmSO2)/(NS-IFIRSTACT+1))*mmSO4/mmSO2
         gas(ngso2) = 0.d0

         DO i=IFIRSTACT,NS
            aerosol(i,na4) = aerosol(i,na4)+dso2
         ENDDO
      ENDIF

C     4) Transfer all H2SO4(g) to aerosol phase
C     ------------------------------------------

      dh2SO4 = (C_GAS(ictmH2SO4)/NS)*mmSO4/mmH2SO4
      DO i=1,NS
         aerosol(i,na4) = aerosol(i,na4)+dh2SO4
      ENDDO
      C_GAS(ictmH2SO4) = 0.d0

C     5) Save initial concentrations (for restart)
C     ----------------------------------------------

      do isect=1,NS
         do isp=1,naers
            aerosav(isect,isp) = aerosol(isect,isp)
         enddo
      enddo

      do i=1,ngas_aq
         gasav(i) = gas(i)
      enddo

C     6) CALCULATION OF TOTAL SULFUR/NITROGEN  MASS BEFORE THE CALL
C     --------------------------------------------------------------

      sulfbef = gas(ngso2)/coefloc
      nitbef  = (gas(nga)+gas(ngn)+gas(ngpan)+gas(ngno3)
     &     +gas(ngno)+gas(ngno2)+gas(nghno2))/coefloc

      do i=1, NS
         nitbef  = nitbef  + aerosol(i,naa)   /mmNH4
     &        + aerosol(i,nan)   /mmNO3
         sulfbef = sulfbef + aerosol(i,na4)   /mmSO4
     &        + aerosol(i,nahso5)/mmHSO5
     &        + aerosol(i,nahmsa)/mmHMSA
      enddo
      nitbef = nitbef  * mmN
      sulfbef= sulfbef * mmS

C     7) Calculation of aerosol number and diameter of each section before aqueous chemistry
C     --------------------------------

      do i = 1,NS
         totmass(i)= aerosol(i,naa)+aerosol(i,na4)+
     &        aerosol(i,nan)+aerosol(i,nac)+aerosol(i,nas)+
     &        aerosol(i,nao)+aerosol(i,nae)+aerosol(i,nar)

         if(INUM.EQ.0) then     ! Case when the number concentration is not followed
            numberconc(i)=totmass(i)/(dsf_aero(i)**3.d0)
     &           /cst_pi6/fixed_rho_aero
         endif
      enddo

C     7bis) Calculation of the bulk distribution factors fdist
      totmasstot = 0.D0
      do i=ifirstact,NS
         totmasstot = totmass(i) + totmasstot
      enddo
      do i=1,ifirstact-1
         fdist(i) = 0.d0
      enddo
      do i=ifirstact,NS
         fdist(i) = totmass(i)/totmasstot
      enddo

C     8) Compute relative humidity
C     ----------------------------

      call COMPUTE_RELATIVE_HUMIDITY(HUMID,temp,press2,rh)
      rh = DMIN1(DMAX1(rh, Threshold_RH_inf), Threshold_RH_sup)

C     ***************************************************
C     Begin computations
C     ***************************************************

      ind_ok = 0

      DO j=1,NITVSRM
         IF (ind_ok.eq.0) THEN

c     10) Simulate cloud/fog
c     ----------------------

            do istep = 1, nistep
               tinit= (istep-1)*deltat
               
               call aqoperator(NS, tinit, deltat, gas, aerosol,
     &              lwc_c, temp, press,
     &              fHCHO,fhcooh,fso2,fh2o2,fNH3,fhno3,fhcl,
     &              ifirstact,
     &              pH,
     &              aqSO2,aqH2O2,
     &              fdist)

            enddo

c     11) Compute TOTAL SULFUR/NITROGEN MASS AFTER THE CALL
C     --------------------------------------------------------------

            sulfaf = gas(ngso2)/coefloc
            nitaf  = (gas(nga)+gas(ngn)+gas(ngpan)+gas(ngno3)
     &           +gas(ngno)+gas(ngno2)+gas(nghno2))/coefloc

            do i=1, NS
               nitaf  = nitaf  + aerosol(i,naa)   /mmNH4
     &              + aerosol(i,nan)   /mmNO3
               sulfaf = sulfaf + aerosol(i,na4)   /mmSO4
     &              + aerosol(i,nahso5)/mmHSO5
     &              + aerosol(i,nahmsa)/mmHMSA
            enddo
            nitaf = nitaf  * mmN
            sulfaf= sulfaf * mmS

            sbal    = DABS(sulfaf/sulfbef-1.D0)

C     12) Check if mass balance is OK
C     Restart otherwise (from step 9)
C     with a timestep divided by 2.
C     -------------------------------

            if ((sulfbef.gt.0.1d0) .and.
     &           (sbal.gt.RTOLSULF) .and.
     &           (initso2.gt.minso2)) then
C               write(*,*) 'VSRM: Sulfate balance not OK: restart'

C     computation has failed
C     ----------------------
               do isect=1,NS
                  do isp=1,naers
                     aerosol(isect,isp) = aerosav(isect,isp)
                  enddo
               enddo

               do i=1,ngas_aq
                  gas(i) = gasav(i)
               enddo

               deltat = deltat/2.d0
               nistep = 2*nistep

            else

C     computation is OK
C     -----------------
               ind_ok = 1
            endif

         ENDIF
      ENDDO

C     13) Compute new distribution by projection to
C     the aerosol population
C     ----------------------------------------------

      do i = 1,NS
         totmass(i)=aerosol(i,naa)+aerosol(i,na4)+
     &        aerosol(i,nan)+aerosol(i,nac)+aerosol(i,nas)+
     &        aerosol(i,nao)+aerosol(i,nae)+aerosol(i,nar)
         ! Calcul of density
         volume(i)=aerosol(i,naa)/liquid_density(naa)
     &        +aerosol(i,na4)/liquid_density(na4)
     &        +aerosol(i,nan)/liquid_density(nan)
     &        +aerosol(i,nac)/liquid_density(nac)
     &        +aerosol(i,nas)/liquid_density(nas)
     &        +aerosol(i,nao)/liquid_density(nao)
     &        +aerosol(i,nae)/liquid_density(nae)
     &        +aerosol(i,nar)/liquid_density(nar)
         if (volume(i) .GT. TINYAQ) then
            rho_aero(i) = totmass(i)/volume(i)
         else
            rho_aero(i) = fixed_rho_aero
         endif

         ! Calcul of diameter for redistaq
         if(numberconc(i).gt.0.d0) then
         dnew(i)=(totmass(i)/numberconc(i)
     &        /rho_aero(i)/cst_pi6)**cst_FRAC3
	 else
	   dnew(i)=dsf_aero(i)
	 endif
      enddo

      NF=NS/NB!Number of composition of each size bin
       !!Mass redistribution
      DO f=1,NF
	DO k=1,NB
	  j=(k-1)*NF+f
	  d(k)=dnew(j)
	  Nub(k)=numberconc(j)
	  DO jesp = 1, naers
	    mass_rdb(k,jesp)=aerosol(j,jesp)
	  ENDDO
	ENDDO

	IF (IREDIST .EQ. 1 .OR. IREDIST .EQ. 2) then
	  call redistaq(NB,DBF_AERO,xbf_aero,fixed_rho_aero,d,mass_rdb)

	  do k = 1,Nb
	      i=(k-1)*NF+f
	      if(INUM.EQ.1) then  ! Case when the number concentration is followed
		totmass(i)=aerosol(i,naa)+aerosol(i,na4)+
     &              aerosol(i,nan)+aerosol(i,nac)+aerosol(i,nas)+
     &              aerosol(i,nao)+aerosol(i,nae)+aerosol(i,nar)

		Nub(k)=totmass(i)/(dsf_aero(i)**3.d0)
     &              /cst_pi6/fixed_rho_aero

	      endif
	  enddo

	ELSEIF (IREDIST .EQ. 3 .OR. IREDIST .EQ. 4 .OR.IREDIST .EQ. 5
     &        .OR.IREDIST .EQ. 6 .OR.IREDIST .EQ. 7 .OR.IREDIST .EQ. 8
     &        .OR.IREDIST .EQ. 9 .OR.IREDIST .EQ. 10) then
	  !! fixed_diameter remplaces dold
	  !! In case dold is out of the box as input ,
	  !! redistribution does not work.
	  call redistribution(NB,naers,naw,DBF_AERO,
     &        fixed_diameter,
     &        IREDIST,section_pass,liquid_density,dqlimit,mass_rdb,
     &        Nub,totQ,with_fixed_density, fixed_rho_aero, dnew)
	ELSE
	  write(*,*) "Choose a redistribution method."
	  stop
	ENDIF

	DO k=1,NB
	  j=(k-1)*NF+f
	  numberconc(j)=Nub(k)
	  totmass(j)=totQ(k)
	  DO jesp = 1, naers
	    aerosol(j,jesp)=mass_rdb(k,jesp)
	  ENDDO
	ENDDO

      ENDDO

      do i = 1,NS
	if(INUM.EQ.1) then ! Case when the number concentration is followed
	  if(numberconc(i).gt.0.d0) then
	    dnew(i) = (totmass(i)/(cst_pi6*RHO_AERO(i)*numberconc(i)))
     &        **cst_FRAC3
	  else
	    dnew(i) = dsf_aero(i)
	  endif
	else
	    numberconc(i)=totmass(i)/(dsf_aero(i)**3.d0)
     &           /cst_pi6/fixed_rho_aero
	    dnew(i)=sqrt(DBF_AERO(i)*DBF_AERO(i+1))

	endif
      enddo


C     14) Update gas & aerosol concentrations
C     ---------------------------------------

      C_GAS(ictmNH3) = gas(nga)   /coefloc2 * mmNH3          !NH3
      C_GAS(ictmHNO3)= gas(ngn)   /coefloc2 * mmHNO3         !HNO3
      C_GAS(ictmHCl) = gas(ngc)   /coefloc2 * mmHCl          !HCl
      C_GAS(ictmSO2) = gas(ngso2) /coefloc2 * mmSO2          !SO2
      C_GAS(ictmH2O2)= gas(ngh2o2)/coefloc2 * mmH2O2         !H2O2
      C_GAS(ictmHCHO)= gas(nghcho)/coefloc2 * mmHCHO         !HCHO
      C_GAS(ictmHNO2)= gas(nghno2)/coefloc2 * mmHNO2         !HNO2
      C_GAS(ictmO3)  = gas(ngo3)  /coefloc2 * mmO3           !O3
      C_GAS(ictmOH)  = gas(ngoh)  /coefloc2 * mmOH           !OH
      C_GAS(ictmHO2) = gas(ngho2) /coefloc2 * mmHO2          !HO2
      C_GAS(ictmNO3) = gas(ngno3) /coefloc2 * mmNO3          !NO3
      C_GAS(ictmNO)  = gas(ngno)  /coefloc2 * mmNO           !NO
      C_GAS(ictmNO2) = gas(ngno2) /coefloc2 * mmNO2          !NO2
      C_GAS(ictmPAN) = gas(ngpan) /coefloc2 * mmPAN          !PAN

      DO i = 1,NS
         C_AER(i,ENa) = aerosol(i,nas) ! Na
         C_AER(i,ESO4)= aerosol(i,na4) ! SO4
         C_AER(i,ENH3)= aerosol(i,naa) ! NH4
         C_AER(i,ENO3)= aerosol(i,nan) ! NO3
         C_AER(i,ECl) = aerosol(i,nac) ! Cl
         C_AER(i,EMD) = aerosol(i,nar) ! DUST
         C_AER(i,EBC) = aerosol(i,nae) ! EC
         DO j = 1,nesp_aec
            jesp = aec_species(j)
            C_AER(i,jesp) = aerosol(i,nao)*foa(j,i)
         ENDDO
         C_AER(i,EH2O)= aerosol(i,naw) ! H2O
      ENDDO

      do i = 1, NS
         C_NUM(i) = numberconc(i)
      enddo

C     15) Compute in-cloud scavenging
C     -------------------------------

      IF (rain_rate.GT.0.d0) THEN
                                ! compute scavenging coefficient, same for all species
         drop_diam = 9.76d-4 * (rain_rate**0.21d0 ) ! m
         scav_coef = 4.17d-7 * collision_eff * rain_rate / drop_diam ! s-1

                                ! compute scavenged fraction of droplet
         scav_frac = 1.D0 - dexp(scav_coef * (t0 - t1) ) ! adim

                                ! compute scavenged quantiy
         DO i = ifirstact,NS
            DO isp = E1,E2
               qscav_aer(i,isp) = C_AER(i,isp) * scav_frac
               C_AER(i,isp) = C_AER(i,isp)- qscav_aer(i,isp)
            ENDDO
            qscav_num(i) = C_NUM(i) * scav_frac
            C_NUM(i) = C_NUM(i) * (1.d0 - scav_frac)
         ENDDO

         qscav_gas(ictmSO2) = aqSO2 * scav_frac * 0.7901d0
         qscav_gas(ictmH2O2)= aqH2O2 * scav_frac

                                ! remove scavenged quantity
         C_GAS(ictmSO2) = C_GAS(ictmSO2) - qscav_gas(ictmSO2)
         C_GAS(ictmH2O2)= C_GAS(ictmH2O2) - qscav_gas(ictmH2O2)

      ENDIF
! 
! C     16) Treshold TINYAQ.
! C----
!       DO i=1,NGAS!problem of dimension
!          C_GAS(i)=dmax1(C_GAS(i),TINYAQ)
!       ENDDO
!       DO i=1,NS!problem of dimension
!         DO jesp= 1, NAER
!          C_AER(i,jesp)=dmax1(C_AER(i,jesp),TINYAQ)
!         ENDDO
!       ENDDO
!       DO i=1,NS
!          C_NUM(i) = dmax1(C_NUM(i), TINYNAQ)
!       ENDDO

      
      return
      end
