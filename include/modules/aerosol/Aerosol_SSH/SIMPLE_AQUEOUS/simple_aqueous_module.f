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

      subroutine simple_aqueous_module(NGAS,NAER,NS,NB,NICD,
     &     RHO_AERO,fixed_rho_aero,DBF_AERO,dsf_aero,xbf_aero,C_GAS,
     &     C_AER,HUMID,press_in_pa,temp,lwc_c,t0,t1,rain_rate,pH,
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
C     PRESS_IN_ATM: pressure             ([Pa]).
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
C    Shupeng ZHU: add composition support for SCRAM (Oct.2014)
C
C------------------------------------------------------------------------
C     
C     -- AUTHOR(S)
C     
C     Marilyne Tombette, CEREA, , on the basis of Pandis/Seinfeld book, 
C     pp. 
C     
C------------------------------------------------------------------------

      IMPLICIT NONE

      include 'param.inc'
      include 'CONST_A.INC'
      include 'data_mass.inc'
      include 'aerpar_simple.inc'
      include 'droppar_simple.inc'
      include 'varp.inc'
      include 'varp_cloud.inc'
      include 'num_aq.inc'

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
      integer with_fixed_density ! YK
      
      double precision d(NB)
      double precision mass_rdb(NB,naers)
      double precision totQ(NB)
      double precision Nub(NB)
      double precision dqlimit
      double precision dold(NS)
      double precision volume(NS)      

      double precision coefloc,tinit,coefloc2, coef_cte
      double precision dh2SO4,dso2
      double precision gas(ngas_aq), aerosol(NS,naers)!naers = 9
      double precision gasav(ngas_aq), aerosav(NS,naers)
      double precision rh,temp,press_in_atm,lwc_c,lwc,t0,t1,initso2
      double precision dnew(NS),deltat,deltatmax,press_in_pa
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
      double precision scav_frac,aqSO2,aqH2O2,aqO3
      double precision DBF_AERO(NB+1),RHO_AERO(NS)
      double precision dsf_aero(NS),xbf_aero(NB+1)
      double precision fixed_rho_aero,totmasstot
      double precision fixed_diameter(NS)
      
      double precision dsivdt_o3,dsivdt_h2o2,dsivdt
      double precision fc_o3, fc_h2o2, fc_siv
      double precision tot_sulf
      double precision tot_nit,tot_amm,fo3,fh2o2
      double precision ptotso2
      double precision conc_in(5), conc_out(10)      

      double precision henrys(ngas_aq)
      double precision cst_dissoc(nreact_dissoc)
      double precision cst_oxydation(nreact_oxydation)
      double precision tnitold, ammonold,sulfateold
      double precision so2g,henrys_eff_so2

      double precision cst_pr,cst_pr2,cst_RT, cst_RTLWC     

      data collision_eff /0.9d0/
      
c     Initialisation
c     --------------
      cst_pr = 1.987d-3         ! perfect gas constant in kcal.K-1.mol-1
      cst_pr2 = 8.206d-2        ! perfect gas constant in atm.L.K-1.mol-1
      cst_RT = cst_pr2 * temp !*temperature
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
    
      !nesp = NGAS+NS*NAER

      press_in_atm  = press_in_pa/101325.d0 ! pressure in atm
      lwc = lwc_c * 1.d-6       ! in L_water/L_air
      cst_RTLWC = cst_pr2 * temp * lwc ! in atm.mol-1.L_water

      call initactiv(NS,dsf_aero,ifirstact)!within VSRM
      !computes the label for the first activated bin [ifirstact,NS]

      nistep = NITSUBAQ
      deltat = (t1 - t0) / nistep !timestep in s

      tot_sulf=0.d0

      do i=1,ngas_aq
         gas(i)=0.d0
      enddo
      
      do isect=1,NS
         totmass(isect)=0.d0
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

      do j=1,NS
         do isp=1,naers
            aerosol(j,isp)=0.d0
         enddo
      enddo

C     Calculation of aerosol number and diameter of each section before aqueous chemistry
C     --------------------------------
      DO i = 1, NB
         fixed_diameter(i) = sqrt(DBF_AERO(i)*DBF_AERO(i+1))
      ENDDO


C     Compute henry's and equilibrium constants.
C     ------------------------------------------
      
      coef_cte= (1.d0/temp - 1.d0/temp_ref)

      do i=1,ngas_aq
         henrys(i) = chenry(i) * exp(-dhhenry(i)*coef_cte / cst_pr)
      enddo

      do i=1,nreact_dissoc
         cst_dissoc(i) = ckdissoc(i) 
     $        * exp(dhkdissoc(i)*coef_cte)
      enddo

      do i=1,nreact_oxydation
         cst_oxydation(i) = ckoxydation(i) 
     $        * exp(-dhkoxydation(i)*coef_cte 
     $        / cst_pr)
      enddo
      
C     Compute sectional diameters and bimodal distribution
C     -------------------------------------------------------

C     convert \mu g.m-3 -> atm
      gas(igso2)  = C_GAS(ictmSO2)  *1.d-9*cst_RT/mmSO2  !SO2
      gas(ignh3)  = C_GAS(ictmNH3)  *1.d-9*cst_RT/mmNH3  !NH3
      gas(ighno3) = C_GAS(ictmHNO3) *1.d-9*cst_RT/mmHNO3 !HNO3
      gas(igh2o2) = C_GAS(ictmH2O2) *1.d-9*cst_RT/mmH2O2 !H2O2
      gas(igo3)   = C_GAS(ictmO3)   *1.d-9*cst_RT/mmO3   !O3
      gas(igCO2) = 350.d0 * 1.d-6 * press_in_atm ! 350 ppm in averaged.

C     totmass en \mu g.m-3 without sulfate
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

         DO j = 1,n_aec
            jesp = aec_species(j)
            aerosol(i,nao) = aerosol(i,nao) + C_AER(i,jesp)
            foa(j,i) = C_AER(i,jesp)
            ftot = ftot + C_AER(i,jesp)
         ENDDO

         DO j=1, n_aec
            if (ftot .GT. 0.d0) then
               foa(j,i)=foa(j,i)/ftot
            else
               foa(j,i)=0.d0
            endif
         ENDDO

         aerosol(i,naw) = C_AER(i,EH2O) ! H2O
      ENDDO
      
C     Transfer all H2SO4(g) to aerosol phase      
C     ------------------------------------------
      dh2SO4 = (C_GAS(ictmH2SO4)/NS) * mmSO4/mmH2SO4
      DO i=1,NS
         aerosol(i,na4) = aerosol(i,na4)+dh2SO4
      ENDDO
      C_GAS(ictmH2SO4) = 0.d0

 !!!core problem of DBF
      DO i = 1,NS
	 b=ID_SIZE(i,1)
         totmass(i)=aerosol(i,naa)+aerosol(i,na4)+
     &        aerosol(i,nan)+aerosol(i,nac)+aerosol(i,nas)+
     &        aerosol(i,nao)+aerosol(i,nae)+aerosol(i,nar)
         if(INUM.EQ.1) then
            numberconc(i) = C_NUM(i)
            ! Calcul of diameter for redistribution euler
	     if (numberconc(i) .GT. 0.d0) then
		dold(i) = (totmass(i)/(cst_pi6*RHO_AERO(i)*numberconc(i)))
     &           **cst_FRAC3
             else
               dold(i)=sqrt(DBF_AERO(b)*DBF_AERO(b+1))
            endif
         else

            if (with_fixed_density .ne. 1) then
               numberconc(i)=totmass(i)/(dsf_aero(i)**3.d0)
     &              /cst_pi6/RHO_AERO(i)
            else
               numberconc(i)=totmass(i)/(dsf_aero(i)**3.d0)
     &              /cst_pi6/fixed_rho_aero
            endif
            dold(i)=sqrt(DBF_AERO(b)*DBF_AERO(b+1))
         endif

         !!! added by SZ
         if (dold(i) .LT. DBF_AERO(b) 
     $        .or. dold(i) .GT. DBF_AERO(b+1)) THEN
            dold(i)=sqrt(DBF_AERO(b)*DBF_AERO(b+1))
         endif

      ENDDO

C     Calculation of the bulk distribution factors fdist
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
     
C     Compute total sulfate, ammonium and nitrate.
C     -------------------------------------------
      tot_sulf = 0.d0
      tot_nit = 0.d0
      tot_amm = 0.d0

      DO i=ifirstact,NS
         tot_sulf = tot_sulf + aerosol(i,na4) 
         tot_nit = tot_nit + aerosol(i,nan) 
         tot_amm = tot_amm + aerosol(i,naa) 
      enddo
      
      sulfateold=tot_sulf
      tnitold=tot_nit
      ammonold=tot_amm

C     ***************************************************
C     Begin computations
C     ***************************************************

c     Simulate cloud/fog
c     ----------------------

      conc_in(1) = gas(igso2)              ! S(IV) in atm
      conc_in(2) = tot_sulf * 1.d-9 * cst_RT /mmSO4 ! S(VI) in atm
      conc_in(3) = tot_amm * 1.d-9 * cst_RT / mmNH4 + gas(ignh3) ! in atm
      conc_in(4) = tot_nit * 1.d-9 * cst_RT / mmNO3 + gas(ighno3) ! in atm
      conc_in(5) = gas(igco2) ! in atm

C     First computation of pH
C     -----------------------
      call compute_ph(conc_in,henrys,cst_dissoc,temp,
     $     lwc,ph,conc_out)

C     First computation of O3 and H2O2 aqueous concentration.
C     An instantaneous equilibrium between the gas phase concentration
C     and the aqueous phase concentration is assumed.
            aqo3 = gas(igo3) * henrys(igo3) /
     $           (1.d0 + cst_RTLWC * henrys(igo3)) ! in mol.L-1_water

            aqh2o2 = gas(igh2o2) * henrys(igh2o2) /
     $           (1.d0 + cst_RTLWC * henrys(igh2o2)) ! in mol.L-1_water


C     Loop on time step for aqueous module.
C     -------------------------------------

      do istep = 1,nistep

C     No aqueous chemistry without S(IV).
         if (conc_in(1).gt.0.d0) then

C     Compute reaction rates of S(IV) with O3 and H2O2
C     ------------------------------------------------

C     rate Seinfel Pandis for O3

         dsivdt_o3 = - (cst_oxydation(1)*conc_out(2) !SO2 aq
     $        + cst_oxydation(2)*conc_out(3) !HSO3-
     $        + cst_oxydation(3)*conc_out(4)) !SO3--
     $        * aqo3 

C     rate VSRM for H2O2

         dsivdt_h2o2 = - 1300000.d0* exp( -4430.d0* coef_cte)
     $        * conc_out(2) * aqh2o2
     $        / (1.d0+16.d0*conc_out(1))

C     total rate

         dsivdt = dsivdt_o3 + dsivdt_h2o2

C     Compute concentrations.
C     -----------------------

C     Forecast of S(IV), O3 and H2O2 concentrations.
            fc_siv = conc_in(1) + dsivdt * deltat * cst_RTLWC ! in atm
            fc_o3 = gas(igo3) + dsivdt_o3 * deltat * cst_RTLWC ! in atm
            fc_h2o2 = gas(igh2o2) + dsivdt_h2o2 * deltat * cst_RTLWC ! in atm

            if ( (fc_siv.lt.0.d0) .or.
     $           ( (fc_o3.lt.0.d0) .or. (fc_h2o2.lt.0.d0) ) ) then

C     Determine the maximum time step for the most limiting reactant among
C     SO2, O3 and H2O2.
               deltatmax = DMIN1( (conc_in(1) / cst_RTLWC / (-dsivdt)),
     $              DMIN1((gas(igo3) / cst_RTLWC / (-dsivdt_o3)),
     $              (gas(igh2o2) / cst_RTLWC / (-dsivdt_h2o2))))

               conc_in(1) = DMAX1( 0.d0,
     $              (conc_in(1) + dsivdt * deltatmax * cst_RTLWC ))
               gas(igo3) = DMAX1( 0.d0,
     $              (gas(igo3) + dsivdt_o3 * deltatmax * cst_RTLWC ))
               gas(igh2o2) = DMAX1( 0.d0,
     $              (gas(igh2o2) + dsivdt_h2o2 * deltatmax * cst_RTLWC))
               conc_in(2) = conc_in(2) - dsivdt * deltatmax * cst_RTLWC
            else
               conc_in(1) = fc_siv
               gas(igo3) = fc_o3
               gas(igh2o2) = fc_h2o2
               conc_in(2) = conc_in(2) - dsivdt * deltat * cst_RTLWC
            endif

c$$$
c$$$C     Update aqueous concentration and pH.
c$$$C     ------------------------------------
c$$$
c$$$C     An instantaneous equilibrium between the gas phase concentration
c$$$C     and the aqueous phase concentration is assumed.
c$$$            aqh2o2 = gas(igh2o2) * henrys(igh2o2) /
c$$$     $           (1.d0 + cst_RTLWC * henrys(igh2o2)) ! in mol.L-1_water
c$$$
c$$$            aqo3 = gas(igo3) * henrys(igo3) /
c$$$     $           (1.d0 + cst_RTLWC * henrys(igo3)) ! in mol.L-1_water
c$$$
c$$$C     Update pH.
c$$$            call compute_ph(conc_in, henrys, cst_dissoc,
c$$$     $           temp, lwc, ph, conc_out)
c$$$
         endif
      enddo

C     New total sulfate, nitrate and ammonium.
C     ----------------------------------------

C     Unlike nitrate and ammonium, sulfate is not stored in conc_out.
      tot_sulf = conc_in(2) * mmSO4 * 1.d9 / cst_RT !  in \mu g.m-3
      tot_nit= conc_out(6) * mmNO3 * 1.d9 *lwc ! in \mu g.m-3
      tot_amm= conc_out(8) * mmNH4 * 1.d9 *lwc ! in \mu g.m-3

C     Dissolved S(IV) is added to gas phase SO2.
      gas(igso2) = conc_in(1) ! in atm
C     Portion of S(IV) concerned by in-cloud scavenging.
      aqso2 = conc_out(2) + conc_out(3) + conc_out(4) ! in mol.L-1

C     Dissolved HNO3 is added to gas phase HNO3.
      gas(ighno3) = conc_out(9) + conc_out(5) * cst_RTLWC ! in atm

C     Dissolved NH3 is added to gas phase NH3.
      gas(ignh3) = conc_out(10) + conc_out(7) * cst_RTLWC ! in atm      

C     Compute new distribution by projection to
C     the aerosol population (as in VSRM).
C     ----------------------------------------------

      do isect=1,NS
         aerosol(isect,nan)    = aerosol(isect,nan)
     &        + fdist(isect) * (tot_nit - tnitold) ! NITRATE (aq) 
         aerosol(isect,naa)    = aerosol(isect,naa) 
     &        + fdist(isect) * (tot_amm - ammonold) ! AMMONIUM (aq)
         aerosol(isect,na4)    = aerosol(isect,na4)
     &        + fdist(isect) * (tot_sulf - sulfateold) ! SULFATE (aq)
      enddo

      do isp=1, naers
         do isect=1,NS
            aerosol(isect,isp)=dmax1(aerosol(isect,isp),tinyaq2)
         enddo
      enddo     

      do i = 1,NS
         totmass(i)=aerosol(i,naa)+aerosol(i,na4)+
     &        aerosol(i,nan)+aerosol(i,nac)+aerosol(i,nas)+
     &        aerosol(i,nao)+aerosol(i,nae)+aerosol(i,nar)
         dnew(i) = (totmass(i) / numberconc(i) / rho_aero(i)
     $        / cst_pi6)**cst_FRAC3
      enddo

      NF=NS/NB!number of fraction
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
     &        .OR.IREDIST .EQ. 9 .OR.IREDIST .EQ. 10.OR.IREDIST .EQ. 11
     &        .OR.IREDIST .EQ. 12.OR.IREDIST.EQ.13) then
	  !! fixed_diameter remplaces dold
	  !! In case dold is out of the box as input ,
	  !! redistribution does not work.
	  call redistribution(NB,naers,naw,DBF_AERO,
     &        fixed_diameter,
     &        IREDIST,section_pass,liquid_density,dqlimit,mass_rdb,
     &        Nub,totQ,with_fixed_density,fixed_rho_aero,dnew)
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

c$$$      do i = 1,NS
c$$$	if(INUM.EQ.1) then ! Case when the number concentration is followed
c$$$	  if(numberconc(i).gt.0.d0) then
c$$$	    dnew(i) = (totmass(i)/(cst_pi6*RHO_AERO(i)*numberconc(i)))
c$$$     &        **cst_FRAC3
c$$$	  else
c$$$	    dnew(i) = dsf_aero(i)
c$$$	  endif
c$$$	else
c$$$	    numberconc(i)=totmass(i)/(dsf_aero(i)**3.d0)
c$$$     &           /cst_pi6/fixed_rho_aero
c$$$	    dnew(i)=sqrt(DBF_AERO(i)*DBF_AERO(i+1))
c$$$
c$$$	endif
c$$$      enddo !! test YK
      
C     14) Update gas & aerosol concentrations in model.
C     -------------------------------------------------

      C_GAS(ictmHNO3)= gas(ighno3) * 1.d9 * mmHNO3 / cst_RT !HNO3
      C_GAS(ictmSO2) = gas(igso2)  * 1.d9 * mmSO2  / cst_RT !SO2
      C_GAS(ictmH2O2)= gas(igh2o2) * 1.d9 * mmH2O2 / cst_RT !H2O2
      C_GAS(ictmO3)  = gas(igo3)   * 1.d9 * mmO3   / cst_RT !O3
      DO i = 1,NS
         C_AER(i,ENa) = aerosol(i,nas) ! Na
         C_AER(i,ESO4)= aerosol(i,na4) ! SO4
         C_AER(i,ENH3)= aerosol(i,naa) ! NH4
         C_AER(i,ENO3)= aerosol(i,nan) ! NO3
         C_AER(i,ECl) = aerosol(i,nac) ! Cl
         C_AER(i,EMD) = aerosol(i,nar) ! DUST
         C_AER(i,EBC) = aerosol(i,nae) ! EC
         DO j = 1,n_aec
            jesp = aec_species(j)
            C_AER(i,jesp) = aerosol(i,nao)*foa(j,i)
         ENDDO
         C_AER(i,EH2O)= aerosol(i,naw) ! H2O
      ENDDO
      
      do i = 1, NS
         C_NUM(i) = numberconc(i)
      enddo

C     Computes in-cloud scavenging.
C     -----------------------------

      IF (rain_rate.GT.0.d0) THEN
! compute scavenging coefficient, same for all species
         drop_diam = 9.76d-4 * (rain_rate**0.21d0 ) ! m
         scav_coef = 4.17d-7 * collision_eff * rain_rate / drop_diam ! s-1

! compute scavenged fraction of droplet
         scav_frac = 1.D0 - dexp(scav_coef * (t0 - t1) ) ! adim

! compute scavenged quantity

         DO i = ifirstact,NS
            DO isp = E1,E2
               qscav_aer(i,isp) = C_AER(i,isp) * scav_frac
            ENDDO
         ENDDO
         qscav_gas(ictmSO2) =  aqSO2 * lwc * 1.d9 * mmSO2
     $        * scav_frac * 0.7901d0 ! Remarks: what is the meaning of the coefficient 0.7901d0?
         qscav_gas(ictmH2O2)= aqH2O2 * lwc * 1.d9 * mmH2O2 * scav_frac

C     Number concentration
         DO i=1,NS
            qscav_num(i) = C_NUM(i) * scav_frac
         ENDDO

! remove scavenged quantity
         DO i = ifirstact,NS
            DO isp = E1,E2
               C_AER(i,isp) = C_AER(i,isp)- qscav_aer(i,isp)
            ENDDO
         ENDDO
         C_GAS(ictmSO2) = C_GAS(ictmSO2) - qscav_gas(ictmSO2)
         C_GAS(ictmH2O2)= C_GAS(ictmH2O2) - qscav_gas(ictmH2O2)

C     Number concentration
          DO i=1,NS
            C_NUM(i) = C_NUM(i) - qscav_num(i)
         ENDDO
      ENDIF

C     Treshold TINYAQ.
C---- 
!       DO i=1,NGAS!problem of dimension
!          C_GAS(i)=dmax1(C_GAS(i),TINYAQ)
!       ENDDO
!       DO i=1,NS!problem of dimension
! 	DO jesp= 1, NAER
!          C_AER(i,jesp)=dmax1(C_AER(i,jesp),TINYAQ)
!         ENDDO
!       ENDDO      
!       DO i=1,NS
!          C_NUM(i) = dmax1(C_NUM(i), TINYNAQ)
!       ENDDO

      return
      end
