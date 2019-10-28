C-----------------------------------------------------------------------
C     Copyright (C) 2001-2012, ENPC - INRIA - EDF R&D
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



      SUBROUTINE chem (ns,nr,nrphot,nreactphot,nemis,nemisspecies,
     $     convers_factor,convers_factor_jac,ts,DLattenuation,DLhumid,
     $     DLtemp,DLpress,DLCsourc,DLCphotolysis_rates,delta_t,
     $     DLattenuationf,DLhumidf,DLtempf, DLpressf,DLCsourcf,
     $     DLCphotolysis_ratesf,ncycle,dlon,dlat,DLconc,
     $     option_adaptive_time_step, ATOL, tstep_min,
     $     option_photolysis, option_chemistry, delta_tmax)

C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     This routine computes one timestep for gas-phase chemistry RACM.
C     Chemical kinetics is solved in a grid cell.
C
C------------------------------------------------------------------------
C
C     -- INPUT VARIABLES
C
C     TS: initial time (GMT, computed from January 1st, [s]).
C     DLATTENUATION: cloud attenuation at initial time.
C     DLHUMID: specific humidity at initial time ([%]).
C     DLTEMP: temperature at initial time ([K]).
C     DLPRESS: pressure at initial time ([Pa]).
C     DLCSOURC: array of chemical volumic emissions at initial time
C     # ([\mu.g/m^3/s]).
C     DLCPHOTOLYSIS_RATES: photochemical kinetic rates
C     # at initial time ([s^{-1}]).
C     DELTA_T: time step ([s]).
C     option_adaptive_time_step: 1 if adaptive time step.
C     ATOL:  relative tolerance for deciding if the time step is kept.
C     tstep_min: minimum time step.
C     delta_tmax: maximum time step.
C     The same variables are defined at final time of the timestep.
C     'f' is then put at the end of the name.
C
C     -- INPUT/OUTPUT VARIABLES
C
C     DLCONC: array of chemical concentrations ([\mu.g/m^3]).
C     # Before entry, it is given at initial time of the timestep.
C     # On exit, it is computed at final time of the timestep.
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
C     2009/01/22: added adaptatif time stepping (K. Sartelet, CEREA)
C
C     2008/06/17: removed computation of conversion factors to speed up
C     computation (Meryem Ahmed de Biasi, INRIA).
C
C     2008/04/02: suppressed the loop on coordinates (I. Korsakissok, CEREA).
C     2002/02/26: new treatment of sources (Jaouad Boutahar, CEREA).
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Denis Qu�lo, CEREA, June 2001.
C
C------------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION ts,delta_t
      DOUBLE PRECISION tschem,tfchem

      integer ns,nr,nrphot,nemis

      DOUBLE PRECISION DLconc(ns),ZC(ns)
      DOUBLE PRECISION DLtemp,DLtempf
      DOUBLE PRECISION DLattenuation
      DOUBLE PRECISION DLattenuationf
      DOUBLE PRECISION DLhumid,DLhumidf
      DOUBLE PRECISION DLCsourc(Nemis)
      DOUBLE PRECISION DLCsourcf(Nemis)
      DOUBLE PRECISION ZCsourc(ns)
      DOUBLE PRECISION ZCsourcf(ns)
      DOUBLE PRECISION DLRki(Nr),DLRkf(Nr)
      DOUBLE PRECISION DLpress,DLpressf
      DOUBLE PRECISION DLCphotolysis_rates(NRphot)
      DOUBLE PRECISION DLCphotolysis_ratesf(NRphot)

      double precision dlon,dlat

      integer ncycle, ncycle_chem
      double precision convers_factor(ns)
      double precision convers_factor_jac(ns,ns)

      DOUBLE PRECISION Zangzen,Zangzenf
      DOUBLE PRECISION Zatt,Zattf

      DOUBLE PRECISION muzero,DLmuzero
      EXTERNAL muzero

      double precision pi

      integer nreactphot(nrphot)
      integer nemisspecies(nemis)

      INTEGER Jt,Jsp,i

      DOUBLE PRECISION DLk1(ns), DLk2(ns),ZC_old(ns)
      DOUBLE PRECISION EPSDLK
      PARAMETER (EPSDLK = 1.D-15)
      DOUBLE PRECISION supEdtstep, Edtstep(ns),ATOL
      DOUBLE PRECISION tstep,tstep_new,tstep_min,tfchem_tmp, delta_tmax

      INTEGER option_adaptive_time_step
      INTEGER option_photolysis, option_chemistry

C     Constants.
      pi = 3.14159265358979323846D0

C     Spatial extraction for volumic sources.

      DO Jsp=1,ns
         ZCsourc(jsp)=0.D0
         ZCsourcf(jsp)=0.d0
      ENDDO

      DO Jsp=1,Nemis
         ZCsourc(nemisspecies(jsp)+1)=DLCsourc(Jsp)
         ZCsourcf(nemisspecies(jsp)+1)=
     $        DLCsourcf(Jsp)
      ENDDO

C     Cloud attenuation.

      Zatt = DLattenuation
      Zattf = DLattenuationf

C     Projection.

      DO Jsp=1,ns
         ZC(Jsp) = DLconc(Jsp)
      ENDDO

C     Integration of chemistry (eventually with subcycling).
      Ncycle_chem=Ncycle
      IF(option_adaptive_time_step.EQ.1) THEN
         IF(delta_t>delta_tmax) THEN
            Ncycle_chem= ceiling(Ncycle * delta_t / delta_tmax)
         ENDIF
      ENDIF

      DO Jt=1,Ncycle_chem
         tschem=ts+(Jt-1)*delta_t/Ncycle_chem
         tfchem=tschem+delta_t/Ncycle_chem
         tstep = delta_t/Ncycle_chem
         tfchem_tmp = tfchem


         Do while (tschem.LT.tfchem)
C     If option_photolysis is 1,
C     photolytic reactions are calculated in kinetic.f
            DLmuzero=muzero(tschem,Dlon,Dlat)
            Zangzen=dabs(DACOS(DLmuzero)*180.D0/PI)
            DLmuzero=muzero(tfchem_tmp,Dlon,Dlat)
            Zangzenf=dabs(DACOS(DLmuzero)*180.D0/PI)


            IF (option_chemistry.eq.1) then
              CALL Kinetic_racm(nr,DLRKi,DLtemp,DLhumid,DLpress,Zangzen,
     s           Zatt,option_photolysis)
              CALL Kinetic_racm(nr, DLRKf,DLtempf,DLhumidf,DLpressf,
     s           Zangzenf,Zattf,option_photolysis)
            ELSE IF (option_chemistry.eq.2) then
             CALL Kinetic_racm2(nr,DLRKi,DLtemp,DLhumid,DLpress,Zangzen,
     s           Zatt,option_photolysis)
             CALL Kinetic_racm2(nr, DLRKf,DLtempf,DLhumidf,DLpressf,
     s           Zangzenf,Zattf,option_photolysis)
            ELSE IF (option_chemistry.eq.3) then
             CALL Kinetic_cb05(nr,DLRKi,DLtemp,DLhumid,DLpress,Zangzen,
     s           Zatt,option_photolysis)
             CALL Kinetic_cb05(nr, DLRKf,DLtempf,DLhumidf,DLpressf,
     s           Zangzenf,Zattf,option_photolysis)
            ELSE IF (option_chemistry.eq.4) then
             CALL Kinetic_leighton(nr,DLRKi,DLtemp,DLhumid,DLpress,
     s              Zangzen,
     s           Zatt,option_photolysis)
             CALL Kinetic_leighton(nr, DLRKf,DLtempf,DLhumidf,DLpressf,
     s           Zangzenf,Zattf,option_photolysis)
            ENDIF


C     If option_photolysis is 2,
C     photolytic reactions may be read.
            IF (option_photolysis.eq.2) then
c               write(*,*) "==== read photolysis constants ===="
            DO i=1,Nrphot
               DLRKi(Nreactphot(i)+1) = Zatt *
     $              DLCphotolysis_rates(i)
               DLRKf(Nreactphot(i)+1) = Zattf *
     $              DLCphotolysis_ratesf(i)
            ENDDO


            ENDIF
                                ! Solve gas-phase chemistry for the time step
            CALL roschem(ns, nr, ZC,ZCsourc,ZCsourcf,
     $           convers_factor, convers_factor_jac,tschem
     $           ,tfchem_tmp,DLRki,DLRkf,ZC_old,DLK1,DLK2,
     $           option_chemistry)

C     Integration of chemistry with adaptive time stepping
            IF(option_adaptive_time_step.EQ.1) then
C     Check that the time step was ok
               supEdtstep = 0.D0
               Do Jsp = 1,Ns
                  If((DLK1(Jsp).GT.EPSDLK
     &                 .OR.DLK2(Jsp).GT.EPSDLK)
     &                 .AND.ZC(Jsp).GT.EPSDLK) then
                                ! Estimate the relative error
                     Edtstep(Jsp) = 0.5D0 *
     &                    dabs(DLk1(Jsp) + DLk2(Jsp))
     &                    / ZC(Jsp)
                     If(Edtstep(Jsp).GT.supEdtstep) then
                        supEdtstep = Edtstep(Jsp)
                     Endif
                  Endif
               Enddo
               supEdtstep = supEdtstep/ATOL * tstep

               If(supEdtstep.GT.1.D0
     &                 .AND.tstep.GT.tstep_min) then
! The time step is rejected and the computation
! is redone with a smaller time step
                  tstep_new = tstep * 0.9d0 /dsqrt(supEdtstep)
                  tstep_new = DMAX1(tstep_new,tstep_min)
                  tfchem_tmp = tschem + tstep_new
                  If(tfchem_tmp.GT.tfchem) then
                     tfchem_tmp = tfchem
                     tstep_new = tfchem_tmp - tschem
                  Endif
                  tstep = tstep_new
                  Do Jsp=1,ns
                     ZC(Jsp) = ZC_old(Jsp)
                  Enddo
               Else
                                ! The time step is accepted and the time is incremented
                  tschem = tfchem_tmp
                  tfchem_tmp = tfchem_tmp + tstep
                  If(tfchem_tmp.GT.tfchem) then
                     tfchem_tmp = tfchem
                     tstep = tfchem_tmp - tschem
                  Endif
               Endif
            ELSE
               tschem = tfchem
            ENDIF

         Enddo                  !End loop Do while for time stepping

      ENDDO

C     Storage in the 3D array of chemical concentrations.

      DO i=1,ns
         DLconc(i) = ZC(i)
      ENDDO

      END
