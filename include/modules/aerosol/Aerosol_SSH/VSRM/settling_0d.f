C-----------------------------------------------------------------------
C     Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
C     Author(s): Kathleen Fahey and Bruno Sportisse
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

      subroutine settling_0d(Nesp_aer,nbin_aer,
     &     dsf_aero,DLconc_aer,lwcavg,t0,t1,heightfog)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This subroutine computes settling for fog droplets on the basis
C     of a column of fog.
C     A fog is defined if LWC in the first cell is greater than a 
C     threshold LWCFOG. The column of fog is defined above the first cell
C     with this criterion.
C     The parameterization for the scavenging coefficient is taken from 
C     the PhD work of Kathleen Fahey as:
C     
C     lambda =0.014*(lwcavg**1.67d0)+0.009d0*(lwcavg**1.08d0)
C     &              /(lwcavg*heightfog)
C     
C     with LWCAVG the average fog liquid water content and HEIGHTFOG the
C     height of the fog column.
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     LWC: liquid water content ([]).
C     T0/T1: initial/final time ([s]).
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     DLCONC: concentrations.
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
C     1) 2005/11/05: rewrite the subroutine (B.Sportisse).
C------------------------------------------------------------------------
C     
C     -- AUTHOR(S)
C     
C     2004: Kathleen Fahey, CEREA.
C     2005/10/3, cleaning, Bruno Sportisse, CEREA.
C     
C------------------------------------------------------------------------

      IMPLICIT NONE
      
      include 'aerpar.inc'
      include 'droppar.inc'

      INTEGER Nx,Ny,Nz,Nesp_aer,Nbin_aer
      double precision DLconc_aer(Nbin_aer,Nesp_aer)
      double precision lwc, lwc_surf, lwc_min
      double precision scavfog
      double precision t0,t1,dtloc
      double precision heightfog,lwcavg
      double precision lambdafog
      integer Jk,Jj,Ji,jspec,jbin,indok
      logical fog

      integer ifirstact

c      DOUBLE PRECISION maillage(Nz+1)
      double precision dsf_aero(Nbin_aer)

      INTEGER i1,i2,j1,j2

      dtloc = t1-t0

c            heightfog = 0.D0
c           nfoglay  = 0
c          lwcavg    = 0.D0

c            IF (lwc.GT.lwcfog) THEN
C     1) Define the fog layers
C     once LWC(1)>LWCFOG
C     ------------------------      

c               indok = 1
c               DO Jk = 1,Nz
c                  IF ((lwc(Ji,Jj,Jk).GT.lwc_min).AND.
c     &                 (indok.eq.1)) THEN 
c                     heightfog = maillage(Jk+1)
c                     nfoglay  = nfoglay+1
c                     lwcavg    = lwcavg+lwc(Ji,Jj,Jk)
c                  ELSE
c                     indok=0
c                  ENDIF
c               ENDDO
c               IF (nfoglay.GE.1) THEN
c                  lwcavg = lwcavg/nfoglay
                  
C     2) Compute the removal rate 
C     parameterized as a function of 
C     the averaged fog LWC (LWCAVG)
C     ------------------------------
c              if (fog) then 
c                 write(*,*) "fog: ", fog
c                  write(*,*) "lwcavg:", lwcavg,
c     &             "heightfog:", heightfog                    
                  lambdafog =0.014d0 * (lwcavg**1.67d0)
     &                 + 0.009d0 * (lwcavg**1.08d0)
     &                 / (lwcavg * heightfog)
                  scavfog = DEXP(-lambdafog*dtloc)
c                  stop
c                  write(*,*) "lambdafog:", lambdafog,
c     &             "scavfog:", scavfog  
C     3) Apply the removal process for the activated bins
C     ---------------------------------------------------
                  
c                  DO Jk = 1,nfoglay
C     Subroutine initactiv called here in the eventuality of a future
C     physical parametrisation of ifirstact 
                     call initactiv(Nbin_aer,dsf_aero,ifirstact)
                     do jspec = 1,Nesp_aer
                        do jbin = ifirstact,Nbin_aer
                           DLCONC_aer(jbin,jspec)= 
     &                          DLCONC_aer(jbin,jspec)*scavfog
                        ENDDO   
                     ENDDO
c               endif
c                  ENDDO
c               ENDIF
c            ENDIF
C     end fog case
C     end loop along columns
      
      return
      end    
      
