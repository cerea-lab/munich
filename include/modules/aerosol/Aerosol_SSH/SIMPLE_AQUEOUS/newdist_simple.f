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

      subroutine newdist_simple(NS,dbf_aero,dsf_aero,fdistx,
     $     fdistx2,ifirstact)
      
C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This routine computes the activated distribution of aerosols 
C     (before the aqueous-phase model). This is projected onto a
C     fixed bimodal lognormal distribution.
C     
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
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
C     Adapted from VSRM to the simple aqueous module.
C     (Marilyne Tombette, 2007).
C
C------------------------------------------------------------------------
C     
C     -- AUTHOR(S)
C     
C     2005/10/3, Bruno Sportisse, CEREA.
C     
C------------------------------------------------------------------------

      IMPLICIT NONE

      include 'CONST.INC'
      include 'aerpar_simple.inc'
      include 'droppar_simple.inc'

      INTEGER NS
      double precision dsf_aero(NS),dbf_aero(NS+1)
      double precision fdistx(NS), fdistx2(NS)
      double precision sumx,sumx1,sumx2
      double precision fracx(NS),fracx1(NS), fracx2(NS)
      double precision rac2pi1,rac2pi2,logsd1,logsd2,den1,den2
      double precision delD
      integer i

      integer ifirstact

      logsd1 = DLOG(SD1_aq)
      logsd2 = DLOG(SD2_aq)
      rac2pi1 = DSQRT(2.D0*PI)*logsd1
      rac2pi2 = DSQRT(2.D0*PI)*logsd2
      den1 = 2.D0*logsd1**2.D0
      den2 = 2.D0*logsd2**2.D0

      do i=1,ifirstact-1
         fdistx(i)=0.d0
         fdistx2(i)=0.d0
      enddo

      do i = ifirstact,NS
         delD = dbf_aero(i+1)-dbf_aero(i)
         fracx1(i) = (delD/(dsf_aero(i)*rac2pi1))
     &        *DEXP(-(DLOG(dsf_aero(i)/DM1_aq))**2.d0/den1)
         fracx2(i) = (delD/(dsf_aero(i)*rac2pi2))
     &        *DEXP(-(DLOG(dsf_aero(i)/DM2_aq))**2.d0/den2)
         fracx(i) = fracx1(i) + fracx2(i)
      enddo

      sumx = 0.d0
      sumx1 = 0.d0
      sumx2 = 0.d0
      
c     Normalize so fdist sums to 1      
      do i = ifirstact,NS
         sumx = sumx + fracx(i)
         if(dsf_aero(i) .lt. dsep) then
            sumx1 = sumx1+fracx1(i)
         else
            sumx2 = sumx2+fracx2(i)
         endif
      enddo

      do i =ifirstact,NS
         fdistx(i) = fracx(i)/sumx
         if( dsf_aero(i) .lt. dsep) then
            fdistx2(i) = fracx1(i)/sumx1
         else
            fdistx2(i) = fracx2(i)/sumx2
         endif
      enddo

      return
      end

