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

        subroutine redistaq(NS,DBF_AERO,XBF_AERO,
     &      fixed_rho_aero,dnew,aerosol) !! ajout ZNA

c       IMPLICIT NONE

        include 'aerpar.inc'
        include 'droppar.inc'
        include 'paraerochem.inc'
        INCLUDE 'CONST.INC'
        INCLUDE 'CONST_A.INC'

        INTEGER NS
        double precision aerosol(NS,naers)
        double precision aerosolnew(NS,naers)
******  
        INTEGER jj,j1lo,j1hi,js1,js2
        INTEGER isp, i
        double precision dbnew(NS+1),dnew(NS)
        double precision x2lo,x2hi,xmin,xmax,frac
        double precision HSDaq(NS),XBDaq(NS+1)
        double precision aeroorig(naers),aeronew(naers)
        double precision DBF_AERO(NS+1),XBF_AERO(NS+1)
        double precision fixed_rho_aero
        
C       Modif: 
C       1) update DACTIV in microm now.
C       2) replace nsect by NS.
******  zero init

        DO isp = 1,naers
           DO jj=1,NS
            aerosolnew(jj,isp)=0.D0
           END DO
        ENDDO

        dbnew(1) = DBF_AERO(1)

        do i = 1,NS-1
         if (dnew(i) .lt. dactiv) then
           dbnew(i+1) = DBF_AERO(i+1)
         else
           dbnew(i+1) = (dnew(i)*dnew(i+1))**0.5d0
         endif
        enddo 
        if (dbnew(NS) .lt. DBF_AERO(NS+1)) then
         dbnew(NS+1) = DBF_AERO(NS+1)
        else
         dbnew(NS+1) = dbnew(NS)*(dbnew(NS)
     &          /dbnew(NS-1))
       endif


	DO i = 1,NS
	   IF (dbnew(i) .gt. dbnew(i+1)) THEN
	      write(6,*) 'section cross'
	   ENDIF
	ENDDO
	
	DO isp = 1,naers
	   aeroorig(isp) = 0.d0
	ENDDO
	
	DO isp = 1,naers
	   DO jj=1,NS
	      aeroorig(isp)=aeroorig(isp) + aerosol(jj,isp)
	   END DO
	ENDDO

	do i = 1,NS+1
	   XBDaq(i) = dlog(cst_pi6*fixed_rho_aero*(dbnew(i)**3.d0))
	enddo
      	
	
	do i=1,NS
	   HSDaq(i) =XBDaq(i+1)-XBDaq(i)
	enddo
	

******  compute redistribution
	x2lo=XBDaq(1)
	
	CALL LOCATE(NS+1,XBF_AERO,x2lo,j1lo)
	
	DO js2=1,NS
	   x2hi=XBDaq(js2+1)
	   CALL LOCATE(NS+1,XBF_AERO,x2hi,j1hi)

                                ! redistribute over fixed sections
	   DO js1=MAX0(j1lo,1),MIN0(j1hi,NS)
	      
	      xmin=DMAX1(x2lo,XBF_AERO(js1))
	      xmax=DMIN1(x2hi,XBF_AERO(js1+1))

	      
	      frac=(xmax-xmin)/HSDaq(js2)
	      

	      DO isp = 1,naers
		 aerosolnew(js1,isp)=aerosolnew(js1,isp)+
     &      aerosol(js2,isp)*frac
	      ENDDO
	   ENDDO

	   
	   x2lo=x2hi
	   j1lo=j1hi
	END DO
	
******  turn back to conc vector
	
	DO isp = 1,naers
	   DO jj=1,NS
	      aerosol(jj,isp)=aerosolnew(jj,isp)
	   END DO
	END DO
	
******  if diameters went above DMAX, and there is a loss of mass outside set distribution, put mass into last bin
	DO isp = 1,naers
	   aeronew(isp) = 0.d0
	ENDDO
	
	DO isp = 1,naers
	   DO jj=1,NS
	      aeronew(isp)=aeronew(isp) + aerosol(jj,isp)
	   END DO
	ENDDO
	
	DO isp = 1,naers
	   IF (aeronew(isp) .LT. aeroorig(isp)) then
	      aerosol(NS,isp) = aerosol(NS,isp) + 
     &        (aeroorig(isp)-aeronew(isp))
	   ENDIF
	ENDDO

**************************************************
	END
**************************************************
************************************************
