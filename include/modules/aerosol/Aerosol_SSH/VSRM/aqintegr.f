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

      subroutine aqintegr(y,t,tout)

      implicit none 
      external aqfex,jac
      include 'aerpar.inc'
      include 'droppar.inc'
      include 'num_aq.inc'

      integer lrw, liw, j
      parameter (lrw = 22 +  9*MEQN1 + 2*MEQN1**2)
      parameter (liw = 30 + MEQN1)

      double precision Y(meqn1)
      double precision RWORK(lrw)
      integer IWORK(liw)
      DOUBLE PRECISION ATOL, RPAR, RTOL, T, TOUT
      integer ISTATE, ITASK, IPAR, MF, IOPT, ITOL
      
      MF = 22
      ITOL = 1
      RTOL = 1.D-5
      ATOL = 5.D-4
      ITASK = 1
      ISTATE = 1
      IOPT = 1
      RPAR = 0.D0
      IPAR = 0

C set all values to 0.0
      DO j=1,liw
         IWORK(j)=0
      ENDDO

      DO j=1,lrw
         RWORK(j)=0.D0
      ENDDO
      
C IWORK(6) is MXSTEP, maximum of steps allowed
      IWORK(6)=100000
c     
c     READY FOR THE CALL TO DVODE
c     
      CALL DVODE (aqfex, MEQN1, Y, T, TOUT, ITOL, RTOL, ATOL, ITASK,
     1     ISTATE, IOPT, RWORK, LRW, IWORK, LIW, JAC, MF,
     2     RPAR, IPAR)

      DO j=1,MEQN1
         y(j)=dmax1(y(j),TINYAQ)
      enddo

      end
