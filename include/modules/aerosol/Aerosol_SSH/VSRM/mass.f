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

      subroutine mass(wl,radius,temp,pres,gcon,con,c,akeq,akhen,fgl,flg,
     &     gfgl,gflg)

C------------------------------------------------------------------------
C     
C     -- DESCRIPTION 
C     
C     This routine computes the mass flux for Henry's reactions between
C     gas-phase and aqueous-phase.
C------------------------------------------------------------------------
C     
C     -- INPUT VARIABLES
C     
C     WL     : liquid water content.
C     RADIUS : droplet radius ([m]).
C     TEMP   : temperature    ([K]).
C     PRES   : pressure       ([atm]).
C     GCON   : gas-phase concentration. 
C     CON    : aqueous-phase concentration.
C     C      : concentration vector of all species.
C     AKEQ   : equilibrium rates.
C     AKHEN  : Henry's rates.
C     
C     -- INPUT/OUTPUT VARIABLES
C     
C     -- OUTPUT VARIABLES
C     
C     F*, G*: flux for mass transfer.
C------------------------------------------------------------------------
C     
C     -- REMARKS
C     
C------------------------------------------------------------------------
C     
C     -- MODIFICATIONS
C     
C     1) Rename N in FUCHS.
C     2) Optimize computation with RTLOC.
C     3) Use function for free mean path.
C     4) Add pres in the inputs.
C------------------------------------------------------------------------
C     
C     -- AUTHOR(S)
C     
C     2005/10/3, Bruno Sportisse, CEREA.
C     
C------------------------------------------------------------------------

      IMPLICIT NONE

      include 'aerpar.inc'
      include 'droppar.inc'

      double precision kn,fuchs,ikn,kmt
      double precision temp,pres,prespa
      
      double precision gcon(21),con(28),akeq(17)
      double precision akhen(21),fgl(21),flg(21)
      double precision c(46),gfgl(21),gflg(21)
      double precision ekhen(21)
      double precision radius,rtloc,wl

      double precision acc, dg
      
      double precision AIR_FREE_MEAN_PATH,VIM

      integer i
      
c     EKHEN(I) IS THE EFFECTIVE HENRY'S LAW CONSTANT

      ekhen(1)=akhen(1)*(1.d0+akeq(1)/c(46)+akeq(1)*
     &     akeq(2)/c(46)**2.d0)
      ekhen(2)=1.0d30
      ekhen(3)=akhen(3)*(1.d0+akeq(7)/c(46))
      ekhen(4)=akhen(4)*(1.d0+akeq(6)/c(46))
      ekhen(5)=akhen(5)*(1.d0+akeq(8)/c(46)+akeq(8)*
     &     akeq(9)/c(46)**2.d0)
      ekhen(6)=akhen(6)*(1.d0+akeq(5)/c(46))
      ekhen(7)=akhen(7)*((1.d0+akeq(12))/akeq(12))
      ekhen(8)=akhen(8)*(1.d0+akeq(13)/c(46))
      ekhen(9)=akhen(9)
      ekhen(10)=akhen(10)
      ekhen(11)=akhen(11)
      ekhen(12)=akhen(12)
      ekhen(13)=akhen(13)
      ekhen(14)=akhen(14)
      ekhen(15)=akhen(15)*(1.d0+akeq(14)/c(46)+(akeq(14)*c(37))/
     &     (akeq(16)*c(46)))
      ekhen(16)=akhen(16)
      ekhen(17)=akhen(17)*(1.d0+akeq(15)/c(46))
      ekhen(18)=akhen(18)
      ekhen(19)=akhen(19)*(1.d0+akeq(10)/c(45))
      ekhen(20)=akhen(20)
      ekhen(21)=akhen(21)

c     COMPUTE_AIR_FREE_MEAN_PATH IS THE MEAN FREE PATH OF AIR ([\mu m]). 
c     KN IS THE KNUDSEN NUMBER

      prespa=pres*101325.d0       
      call COMPUTE_AIR_FREE_MEAN_PATH(TEMP,
     &     prespa,air_free_mean_path,VIM)
      kn=air_free_mean_path*1.D-6/radius

      ikn=1.d0/kn

c     ACC IS THE ACCOMODATION COEFFICIENT ASSUMED THE SAME HERE FOR
c     ALL THE SPECIES

      acc=ACC_AQ

c     N IS THE COEFFICIENT ENTERING THE FLUX EXPRESSION

      fuchs =1.d0/(1.d0+((1.33d0+0.71d0*ikn)/(1.d0+ikn)+4.d0*(1.d0-acc)
     &     /(3.d0*acc))*kn)

c     DG IS THE GAS PHASE DIFFUSIVITY ASSUMED HERE THE SAME FOR ALL
c     THE GASES. WE SHALL PROBABLY HAVE TO CHANGE IT LATER.
c     DG=1.x10-5 m**2/sec

      dg=DG_AQ

c     THE GAS CONSTANT =0.082 (LT.ATM/MOL K)

      rtloc = 0.082058d0*temp
      kmt=(3.0d0*fuchs*dg)/(radius*radius)

      do i=1,21
         fgl(i)=kmt*gcon(i)
         flg(i)=(kmt*con(i))/(ekhen(i)*rtloc)
         gfgl(i)=fgl(i)*wl
         gflg(i)=flg(i)*wl
      enddo

      return
      end
