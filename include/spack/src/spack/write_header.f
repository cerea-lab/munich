      subroutine write_header
C------------------------------------------------------------------------
C
C     -- DESCRIPTION
C
C     Write headers of subroutines
C
C     FEXCHEM.F    : chemical production terms (P-Lc).
C     JACDCHEMDC.F : Jacobian Matrix (d/dc(P-Lc)).
C     KINETIC.F    : kinetic rate k.
C     RATES.F      : reaction rates w.
C     DRATEDC.F    : derivatives dw/dc.
C     FEXPROD.F    : production terms P.
C     FEXLOSS.F    : loss terms L.
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
C     Label 300: blank line.
C     Label 150: comments.
C     Label 200: C-------...
C------------------------------------------------------------------------
C
C     -- MODIFICATIONS
C
C     2008/12/23: correction of input variables in the generated files
C                (Youngseob KIM)
C     2009/07/16: Modification of headers for SIREAM
C     2010/03/04: Addition of "dimensions.f"
C
C------------------------------------------------------------------------
C
C     -- AUTHOR(S)
C
C     Bruno Sportisse, CEREA, 2003.
C
C------------------------------------------------------------------------
      implicit none
      include 'nficfort'
      include 'parametre'
      integer nwrite

C     Setup of Heterogeneous reactions no.
      integer HR1,HR2,HR3,HR4
      common /HR/ HR1,HR2,HR3,HR4

c---------- for RACM
c      if (chem_mechanism.eq.1) then
      if (mechanism_name .eq. "racm  " .or.
     &    mechanism_name .eq. "radm  ") then
      HR1=238
      HR2=239
      HR3=240
      HR4=241
c      mechanism_name = 'racm'
c---------- for CB05
c      elseif (chem_mechanism.eq.2) then
      elseif (mechanism_name .eq. "cb05  ") then
      HR1=156
      HR2=157
      HR3=158
      HR4=159
c      mechanism_name = 'cb05'
c---------- for RACM2
c      elseif (chem_mechanism.eq.3) then
      elseif (mechanism_name .eq. "racm2 ") then
      HR1=352
      HR2=353
      HR3=354
      HR4=355
c      mechanism_name = 'racm2'
c---------- for MELCHIOR2
c      elseif (chem_mechanism.eq.4) then
      elseif (mechanism_name .eq. "melchior2") then
      HR1=1
      HR2=1
      HR3=1
      HR4=1
c      mechanism_name = 'melchior2'
      else
      write(*,*) "Error:"
      write(*,*) "wrong number of chemical mechanism in parametre"
      endif

C     Routine dimensions.f
      nficdim=ipiste
      ipiste=ipiste+1
      open(nficdim,file='dimensions.f',status='new')
      nwrite=nficdim

      call write_common_header(nwrite)

C     Routine kinetic.f
      nfick=ipiste
      ipiste=ipiste+1
      open(nfick,file='kinetic.f',status='new')
      nwrite=nfick

      call write_common_header(nwrite)

      if (aerosol_formation) then
      write(*,*)
      write(*,*) "###################################################"
      write(*,*) "Option activated: aerosol_formation (see parametre)"
      write(*,*) "###################################################"
      write(*,*)
C=========== aerosol formation =============================
C   Heterogeneous reactions and SVOC formations are included
C   in the gas-phase reactions
C===========

      write(nwrite,401) function_suffix
 401  format(6x,'subroutine kinetic',a,'(')
      write(nwrite,4011)
 4011 format(5x,'$',4x,'Ns,Nbin_aer,nr,IHETER,ICLD,rk,temp,xlw,')
      write(nwrite,4012)
 4012 format(5x,'$',4x,'Press,azi,att,lwctmp,granulo,WetDiam,',
     2     'dsf_aero,ispeclost,')
      write(nwrite,4013)
 4013 format(5x,'$',4x,'Wmol,LWCmin,option_photolysis)')

C===========
C     Only gas-phase reactions
C==========
      else
      write(nwrite,4014) function_suffix
 4014 format(6x,'subroutine kinetic',a,'(')
      write(nwrite,4015)
 4015 format(5x,'$',4x,'nr,rk,temp,xlw,Press,azi,att,')
      write(nwrite,4016)
 4016 format(5x,'$',4x,'option_photolysis)')
      endif
C=====================
      write(nwrite,300)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,402)
 402  format('C',5x,'-- DESCRIPTION')
      write(nwrite,150)
      write(nwrite,403)
 403  format('C',5x,'This routine computes the kinetic rates',
     2     ' for the gas-phase.')
      write(nwrite,404)
 404  format('C     This routine is automatically generated by SPACK.')
      write(nwrite,250)filemeca
      write(nwrite,251)filespecies
 250  format('C     Mechanism: ',a20)
 251  format('C     Species: ',a20)
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,405)
 405  format('C     -- INPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,4061)
 4061 format('C     NR: reaction number.')
      write(nwrite,406)
 406  format('C     TEMP: temperature ([K]).')
      write(nwrite,407)
 407  format('C     XLW: water massic fraction.')
      write(nwrite,408)
 408  format('C     PRESS: pressure ([Pa]).')
      write(nwrite,409)
 409  format('C     AZI: zenithal angle ([degree]).')
      write(nwrite,410)
 410  format('C     ATT: attenuation variable.')
C     aerosol formatios =================
      if (aerosol_formation) then
      write(nwrite,4101)
 4101 format('C     lwctmp: Liquid Water Content ([g.m-3])')
      write(nwrite,4102)
 4102 format('C     WetDiam: Bins wet Diameter ([mu m])')
      write(nwrite,4103)
 4103 format('C     granulo: Particles number in each bin.')
      endif
C ========================================
      write(nwrite,150)
      write(nwrite,413)
 413  format('C     -- INPUT/OUTPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,414)
 414  format('C     -- OUTPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,415)
 415  format('C     RK: kinetic rates.')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,416)
 416  format('C     -- REMARKS')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,417)
 417  format('C     -- MODIFICATIONS')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,418)
 418  format('C     -- AUTHOR(S)')
      write(nwrite,150)
      write(nwrite,419)
 419  format('C     SPACK.')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,300)
      write(nwrite,420)
 420  format('      implicit none')
      write(nwrite,300)
C     aerosol formation ==============
      if (aerosol_formation) then
      write(nwrite,421)
 421  format('      integer Ns,Nbin_aer,nr')
      write(nwrite,422) HR1,HR2,HR3,HR4
 422  format('      double precision rk',I3,',rk',I3,',rk',I3,',rk',I3)
      write(nwrite,423)
 423  format('      double precision lwctmp')
      write(nwrite,424)
 424  format('      double precision WetDiam(Nbin_aer)')
      write(nwrite,425)
 425  format('      double precision granulo(Nbin_aer)')
      write(nwrite,426)
 426  format('      double precision dsf_aero(Nbin_aer)')
      write(nwrite,427)
 427  format('      integer ICLD,IHETER')
      write(nwrite,428)
 428  format('      integer ispeclost(4)')
      write(nwrite,429)
 429  format('      double precision Wmol(Ns),LWCmin')
      else
      write(nwrite,430)
 430  format('      integer nr')
      endif
C ====================================
      write(nwrite,432)
 432  format('      double precision rk(nr),temp,xlw,Press')
      write(nwrite,433)
 433  format('      double precision Effko,Rapk,facteur,SumM,',
     2     'azi,att')
      write(nwrite,434)
 434  format('      double precision YlH2O')
      write(nwrite,4342)
 4342 format('      double precision VO2, VN2, VH2O')
      write(nwrite,4343)
 4343 format('      double precision vd_NO2, h_mixing')
      write(nwrite,435)
 435  format('      integer option_photolysis')
      write(nwrite,300)
      write(nwrite,438)
 438  format('C     Compute third body.')
      write(nwrite,439)
 439  format('C     Conversion = Avogadro*1d-6/Perfect gas constant.')
      write(nwrite,440)
 440  format('C     PRESS in Pascal, SUMM in molecules/cm3,',
     2     ' TEMP in Kelvin')
      write(nwrite,300)
      write(nwrite,441)
 441  format('      SumM = Press * 7.243D16 / temp')
      write(nwrite,300)
      write(nwrite,442)
 442  format('C     Number of water molecules computed from',
     2     ' the massic fraction')
      write(nwrite,443)
 443  format('C     (absolute humidity)')
      write(nwrite,300)
      write(nwrite,444)
 444  format('      YlH2O = 29.d0*SumM*xlw/(18.d0+11.d0*xlw)')
      write(nwrite,300)
      write(nwrite,445)
 445  format('C     For the zenithal angle at tropics')
      write(nwrite,300)
      write(nwrite,446)
 446  format('      azi=abs(azi)')
      write(nwrite,300)

C     Routine fexchem.f
      nficf=ipiste
      ipiste=ipiste+1
      open(nficf,file='fexchem.f',status='new')
      nwrite=nficf

      call write_common_header(nwrite)
C     aerosol formation =====================
      if (aerosol_formation) then
      write(nwrite,501) function_suffix
 501  format(6x,'subroutine fexchem',a,'(')
      write(nwrite,5011)
 5011 format(5x,'$',4x,'NS,Nr,nemis,y,rk,ZCsourc,convers_factor,chem)')
      else
      write(nwrite,5012) function_suffix
 5012 format(6x,'subroutine fexchem',a,'(')
      write(nwrite,5013)
 5013 format(5x,'$',4x,'ns,nr,y,rk,ZCsourc,convers_factor,chem)')
      endif
C ====================================
      write(nwrite,300)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,502)
 502  format('C',5x,'-- DESCRIPTION')
      write(nwrite,150)
      write(nwrite,503)
 503  format('C',5x,'This routine computes the chemical production',
     2     ' term.')
      write(nwrite,504)
 504  format('C     This routine is automatically generated by SPACK.')
      write(nwrite,250)filemeca
      write(nwrite,251)filespecies
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,505)
 505  format('C     -- INPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,5061)
 5061 format('C     Ns: chemical species number.')
      write(nwrite,506)
 506  format('C     Y: chemical concentrations.')
      write(nwrite,415)
C     415  format('C     RK: kinetic rates.')
 507  format('C     W: reaction rates.')
      write(nwrite,508)
 508  format('C     ZCSOURC: volumic emissions.')
      write(nwrite,150)
      write(nwrite,513)
 513  format('C     -- INPUT/OUTPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,514)
 514  format('C     -- OUTPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,515)
 515  format('C     CHEM: array of chemical production terms.')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,516)
 516  format('C     -- REMARKS')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,517)
 517  format('C     -- MODIFICATIONS')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,518)
 518  format('C     -- AUTHOR(S)')
      write(nwrite,150)
      write(nwrite,519)
 519  format('C     SPACK.')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,300)
      write(nwrite,520)
 520  format('      implicit none')
      write(nwrite,300)
C     aerosol formation ===============
      if (aerosol_formation) then
      write(nwrite,532)
 532  format('      integer nemis')
      endif
C ===================================
      write(nwrite,535)
 535  format('      integer nr,ns,i')
      write(nwrite,536)
 536  format('      double precision w(nr),rk(nr),y(ns),chem(ns),',
     2     'ZCsourc(ns)')
      write(nwrite,537)
 537  format('      double precision conc(ns), convers_factor(ns)')
      write(nwrite,300)
      write(nwrite,538)
 538  format('      do i=1,ns')
      write(nwrite,539)
 539  format('       chem(i)=0.d0')
      write(nwrite,540)
 540  format('      enddo')
      write(nwrite,300)
      write(nwrite,541)
 541  format('C     Conversion mug/m3 to molecules/cm3.')
      write(nwrite,300)
      write(nwrite,542)
 542  format('      do i = 1, ns')
      write(nwrite,543)
 543  format('         conc(i) = y(i) * convers_factor(i)')
      write(nwrite,544)
 544  format('      enddo')
      write(nwrite,300)
      write(nwrite,545)
 545  format('C     Compute reaction rates.')
      write(nwrite,300)
      write(nwrite,546) function_suffix
 546  format('      call rates',a,'(')
      write(nwrite,5461)
 5461 format(5x,'$',4x,'ns,nr,rk,conc,w)')
      write(nwrite,300)
      write(nwrite,547)
 547  format('C     Chemical production terms.')
      write(nwrite,300)

C     Routine jacdchemdc.f
      nficj=ipiste
      ipiste=ipiste+1
      open(nficj,file='jacdchemdc.f',status='new')
      nwrite=nficj

      call write_common_header(nwrite)

      write(nwrite,601) function_suffix
 601  format(6x,'subroutine jacdchemdc',a,'(')
      write(nwrite,6011)
 6011 format(5x,'$',4x,
     2       'ns,nr,y,convers_factor,convers_factor_jac,rk,JacC)')
      write(nwrite,300)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,602)
 602  format('C',5x,'-- DESCRIPTION')
      write(nwrite,150)
      write(nwrite,603)
 603  format('C',5x,'This routine computes the Jacobian matrix',
     2     ' for the gas-phase.')
      write(nwrite,604)
 604  format('C     This routine is automatically generated by SPACK.')
      write(nwrite,250)filemeca
      write(nwrite,251)filespecies
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,605)
 605  format('C     -- INPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,5061)
c 5061 format('C     NS: chemical species number.')
      write(nwrite,4061)
c 4061 format('C     NR: reaction number.')
      write(nwrite,606)
 606  format('C     Y: chemical concentrations.')
      write(nwrite,415)
C     415  format('C     RK: kinetic rates.')
 607  format('C     DW: derivative of reaction rates wrt Y.')
      write(nwrite,611)
 611  format('C     convers_factor: unity conversion factor of ',
     2     'mug/m3 to molecules/cm3.')
      write(nwrite,612)
 612  format('C     convers_factor_jac: Wmol(i)/Wmol(j) for JacC(i,j)')
      write(nwrite,150)
      write(nwrite,613)
 613  format('C     -- INPUT/OUTPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,614)
 614  format('C     -- OUTPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,615)
 615  format('C     JACC: Jacobian matrix.')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,616)
 616  format('C     -- REMARKS')
      write(nwrite,150)
      write(nwrite,609)
 609  format('C     The matrix JACC could be stored in a ',
     2     'low-dimensional vector.')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,617)
 617  format('C     -- MODIFICATIONS')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,618)
 618  format('C     -- AUTHOR(S)')
      write(nwrite,150)
      write(nwrite,619)
 619  format('C     SPACK.')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,300)
      write(nwrite,620)
 620  format('      implicit none')
      write(nwrite,300)
      write(nwrite,633)
 633  format('      integer i,j,nr,ns')
      write(nwrite,634)
 634  format('      double precision dw(nr,ns),rk(nr),y(ns),',
     2     'JacC(ns,ns)')
      write(nwrite,6341)
 6341 format('      double precision convers_factor(ns)')
      write(nwrite,6342)
 6342 format('      double precision convers_factor_jac(ns,ns)')
      write(nwrite,6343)
 6343 format('      double precision conc(ns)')
      write(nwrite,300)
      write(nwrite,635)
 635  format('      do j=1,ns')
      write(nwrite,636)
 636  format('       do i=1,ns')
      write(nwrite,637)
 637  format('        JacC(i,j)=0.d0')
      write(nwrite,638)
 638  format('       enddo')
      write(nwrite,639)
 639  format('      enddo')
      write(nwrite,300)
      write(nwrite,541)
c 541 format('C     Conversion mcg/m3 to molecules/cm3.')
      write(nwrite,300)
      write(nwrite,542)
c 542 format('      do i = 1, ns')
      write(nwrite,543)
c 543 format('         conc(i) = y(i) * convers_factor(i)')
      write(nwrite,544)
c 544 format('      enddo')
      write(nwrite,300)
      write(nwrite,640) function_suffix
 640  format('      call dratedc',a,'(')
      write(nwrite,6401)
 6401 format(5x,'$',4x,'ns,nr,rk,conc,dw)')
      write(nwrite,300)

C==========================================
C     Routine fexloss.f
      nficloss=ipiste
      ipiste=ipiste+1
      open(nficloss,file='fexloss.f',status='new')
      nwrite=nficloss

      call write_common_header(nwrite)

      write(nwrite,701) function_suffix
 701  format(6x,'subroutine fexloss',a,'(')
      write(nwrite,7011)
 7011 format(5x,'$',4x,'NS,dw,loss)')
      write(nwrite,300)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,502)
C     502  format('C',5x,'-- DESCRIPTION')
      write(nwrite,150)
      write(nwrite,703)
 703  format('C',5x,'This routine computes the chemical loss ',
     2     ' term L in a P-Lc formulation.')
      write(nwrite,504)
C     504  format('C     This routine is automatically generated by SPACK.')
      write(nwrite,250)filemeca
      write(nwrite,251)filespecies
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,505)
C     505  format('C     -- INPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,607)
C     607  format('C     DW: derivative of reaction rates wrt Y.')
      write(nwrite,150)
      write(nwrite,513)
C     513  format('C     -- INPUT/OUTPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,514)
C     514  format('C     -- OUTPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,715)
 715  format('C     LOSS: array of chemical loss terms.')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,516)
C     516  format('C     -- REMARKS')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,517)
C     517  format('C     -- MODIFICATIONS')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,518)
C     518  format('C     -- AUTHOR(S)')
      write(nwrite,150)
      write(nwrite,519)
C     519  format('C     SPACK.')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,300)
      write(nwrite,520)
C     520  format('      implicit none')
      write(nwrite,300)
      write(nwrite,737)
 737  format('      integer i,nr,NS')
      write(nwrite,733)
 733  format('      double precision dw(nr,NS),loss(NS)')
      write(nwrite,300)
      write(nwrite,738)
 738  format('      do i=1,NS')
      write(nwrite,739)
 739  format('       loss(i)=0.d0')
      write(nwrite,740)
 740  format('      enddo')
      write(nwrite,300)
      write(nwrite,741)
 741  format('C     Chemical loss terms.')
      write(nwrite,300)

C===============================================
C     Routine fexprod.f
      nficprod=ipiste
      ipiste=ipiste+1
      open(nficprod,file='fexprod.f',status='new')
      nwrite=nficprod

      call write_common_header(nwrite)

      write(nwrite,801) function_suffix
 801  format(6x,'subroutine fexprod',a,'(')
      write(nwrite,8011)
 8011 format(5x,'$',4x,'NS,w,prod)')
      write(nwrite,300)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,502)
C     502  format('C',5x,'-- DESCRIPTION')
      write(nwrite,150)
      write(nwrite,803)
 803  format('C',5x,'This routine computes the production ',
     2     ' term P in a P-Lc formulation.')
      write(nwrite,504)
C     504  format('C     This routine is automatically generated by SPACK.')
      write(nwrite,250)filemeca
      write(nwrite,251)filespecies
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,505)
C     505  format('C     -- INPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,507)
C     507  format('C     W: reaction rates.')
      write(nwrite,150)
      write(nwrite,513)
C     513  format('C     -- INPUT/OUTPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,514)
C     514  format('C     -- OUTPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,815)
 815  format('C     PROD: array of chemical production terms.')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,516)
C     516  format('C     -- REMARKS')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,517)
C     517  format('C     -- MODIFICATIONS')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,518)
C     518  format('C     -- AUTHOR(S)')
      write(nwrite,150)
      write(nwrite,519)
C     519  format('C     SPACK.')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,300)
      write(nwrite,520)
C     520  format('      implicit none')
      write(nwrite,300)
      write(nwrite,737)
C     737  format('      integer i,nr,NS')
      write(nwrite,833)
 833  format('      double precision w(nr),prod(NS)')
      write(nwrite,300)
      write(nwrite,738)
C     538  format('      do i=1,NS')
      write(nwrite,839)
 839  format('       prod(i)=0.d0')
      write(nwrite,740)
C     540  format('      enddo')
      write(nwrite,300)
      write(nwrite,841)
 841  format('C     Chemical production terms.')
      write(nwrite,300)

C====================================================
C     Routine rates.f
      nficw=ipiste
      ipiste=ipiste+1
      open(nficw,file='rates.f',status='new')
      nwrite=nficw

      call write_common_header(nwrite)

      write(nwrite,901) function_suffix
 901  format(6x,'subroutine rates',a,'(')
      write(nwrite,9011)
 9011 format(5x,'$',4x,'ns,nr,rk,y,w)')
      write(nwrite,300)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,502)
C     502  format('C',5x,'-- DESCRIPTION')
      write(nwrite,150)
      write(nwrite,902)
 902  format('C',5x,'This routine computes the reaction rates.')
      write(nwrite,504)
C     504  format('C     This routine is automatically generated by SPACK.')
      write(nwrite,250)filemeca
      write(nwrite,251)filespecies
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,505)
C     505  format('C     -- INPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,5061)
c 5061 format('C     NS: chemical species number.')
      write(nwrite,4061)
c 4061 format('C     NR: reaction number.')
      write(nwrite,903)
 903  format('C     RK: kinetic rates.')
      write(nwrite,506)
C     506  format('C     Y: chemical concentrations.')
      write(nwrite,150)
      write(nwrite,513)
C     513  format('C     -- INPUT/OUTPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,514)
C     514  format('C     -- OUTPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,507)
C     507  format('C     W: reaction rates.')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,516)
C     516  format('C     -- REMARKS')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,517)
C     517  format('C     -- MODIFICATIONS')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,518)
C     518  format('C     -- AUTHOR(S)')
      write(nwrite,150)
      write(nwrite,519)
C     519  format('C     SPACK.')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,300)
      write(nwrite,520)
C     520  format('      implicit none')
      write(nwrite,300)
      write(nwrite,905)
 905  format('      integer nr,ns')
      write(nwrite,904)
 904  format('      double precision rk(nr),y(ns),w(nr)')
      write(nwrite,300)
      write(nwrite,300)
      write(nwrite,300)

C     Routine dratedc.f
      nficdw=ipiste
      ipiste=ipiste+1
      open(nficdw,file='dratedc.f',status='new')
      nwrite=nficdw

      call write_common_header(nwrite)

      write(nwrite,910) function_suffix
 910  format(6x,'subroutine dratedc',a,'(')
      write(nwrite,9101)
 9101 format(5x,'$',4x,'ns,nr,rk,y,dw)')
      write(nwrite,300)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,502)
C     502  format('C',5x,'-- DESCRIPTION')
      write(nwrite,150)
      write(nwrite,911)
 911  format('C',5x,'This routine computes the derivative of reaction ',
     2     ' rates.')
      write(nwrite,504)
C     504  format('C     This routine is automatically generated by SPACK.')
      write(nwrite,250)filemeca
      write(nwrite,251)filespecies
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,505)
C     505  format('C     -- INPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,5061)
c 5061 format('C     NS: chemical species number.')
      write(nwrite,4061)
c 4061 format('C     NR: reaction number.')
      write(nwrite,903)
C     903  format('C     RK: kinetic rates.')
      write(nwrite,506)
C     506  format('C     Y: chemical concentrations.')
      write(nwrite,150)
      write(nwrite,513)
C     513  format('C     -- INPUT/OUTPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,514)
C     514  format('C     -- OUTPUT VARIABLES')
      write(nwrite,150)
      write(nwrite,607)
C     607  format('C     DW: derivative of reaction rates wrt Y.')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,516)
C     516  format('C     -- REMARKS')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,517)
C     517  format('C     -- MODIFICATIONS')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,150)
      write(nwrite,518)
C     518  format('C     -- AUTHOR(S)')
      write(nwrite,150)
      write(nwrite,519)
C     519  format('C     SPACK.')
      write(nwrite,150)
      write(nwrite,200)
      write(nwrite,300)
      write(nwrite,520)
C     520  format('      implicit none')
      write(nwrite,300)
      write(nwrite,905)
c 905 format('      integer nr,ns')
      write(nwrite,912)
 912  format('      double precision rk(nr),y(ns),dw(nr,ns)')
      write(nwrite,300)
      write(nwrite,300)
      write(nwrite,300)

C     File non_zero.dat
      nficnz=ipiste
      ipiste=ipiste+1
      open(nficnz,file='non_zero.dat',status='new')
C
 150  format('C',6x,a65)
 200  format('C-----------------------------------',
     2     '-------------------------------------')
 300  format(' ')

      RETURN
      END
