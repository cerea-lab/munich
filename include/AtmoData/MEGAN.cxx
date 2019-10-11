// Copyright (C) 2003-2007, ENPC - INRIA - EDF R&D
// Author(s): Vivien Mallet
//
// This file is part of AtmoData library, a tool for data processing in
// atmospheric sciences.
//
// AtmoData is developed in the INRIA - ENPC joint project-team CLIME and in
// the ENPC - EDF R&D joint laboratory CEREA.
//
// AtmoData is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// AtmoData is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
// FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
// details.
//
// For more information, visit the AtmoData home page:
//      http://cerea.enpc.fr/polyphemus/atmodata.html


#ifndef MEGAN_FILE_EMISSIONS_CXX

namespace AtmoData
{
  void defined_species(string species)
  {
	if(species!="ISOP" && species!="MBO" && species!="FORM" &&
	   species!="CO" && species!="MYRC" && species!="SABI" &&
	   species!="LIMO" && species!="3CAR" && species!="OCIM" &&
	   species!="BPIN" && species!="APIN" && species!="OMTP" &&
	   species!="FARN" && species!="BCAR" && species!="OSQT" &&
	   species!="MEOH" && species!="ACTO" && species!="NO" &&
	   species!="CH4" && species!="ACTA")
	  throw string("ERROR: Species ")
		+ species
		+ string(" not defined in MEGAN.");
  }
  
  
  template <class real>
  void GAMMA_T_ISOP(const real& temperature, const real& daily_temperature, real& gamma)
  {
	real CT1,CT2,Eopt,Topt,numerator,divisor,x;
	CT1=80.0;
	CT2=200;
	Eopt=1.75*exp(0.08*(daily_temperature-297));
	Topt=313+0.6*(daily_temperature-297);
	x=((1/Topt)-(1/temperature))/0.00831;
	numerator=Eopt*CT2*exp(CT1*x);
	divisor=CT2-CT1*(1-exp(CT2*x));
	gamma=numerator/divisor;

	/*const double alpha(0.0027), c_l1(1.066);
    const double Ts(303), Tm(314), R(8.314), c_t1(95000), c_t2(230000);
    const double alpha2(alpha * alpha), ratio1(c_t1 / (R * Ts)),
      ratio2(c_t2 / (R * Ts));

	real T = temperature;
	
	gamma = exp(ratio1 * (T-Ts) / T)
	/ (1. + exp(ratio2 * (T-Tm) / T));*/
  }

  template <class real>
  void GAMMA_T_non_ISOP(string species, real& temperature, real& gamma)
  {
	real temperature_factor,Ts;

	Ts=303.0;
	
	if(species=="ISOP" || species=="MBO" || species=="FORM" || species=="CO")
	  temperature_factor=0.09;
	else if(species=="MYRC" || species=="SABI" || species=="LIMO"
			|| species == "3CAR" || species == "OCIM" || species=="BPIN"
			|| species=="APIN" || species == "OMTP")
	  temperature_factor=0.1;
	else if(species=="FARN" || species=="BCAR" || species=="OSQT")
	  temperature_factor=0.17;
	else if(species=="MEOH")
	  temperature_factor=0.08;
	else if(species=="ACTO" || species=="NO")
	  temperature_factor=0.11;
	else if(species=="CH4")
	  temperature_factor=0.05;
	else if(species=="ACTA")
	  temperature_factor=0.13;
	else
	  throw string("ERROR: Species ")
		+ species
		+ string(" temperature factor not defined.");
	
	gamma=exp(temperature_factor*(temperature-Ts));
  }

  template <class real>
  void GAMMA_LAI(real& LAI, real& gamma)
  {
	gamma=0.49*LAI/(pow((1.0+0.2*pow(LAI,2)),0.5));
  }
  
  template <class real>
  void solarangle(int iday, real& localhour, real& lat, real& sinangle)
  {
	real pi,sindelta,cosdelta,A,B;
	pi = 3.14159;
	sindelta = -sin(0.40907)*cos(6.28*float((iday+10))/365.0);
	cosdelta = pow((1-pow(sindelta,2)),0.5);
	A = sin(lat*pi/180.0)*sindelta;
	B = cos(lat*pi/180.0)*cosdelta;
	sinangle=A+B*cos(2.0*pi*(localhour-12.0)/24.0);
  }

  template <class real>
  void GAMMA_LIGHT(int iday, real& sinangle, real& PAR,
				   real& daily_PAR, real& gamma)
  {
	real A,B,top_PAR, phi;

	if (sinangle < 0.0)
	  gamma=0.0;
	else
	  {
		top_PAR=3000.0+99.0*cos(2.0*3.14-float((iday-10))/365.0);
		phi=PAR/(sinangle*top_PAR);
		B=1.0+0.0005*(daily_PAR-400.0);
		A=2.46*B*phi-0.9*pow(phi,2);
		gamma=sinangle*A; 
		
		/*real alpha=0.0027;
		real beta=1.066;
		gamma=alpha*beta*PAR/(pow((1.0+pow(alpha*PAR,2)),0.5));*/

		if (gamma < 0.0)
		  gamma=0.0;
	  }
  }
  
  template <class real>
  void GAMMA_AGE(string species, real& LAIp, real& LAIc,
				 real& daily_temperature, real& t, real& gamma)
  {
	real Fnew,Fgro,Fmat,Fold;
	real ti,tm;
	real Anew,Agro,Amat,Aold;

	if(species=="ISOP" || species=="MBO")
	  {
		Anew=0.05;
		Agro=0.6;
		Amat=1.125;
		Aold=1.0;
	  }
	else if(species=="ACTO" || species=="ACTA" || species=="FORM"
			|| species=="CH4" || species=="NO" || species=="CO")
	  {
		Anew=1.0;
		Agro=1.0;
		Amat=1.0;
		Aold=1.0;
	  }
	else if(species=="MYRC" || species=="SABI" || species=="LIMO"
			|| species =="3CAR" || species=="OCIM" || species=="BPIN"
			|| species=="APIN" || species == "OMTP")
	  {
		Anew=2.0;
		Agro=1.8;
		Amat=0.95;
		Aold=1.0;
	  }
	else if(species=="FARN" || species =="BCAR" || species =="OSQT")
	  {
		Anew=0.4;
		Agro=0.6;
		Amat=1.075;
		Aold=1.0;
	  }
	else if(species=="MEOH")
	  {
		Anew=3.0;
		Agro=2.6;
		Amat=0.85;
		Aold=1.0;
	  }
	else
	  throw string("ERROR: Species ")
		+ species
		+ string(" canopies factors not defined.");
	
	if (daily_temperature > 303.0)
	  ti=2.9;
	else
	  ti=5.0+(0.7*(300.0-daily_temperature));
	tm=2.3*ti;

	if(LAIp<LAIc)
	  {
		if(t<=ti)
		  Fnew=1-LAIp/LAIc;
		else
		  Fnew=(ti/t)*(1-LAIp/LAIc);
		if(t<=tm)
		  Fmat=LAIp/LAIc;
		else
		  Fmat=LAIp/LAIc+(t-tm)/t*(1-LAIp/LAIc);
		Fgro=1-Fnew-Fmat;
		Fold=0.0;
	  }
	else if (LAIp==LAIc)
	  {
		Fmat=0.8;
		Fnew=0.0;
		Fgro=0.1;
		Fold=0.1;
	  }
	else
	  {
		Fnew=0.0;
		Fgro=0.0;
		Fold=(LAIp-LAIc)/LAIp;
		Fmat=1-Fold;
	  }
	gamma=Fold*Aold+Fnew*Anew+Fgro*Agro+Fmat*Amat;
  }

  template <class real>
  void GAMMA_SM(real& WP, real& SoilMoisture, real& gamma)
  {
	real deltaSM=0.06;
	real thetal;
	thetal=WP+deltaSM;
	if (SoilMoisture > thetal)
	  gamma=1.0;
	else if (SoilMoisture < WP)
	  gamma=0.0;
	else
	  gamma=(SoilMoisture-WP)/deltaSM;
  }

  template <class real>
  void LightDependentFactor(string species, real& LDF)
  {
	if (species=="ISOP" || species=="MBO")
	  LDF=0.9999;
	else if (species=="MYRC" || species=="LIMO" || species =="3CAR")
	  LDF=0.05;
	else if (species=="SABI" || species=="BPIN" || species == "APIN" ||
			 species=="OMTP")
	  LDF=0.1;
	else if (species=="OCIM")
	  LDF=0.8;
	else if (species=="MEOH" || species=="CH4")
	  LDF=0.75;
	else if (species=="ACTO")
	  LDF=0.25;
	else if (species=="NO")
	  LDF=0.0;
	else if (species=="ACTA" || species=="FARN" || species=="BCAR" ||
			 species =="OSQT" || species=="FORM" || species=="CO")
	  LDF=0.5;
	else
	  throw string("ERROR: Species ")
		+ species
		+ string(" Light Dependent Factor not defined.");
  }

  template <class real>
  void aggregation_MEGAN(const Data<real, 4>& emis_out,
  						 string file_aggregation,
  						 Data<real, 4>& emissions,
  						 vector<string>& biogenic_names,
  						 vector<string>& output_names)
  {
	int Nsp_in = emis_out.GetLength(0);
	int Ny = emis_out.GetLength(1);
	int Nx = emis_out.GetLength(2);
	int Nt = emis_out.GetLength(3);
	int Nsp_out = emissions.GetLength(0);
	int Nmegan = 20;

	string name_in,name_out;
	RegularGrid<int> GridInSpecies(Nsp_in);
	Data<int, 1> indices(GridInSpecies);

	int i,j,k,x,y,t;
	real coefficient;

	emissions.Fill(0.0);

	ExtStream aggregation_stream(file_aggregation);
    if (!aggregation_stream.is_open())
      throw string(" File ") + file_aggregation + " doesn't exist.";

	while (!aggregation_stream.IsEmpty())
	  {
		aggregation_stream >> name_out;
		if (name_out=="SPECIES")
		  for(i=0;i<Nmegan;i++)
			{
			  aggregation_stream >> name_in;
			  for(j=0;j<Nsp_in;j++)
				if(name_in==biogenic_names[j])
					indices(j)=i;
			}
		else
		  {
			for(i=0;i<Nsp_out;i++)
			  if(name_out==output_names[i])
				for(j=0;j<Nmegan;j++)
				  {
					aggregation_stream >> coefficient;
					if(coefficient!=0.0)
					  for(k=0;k<Nsp_in;k++)
						if(indices(k)==j)
						  for(x=0;x<Nx;x++)
							for(y=0;y<Ny;y++)
							  for(t=0;t<Nt;t++)
								emissions(i,t,y,x)+=coefficient
								  *emis_out(k,y,x,t);
						  
				  }
			
		  }
		
	  }
  }
  
  
		
  
}  // namespace AtmoData.

#define MEGAN_FILE_EMISSIONS_CXX
#endif
