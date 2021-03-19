#include "libRs.hxx"


// Renvoi la résistance de surface Rs (modèle de Zhang).
/*
  \param DynamicViscosity Viscosité dyanmique du milieu (air généralement)
    dans lequel se dépose la particule.
  \param Pressure Pression atmosphérique.
  \param Temperature Température du milieu (air).
  \param ParticleDiameter Diamètre de la particule qui se dépose.
  \param ParticleDensity Densité de la particule qui se dépose.
  \param AirDensity Densité du milieu (air).
  \param FrictionVelocity Vitesse de friction.
  \param CollectorRadius Rayon de collection des éléments du sol
    (valeur tabulée dans le modèle de Zhang).
  \param Alpha Constante du modèle de Zhang
    (dépend du LUC et de la saison SC).
  \param Beta Constante du modèle de Zhang
    (dépend du LUC et de la saison SC).
  \param Gamma Constante du modèle de Zhang
    (dépend du LUC et de la saison SC).

  \return La résistance de surface.
*/
float ComputeSurfaceResistance(float DynamicViscosity, float Pressure,
			       float Temperature, float ParticleDiameter,
			       float ParticleDensity, float AirDensity,
			       float FrictionVelocity, float CollectorRadius,
			       float Alpha, float Beta, float Gamma)
{
  // The resistance model for dry deposition postulates that adjacent
  // to the surface exists a quasi-laminar layer.

  // To represent this resistance for particles, an expression developped
  // by Zhang et al. (2001) has been used.
  
  /////////// Collection efficiency coefficient ///////////
  // Impaction
  float Eim = ComputeImpactionEfficiency(DynamicViscosity, Pressure,
					 Temperature, ParticleDiameter,
					 ParticleDensity, AirDensity,
					 FrictionVelocity, CollectorRadius,
					 Alpha, Beta);
  //cout<<"Eim = "<<Eim<<endl;
  // Interception
  float Ein = ComputeInterceptionEfficiency(ParticleDiameter, CollectorRadius);
  //cout<<"Ein = "<<Ein<<endl;
  // Brownian
  float Eb = ComputeBrownianEfficiency(DynamicViscosity, Pressure, Temperature,
				       ParticleDiameter, AirDensity, Gamma);
  //cout<<"Eb = "<<Eb<<endl;
  // Sticking correction factor
  float Rr = StickingCorrectionFactor(DynamicViscosity, Pressure, Temperature,
				      ParticleDiameter, ParticleDensity,
				      AirDensity, FrictionVelocity,
				      CollectorRadius);
  //cout<<"Rr = "<<Rr<<endl;
  // Surface resistance
	float Rb;
	if(ParticleDiameter > 5e-6)
	  {
	    Rb = 1. / ( Epsilon * FrictionVelocity * ( Eim + Ein + Eb ) * Rr );
	    //cout<<"ParticleDiameter > 5e-6: Rb = "<<Rb<<endl;
	  }
	else
	  {
	    //cout<<"ParticleDiameter < 5e-6: Rb = "<<Rb<<endl;
	    Rb = 1. / ( Epsilon * FrictionVelocity * ( Eim + Ein + Eb ) );
	  }
  return Rb;
}
/*
//     SmallRadius: Characteristic radius of small receptors. ([m])
//     LargeRadius: Characteristic radius of large receptors. ([m])

float ComputeSurfaceResistance(float DynamicViscosity, float Pressure,
			       float Temperature, float ParticleDiameter,
			       float ParticleDensity, float AirDensity,
			       float FrictionVelocity, float LargeRadius,
			       float SmallRadius, float Alpha, float Beta,
			       float Gamma)
{
  // The resistance model for dry deposition postulates that adjacent
  // to the surface exists a quasi-laminar layer.

  // To represent this resistance for particles, an expression developped
  // by Zhang et al. (2001) has been used.
  
  /////////// Collection efficiency coefficient ///////////
  // Impaction
  float Eim = ComputeImpactionEfficiency(DynamicViscosity, Pressure,
					 Temperature, ParticleDiameter,
					 ParticleDensity, AirDensity,
					 FrictionVelocity, LargeRadius,
					 Alpha, Beta);
  // Interception
  float Ein = ComputeInterceptionEfficiency(ParticleDiameter, LargeRadius,
					    SmallRadius);
  // Brownian
  float Eb = ComputeBrownianEfficiency(DynamicViscosity, Pressure, Temperature,
				       ParticleDiameter, AirDensity, Gamma);

  // Sticking correction factor
  float Rr = StickingCorrectionFactor(DynamicViscosity, Pressure, Temperature,
				      ParticleDiameter, ParticleDensity,
				      AirDensity, FrictionVelocity,
				      LargeRadius);

  // Surface resistance
  float Rb = 1. / ( Epsilon * FrictionVelocity * ( Eim + Ein + Eb ) * Rr );

  return Rb;
}
*/

// Renvoi le coefficient de diffusion brownienne calculé via la relation
// de Stokes-Einstein, corrigé par le facteur de Cunningham (1910).
/*
  \param DynamicViscosity Viscosité dyanmique du milieu (air généralement)
    dans lequel se dépose la particule.
  \param Pressure Pression atmosphérique.
  \param Temperature Température du milieu (air).
  \param ParticleDiameter Diamètre de la particule qui se dépose.

  \return Le coefficient de diffusion brownienne pour la particule.
*/
float ComputeBrownianCoeff(float DynamicViscosity, float Pressure,
			   float Temperature, float ParticleDiameter)
{
  //Compute brownian diffusion coefficient

  float CC = ComputeCunninghamFactor(DynamicViscosity, Pressure,
				     Temperature, ParticleDiameter);

  //Einstein equation of brownian diffusion (1905) + Cunningham correction factor
  float D = Kb * Temperature * CC /
    ( 3. * PI * DynamicViscosity * ParticleDiameter );

  return D;
}


// Renvoi le nombre de Schmidt.
/*
  \param DynamicViscosity Viscosité dyanmique du milieu (air généralement)
    dans lequel se dépose la particule.
  \param Pressure Pression atmosphérique.
  \param Temperature Température du milieu (air).
  \param ParticleDiameter Diamètre de la particule qui se dépose.
  \param AirDensity Densité du milieu (air).

  \return Le nombre de Schmidt.
*/
float ComputeSchmidtNumber(float DynamicViscosity, float Pressure,
			   float Temperature, float ParticleDiameter,
			   float AirDensity)
{
//Compute Schmidt Number

  float KinematicViscosity = DynamicViscosity / AirDensity;
  float D = ComputeBrownianCoeff(DynamicViscosity, Pressure,
				 Temperature, ParticleDiameter);

  float Sc = KinematicViscosity / D;

  return Sc;
}


// Renvoi l'effacité de dépôt par diffusion brownienne.
/*
  \param DynamicViscosity Viscosité dyanmique du milieu (air généralement)
    dans lequel se dépose la particule.
  \param Pressure Pression atmosphérique.
  \param Temperature Température du milieu (air).
  \param ParticleDiameter Diamètre de la particule qui se dépose.
  \param AirDensity Densité du milieu (air).
  \param Gamma Constante du modèle de Zhang (dépend du LUC et de la saison SC).

  \return L'effacité de dépôt par diffusion brownienne.
*/
float ComputeBrownianEfficiency(float DynamicViscosity, float Pressure,
				float Temperature, float ParticleDiameter,
				float AirDensity, float Gamma)
{
  //Compute collection efficiency from Brownian diffusion

  float Sc = ComputeSchmidtNumber(DynamicViscosity, Pressure, Temperature,
				  ParticleDiameter, AirDensity);

  float Eb = pow(Sc , -Gamma);
//  float Eb = (1./3.) * pow(Sc , -Gamma);
//  float Eb = pow(Sc , -2./3.);

  return Eb;
}


// Renvoi l'effacité de dépôt par interception.
/*
  \param ParticleDiameter Diamètre de la particule qui se dépose.
  \param CollectorRadius Rayon de collection des éléments du sol
    (valeur tabulée dans le modèle de Zhang).

  \return L'effacité de dépôt par interception.
*/
float ComputeInterceptionEfficiency(float ParticleDiameter,
				    float CollectorRadius)
{
//Collection efficiency from interception
  float Ein, Temp;
  if( CollectorRadius == 0. )
    Ein = 0.;
  else
    {
      Temp = ParticleDiameter / CollectorRadius;
      Ein = 0.5 * pow(Temp, 2.);
    }

  return Ein;
}

// Renvoi l'effacité de dépôt par interception.
/*
  \param ParticleDiameter Diamètre de la particule qui se dépose.
  \param LargeRadius Rayon de collection des éléments du sol
    (valeur tabulée dans le modèle de Zhang).
  \param SmallRadius

  \return L'effacité de dépôt par interception.
*/
/*
float ComputeInterceptionEfficiency(float ParticleDiameter,
				    float LargeRadius,
				    float SmallRadius)

{
//Collection efficiency from interception
  float Ein;

  if( LargeRadius == 0 )
    Ein = 0.;
  else
    {
      Ein = 0.99 * ParticleDiameter / ( LargeRadius + ParticleDiameter )
	+ 0.01 * ParticleDiameter / ( SmallRadius + ParticleDiameter );
    }

  return Ein;
}
*/

// Renvoi le nombre de Stokes (cas d'une surface lisse).
/*
  \param DynamicViscosity Viscosité dyanmique du milieu (air généralement)
    dans lequel se dépose la particule.
  \param Pressure Pression atmosphérique.
  \param Temperature Température du milieu (air).
  \param ParticleDiameter Diamètre de la particule qui se dépose.
  \param ParticleDensity Densité de la particule qui se dépose.
  \param AirDensity Densité du milieu (air).
  \param FrictionVelocity Vitesse de friction.

  \return Le nombre de Stokes (cas d'une surface lisse).
*/

float ComputeStokesNumberSmooth(float DynamicViscosity, float Pressure,
				float Temperature, float ParticleDiameter,
				float ParticleDensity, float AirDensity,
				float FrictionVelocity)
{
  //Compute the Stokes number on smooth surfaces (water, ice...)

  float vs = ComputeSedimentationVelocity(DynamicViscosity, Pressure,
					  Temperature, ParticleDiameter,
					  ParticleDensity);
  float KinematicViscosity = DynamicViscosity / AirDensity;

  float St = vs * pow(FrictionVelocity, 2.) / ( g * KinematicViscosity );

  return St;
}

// Renvoi le nombre de Stokes (cas d'une surface rugueuse).
/*
  \param DynamicViscosity Viscosité dyanmique du milieu (air généralement)
    dans lequel se dépose la particule.
  \param Pressure Pression atmosphérique.
  \param Temperature Température du milieu (air).
  \param ParticleDiameter Diamètre de la particule qui se dépose.
  \param ParticleDensity Densité de la particule qui se dépose.
  \param AirDensity Densité du milieu (air).
  \param FrictionVelocity Vitesse de friction.
  \param CollectorRadius Rayon de collection des éléments du sol
    (valeur tabulée dans le modèle de Zhang).

  \return Le nombre de Stokes (cas d'une surface rugueuse).
*/
float ComputeStokesNumberRough(float DynamicViscosity, float Pressure,
			       float Temperature, float ParticleDiameter,
			       float ParticleDensity, float AirDensity,
			       float FrictionVelocity, float CollectorRadius)
{
  //Compute the Stokes number on rough surfaces (wood, urban...)

  float vs = ComputeSedimentationVelocity(DynamicViscosity, Pressure,
					  Temperature, ParticleDiameter,
					  ParticleDensity);
  float St = vs * FrictionVelocity / ( g * CollectorRadius );

  return St;
}


// Renvoi l'effacité de dépôt par impaction.
/*
  \param DynamicViscosity Viscosité dyanmique du milieu (air généralement)
    dans lequel se dépose la particule.
  \param Pressure Pression atmosphérique.
  \param Temperature Température du milieu (air).
  \param ParticleDiameter Diamètre de la particule qui se dépose.
  \param ParticleDensity Densité de la particule qui se dépose.
  \param AirDensity Densité du milieu (air).
  \param FrictionVelocity Vitesse de friction.
  \param CollectorRadius Rayon de collection des éléments du sol
    (valeur tabulée dans le modèle de Zhang).
  \param Alpha Constante du modèle de Zhang (dépend du LUC et de la saison SC).
  \param Beta Constante du modèle de Zhang (dépend du LUC et de la saison SC).

  \return L'effacité de dépôt par impaction.
*/
float ComputeImpactionEfficiency(float DynamicViscosity, float Pressure,
				 float Temperature, float ParticleDiameter,
				 float ParticleDensity, float AirDensity,
				 float FrictionVelocity, float CollectorRadius,
				 float Alpha, float Beta)
{
  //Collection efficiency from impaction

  float St ;
// float Temp;	
  if( CollectorRadius == 0 )
  {
    St = ComputeStokesNumberSmooth(DynamicViscosity, Pressure, Temperature,
				   ParticleDiameter, ParticleDensity,
				   AirDensity, FrictionVelocity);
//    Temp = pow(10.,(-3./St));
  }
  else
  {
    St = ComputeStokesNumberRough(DynamicViscosity, Pressure, Temperature,
				  ParticleDiameter, ParticleDensity, AirDensity,
				  FrictionVelocity, CollectorRadius);
//    Temp = St / (1. + pow(St,2.));
  }
  float Temp = St / ( St + Alpha );
  float Eim = pow(Temp, Beta);
//	float Eim = Temp;
  
  return Eim;
}

// Renvoi le libre parcours moyen dans l'air.
/*
  \param DynamicViscosity Viscosité dyanmique du milieu (air généralement)
    dans lequel se dépose la particule.
  \param Pressure Pression atmosphérique.
  \param Temperature Température du milieu (air).

  \return Le libre parcours moyen dans l'air.
*/

float ComputeLambda(float DynamicViscosity, float Pressure, float Temperature)
{
  // Compute lamba, which is the air mean free path

  // Air mean free path [Sportisse (2006)]
//  float Lambda = 2. * DynamicViscosity /
//    ( Pressure * pow( 8. / ( PI * R * Temperature), 0.5 ) );
  float Lambda = (2. * DynamicViscosity / Pressure) * 
                 pow( 8. / ( PI * Rair * Temperature), -0.5 );

  return Lambda;
}

// Renvoi le facteur de correction de Cunningham (1910).
/*
  \param DynamicViscosity Viscosité dyanmique du milieu (air généralement)
    dans lequel se dépose la particule.
  \param Pressure Pression atmosphérique.
  \param Temperature Température du milieu (air).
  \param ParticleDiameter Diamètre de la particule qui se dépose.

  \return Le facteur de correction de Cunningham (1910).
*/

float ComputeCunninghamFactor(float DynamicViscosity, float Pressure,
			      float Temperature, float ParticleDiameter)
{
  //Compute the Cunningham correction factor, which is a correction
  // to the drag coefficient

  float Lambda = ComputeLambda(DynamicViscosity, Pressure, Temperature);

//  float CC = 1. + 2. * Lambda / ParticleDiameter *
//    ( A1 + A2 * exp( -A3 * ParticleDiameter / ( 2. * Lambda ) ) );
  float CC = 1. + (2. * Lambda / ParticleDiameter) *
    ( A1 + A2 * exp( -A3 * ParticleDiameter / Lambda ) );

  return CC;
}

// Renvoi le coefficient de rebond.
/*
  \param DynamicViscosity Viscosité dyanmique du milieu (air généralement)
    dans lequel se dépose la particule.
  \param Pressure Pression atmosphérique.
  \param Temperature Température du milieu (air).
  \param ParticleDiameter Diamètre de la particule qui se dépose.
  \param ParticleDensity Densité de la particule qui se dépose.
  \param AirDensity Densité du milieu (air).
  \param FrictionVelocity Vitesse de friction.
  \param CollectorRadius Rayon de collection des éléments du sol
    (valeur tabulée dans le modèle de Zhang).

  \return Le coefficient de rebond.
*/

float StickingCorrectionFactor(float DynamicViscosity, float Pressure,
			       float Temperature, float ParticleDiameter,
			       float ParticleDensity, float AirDensity,
			       float FrictionVelocity, float CollectorRadius)
{
  // This factor represents the fraction of particles, once in contact,
  // that stick to surface.

  //This expression was suggested by Slinn (1982)
  float St;
  if( CollectorRadius == 0 )
    St = ComputeStokesNumberSmooth(DynamicViscosity, Pressure, Temperature,
				   ParticleDiameter, ParticleDensity,
				   AirDensity, FrictionVelocity);
  else
    St = ComputeStokesNumberRough(DynamicViscosity, Pressure, Temperature,
				  ParticleDiameter, ParticleDensity,
				  AirDensity, FrictionVelocity, CollectorRadius);

//  float Rr = exp (- sqrt ( St ) );
    float Rr = exp ( -2.0 * sqrt ( St ) );

  return Rr;
}


// Renvoi la vitesse de sédimentation.
/*
  \param DynamicViscosity Viscosité dyanmique du milieu (air généralement)
    dans lequel se dépose la particule.
  \param Pressure Pression atmosphérique.
  \param Temperature Température du milieu (air).
  \param ParticleDiameter Diamètre de la particule qui se dépose.
  \param ParticleDensity Densité de la particule qui se dépose.

  \return La vitesse de sédimentation.
*/

float ComputeSedimentationVelocity(float DynamicViscosity, float Pressure,
				   float Temperature, float ParticleDiameter,
				   float ParticleDensity)
{
  // Compute the sedimentation velocity for particules which the radius is
  // smaller than 20micrometers [Seinfield and Pandis (1998), p.909]

  float CC = ComputeCunninghamFactor(DynamicViscosity, Pressure,
				     Temperature, ParticleDiameter);

  float vs = pow(ParticleDiameter, 2.) * g * CC * ParticleDensity /
    ( 18. * DynamicViscosity );
  
  return vs;
}

