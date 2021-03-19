#ifndef LIBRS_HXX
#define LIBRS_HXX

#include<cmath>

// Pi
const float PI = 3.14159265358979323846f;
// Constante de Boltzmann (in J/K or m^2.kg/s^2/K)
const float Kb = 1.3806503e-23;
// Constante issue du modèle de Zhang
//const float Epsilon = 3.0;
const float Epsilon = 2.0;
// Constante gravitationelle
const float g = 9.81;
// Constante des gaz parfaits (J/K/mol)
const float R = 8.31447;
// Constante molaire de lair (J/K/kg)
const float Rair = 287.05;
// Masse volumique de l'air (kg/m3)
const float AirDensity = 1.2;
// Cunningham factor depends on 3 coefficient [Knudsen and Weber (1911)],
//  which are experimentally determined [Davies (1945)] :
const float A1 = 1.257;
const float A2 = 0.400;
const float A3 = 0.55;

float ComputeSurfaceResistance(float DynamicViscosity, float Pressure,
			       float Temperature, float ParticleDiameter,
			       float ParticleDensity, float AirDensity,
			       float FrictionVelocity, float CollectorRadius,
			       float Alpha, float Beta, float Gamma);

float ComputeSurfaceResistance(float DynamicViscosity, float Pressure,
			       float Temperature, float ParticleDiameter,
			       float ParticleDensity, float AirDensity,
			       float FrictionVelocity, float LargeRadius,
			       float SmallRadius, float Alpha, float Beta,
			       float Gamma);

float ComputeBrownianCoeff(float DynamicViscosity, float Pressure,
			   float Temperature, float ParticleDiameter);

float ComputeSchmidtNumber(float DynamicViscosity, float Pressure,
			   float Temperature, float ParticleDiameter,
			   float AirDensity);

float ComputeBrownianEfficiency(float DynamicViscosity, float Pressure,
				float Temperature, float ParticleDiameter,
				float AirDensity, float Gamma);

float ComputeInterceptionEfficiency(float ParticleDiameter,
				    float CollectorRadius);

float ComputeInterceptionEfficiency(float ParticleDiameter,
				    float LargeRadius,
				    float SmallRadius);

float ComputeStokesNumberSmooth(float DynamicViscosity, float Pressure,
				float Temperature, float ParticleDiameter,
				float ParticleDensity, float AirDensity,
				float FrictionVelocity);

float ComputeStokesNumberRough(float DynamicViscosity, float Pressure,
			       float Temperature, float ParticleDiameter,
			       float ParticleDensity, float AirDensity,
			       float FrictionVelocity, float CollectorRadius);

float ComputeImpactionEfficiency(float DynamicViscosity, float Pressure,
				 float Temperature, float ParticleDiameter,
				 float ParticleDensity, float AirDensity,
				 float FrictionVelocity, float CollectorRadius,
				 float Alpha, float Beta);

float ComputeLambda(float DynamicViscosity, float Pressure, float Temperature);

float ComputeCunninghamFactor(float DynamicViscosity, float Pressure,
			      float Temperature, float ParticleDiameter);

float StickingCorrectionFactor(float DynamicViscosity, float Pressure,
			       float Temperature, float ParticleDiameter,
			       float ParticleDensity, float AirDensity,
			       float FrictionVelocity, float CollectorRadius);

float ComputeSedimentationVelocity(float DynamicViscosity, float Pressure,
				   float Temperature, float ParticleDiameter,
				   float ParticleDensity);

#endif
