#ifndef LIBRA_LOUIS_HXX
#define LIBRA_LOUIS_HXX

#include<cmath>

float ComputeSurfacePotentialTemp(float SurfaceTemperature,
                                  float SurfacePressure,
                                  float Rair, float Cp);
float ComputePotentialTemp(float Temperature, float Pressure,
                           float Rair, float Cp);
float ComputeRichardson(float SurfacePotentialTemp, float PotentialTemp,
                        float g, float zref, float FirstLevelWindModule);
float ComputeZ0t(float z0);
float ComputeDynamicalStability(float zref, float z0,
				float z0t, float Richardson);
float ComputeDynamicalStability(float zref, float z0,
				float SurfacePotentialTemp,
				float PotentialTemp, float g,
				float FirstLevelWindModule);
float ComputeUstar(float zref, float z0, float Richardson,
		   float FirstLevelWindModule);
float ComputeUstar(float zref, float z0, float FirstLevelWindModule);
float ComputeUstar(float zref, float z0, float SurfacePotentialTemp,
		   float PotentialTemp, float g,
		   float FirstLevelWindModule);
float ComputeRaLouis(float zref, float z0, float Richardson,
                     float FirstLevelWindModule);
float ComputeRaLouis(float zref, float z0, float SurfacePotentialTemp,
                     float PotentialTemp, float g,
                     float FirstLevelWindModule);

#endif
