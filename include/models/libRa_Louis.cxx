#include"libRa_Louis.hxx"

/******************************************************/
/*  Compute the surface potential temperature         */
/*                                                    */
/* param:  SurfaceTemperature   - surface temperature */
/*         SurfacePressure      - surface pressure    */
/*         Rair                 - perfect gas         */
/*                                 constant for air   */
/*         Cp                   - specific heat to    */
/*                                 constant pressure  */
/*                                                    */
/* return: SurfacePotentialTemp - surface potential   */
/*                                 temperature        */
/******************************************************/
float ComputeSurfacePotentialTemp(float SurfaceTemperature,
                                  float SurfacePressure,
                                  float Rair, float Cp)
{
  float SurfPotentialTemp =
    SurfaceTemperature * pow(
                             SurfacePressure / pow(10, 5),
                             -Rair / Cp
                             );
  return SurfPotentialTemp;
}


/*****************************************************/
/*  Compute the potential temperature                */
/*                                                   */
/* param:  Temperature    - temperature              */
/*         Pressure       - pressure                 */
/*         Rair           - perfect gas constant for */
/*                             air                   */
/*         Cp             - specific heat to         */
/*                             constant pressure     */
/*                                                   */
/* return: PotentialTemp  - potential temperature    */
/*****************************************************/
float ComputePotentialTemp(float Temperature, float Pressure,
                           float Rair, float Cp)
{
  float PotentialTemp = Temperature * pow(
                                          Pressure / pow(10, 5),
                                          -Rair / Cp
                                          );
  return PotentialTemp;
}


/*****************************************************/
/*  Compute the Richardson number                    */
/*                                                   */
/* param: SurfacePotentialTemp - surface potential   */
/*                                 temperature       */
/*        PotentialTemp        - potential           */
/*                                 temperature       */
/*        g                    - gravity             */
/*                                 acceleration      */
/*        zref                 - height of the first */
/*                                 vertical node     */
/*        FirstLevelWindModule - wind module at the  */
/*                               first vertical node */
/*                                                   */
/* return: Richardson          - richardson number   */
/*****************************************************/
float ComputeRichardson(float SurfacePotentialTemp, float PotentialTemp,
                        float g, float zref, float FirstLevelWindModule)
{
  float Ri = 2 * g * zref * ( PotentialTemp - SurfacePotentialTemp )
    / ( PotentialTemp + SurfacePotentialTemp )
    / pow( FirstLevelWindModule, 2 );

  return Ri;
}


/*****************************************************/
/*  Compute the thermal roughness length             */
/*                                                   */
/* param:  z0          - dynamical roughness length  */
/*                                                   */
/* return: z0t         - thermal roughness length    */
/*****************************************************/
float ComputeZ0t(float z0)
{
  float z0t = z0 * exp(-2.);

  return z0t;
}

/*****************************************************/
/*  Compute the dynamic stability                    */
/*                                                   */
/* param:  zref        - height of the first         */
/*                          vertical node            */
/*         z0          - dynamical roughness length  */
/*         z0t         - thermal roughness length    */
/*         Richardson  - richardson number           */
/*                                                   */
/* return: Stab  - thermic stability                 */
/*****************************************************/
float ComputeDynamicalStability(float zref, float z0,
				float z0t, float Richardson)
{
  // Von Karman constant
  float Ka = 0.41;

  // Local variables
  float alu = Ka / log( (zref + z0) / z0 );
  float ZsZ0t = (zref + z0t) / z0t;
  float alt = Ka / log( ZsZ0t );

  // Constants of Louis' formulas
  float b = 5., c = 5., d = 5.;

  float Stab;

  if( Richardson < 0 )
    {
      Stab = 1. + 2. * b * fabs(Richardson) /
	( 1. + 3. * b * c * alu * alt *
	  sqrt( fabs(Richardson) ) * sqrt( 1. - 1. / ZsZ0t )
	  * pow( pow(ZsZ0t, float(1./3.) ) - 1., float(3./2.) ) );
    }
  else
    {
      Stab = 1. / ( 1. + 2. * b * fabs(Richardson)
		    / sqrt( 1. + d * fabs(Richardson) ) );
    }

  return Stab;
}

/*****************************************************/
/*  Compute the dynamic stability                    */
/*                                                   */
/* param: zref                 - height of the first */
/*                                 vertical node     */
/*        z0                   - dynamical roughness */
/*                                 length            */
/*        SurfacePotentialTemp - surface potential   */
/*                                 temperature       */
/*        PotentialTemp        - potential           */
/*                                 temperature       */
/*        g                    - gravity             */
/*                                 acceleration      */
/*        FirstLevelWindModule - wind module at the  */
/*                               first vertical node */
/*                                                   */
/* return: Stab                - thermic stability   */
/*****************************************************/
float ComputeDynamicalStability(float zref, float z0, float SurfacePotentialTemp,
				float PotentialTemp, float g,
				float FirstLevelWindModule)
{
  // Von Karman constant
  float Ka = 0.41;

  // Thermal roughness length
  float z0t = ComputeZ0t(z0);

  // Local variables
  float alu = Ka / log( (zref + z0) / z0 );
  float ZsZ0t = (zref + z0t) / z0t;
  float alt = Ka / log( ZsZ0t );

  // Constants of Louis' formulas
  float b = 5., c = 5., d = 5.;

  float Richardson = ComputeRichardson(SurfacePotentialTemp, PotentialTemp,
                                       g, zref, FirstLevelWindModule);

  float Stab;

  if( Richardson < 0 )
    {
      Stab = 1. + 2. * b * fabs(Richardson)
	/ ( 1. + 3. * b * c * alu * alt * sqrt( fabs(Richardson) )
	    * sqrt( 1. - 1. / ZsZ0t ) * pow( pow(ZsZ0t, 1./3.) - 1., 1.5 ) );
    }
  else
    {
      Stab = 1. / ( 1. + 2. * b * fabs(Richardson) /  sqrt( 1. + d * fabs(Richardson) ) );
    }

  return Stab;
}

/*****************************************************/
/*  Compute the thermic stability                    */
/*                                                   */
/* param:  zref        - height of the first         */
/*                          vertical node            */
/*         z0          - dynamical roughness length  */
/*         z0t         - thermal roughness length    */
/*         Richardson  - richardson number           */
/*                                                   */
/* return: Stab  - thermic stability                 */
/*****************************************************/
float ComputeThermalStability(float zref, float z0,
                              float z0t, float Richardson)
{
  // Von Karman constant
  float Ka = 0.41;

  // Local variables
  float alu = Ka / log( (zref + z0) / z0 );
  float ZsZ0t = (zref + z0t) / z0t;
  float alt = Ka / log( ZsZ0t );

  // Constants of Louis' formulas
  float b = 5., c = 5., d = 5.;

  float Stab;

  if( Richardson < 0 )
    {
      Stab = 1. + 3. * b * fabs(Richardson) /
	( 1. + 3. * b * c * alu * alt *
	  sqrt( fabs(Richardson) ) * sqrt( 1. - 1. / ZsZ0t )
	  * pow( pow(ZsZ0t, float(1./3.) ) - 1., float(3./2.) ) );
    }
  else
    {
      Stab = 1. / ( 1. + 3. * b * fabs(Richardson)
		    * sqrt( 1. + d * fabs(Richardson) ) );
    }

  return Stab;
}

/*****************************************************/
/*  Compute the thermic stability                    */
/*                                                   */
/* param: zref                 - height of the first */
/*                                 vertical node     */
/*        z0                   - dynamical roughness */
/*                                 length            */
/*        SurfacePotentialTemp - surface potential   */
/*                                 temperature       */
/*        PotentialTemp        - potential           */
/*                                 temperature       */
/*        g                    - gravity             */
/*                                 acceleration      */
/*        FirstLevelWindModule - wind module at the  */
/*                               first vertical node */
/*                                                   */
/* return: Stab                - thermic stability   */
/*****************************************************/
float ComputeThermalStability(float zref, float z0, float SurfacePotentialTemp,
                              float PotentialTemp, float g,
                              float FirstLevelWindModule)
{
  // Von Karman constant
  float Ka = 0.41;

  // Thermal roughness length
  float z0t = ComputeZ0t(z0);

  // Local variables
  float alu = Ka / log( (zref + z0) / z0 );
  float ZsZ0t = (zref + z0t) / z0t;
  float alt = Ka / log( ZsZ0t );

  // Constants of Louis' formulas
  float b = 5., c = 5., d = 5.;

  float Richardson = ComputeRichardson(SurfacePotentialTemp, PotentialTemp,
                                       g, zref, FirstLevelWindModule);

  float Stab;

  if( Richardson < 0 )
    {
      Stab = 1. + 3. * b * fabs(Richardson) / ( 1. + 3. * b * c * alu * alt * sqrt( fabs(Richardson) ) * sqrt( 1. - 1. / ZsZ0t ) * pow( pow(ZsZ0t, 1./3.) - 1., 1.5 ) );
    }
  else
    {
      Stab = 1. / ( 1. + 3. * b * fabs(Richardson) *  sqrt( 1. + d * fabs(Richardson) ) );
    }

  return Stab;
}

/*****************************************************/
/*  Compute the friction velocity                    */
/*         (from Louis' formula)                     */
/*                                                   */
/* param:  zref        - height of the first         */
/*                          vertical node            */
/*         z0          - dynamical roughness length  */
/*         Richardson  - richardson number           */
/*         FirstLevelWindModule                      */
/*                     - wind module at the          */
/*                          first vertical node      */
/*                                                   */
/* return: ustar       - friction velocity           */
/*****************************************************/
float ComputeUstar(float zref, float z0, float Richardson,
		   float FirstLevelWindModule)
{
  // Von Karman constant
  float Ka = 0.41;

  // Thermic roughness height
  float z0t = ComputeZ0t(z0);

  // Local variables
  float alu = Ka / log( (zref + z0) / z0 );

  // Dynamic stability
  float Stab = ComputeDynamicalStability(zref, z0, z0t, Richardson);

  float ustar = alu * sqrt(Stab) * FirstLevelWindModule;

  return ustar;
}

float ComputeUstar(float zref, float z0, float FirstLevelWindModule)
{
  // Von Karman constant
  float Ka = 0.41;

  float ustar = Ka / log( (zref + z0) / z0 ) * FirstLevelWindModule;

  return ustar;
}

/*****************************************************/
/*  Compute the friction velocity                    */
/*         (from Louis' formula)                     */
/*                                                   */
/* param: zref                 - height of the first */
/*                                 vertical node     */
/*        z0                   - dynamical roughness */
/*                                 length            */
/*        SurfacePotentialTemp - surface potential   */
/*                                 temperature       */
/*        PotentialTemp        - potential           */
/*                                 temperature       */
/*        g                    - gravity             */
/*                                 acceleration      */
/*        FirstLevelWindModule - wind module at the  */
/*                               first vertical node */
/*                                                   */
/* return: ustar               - friction velocity   */
/*****************************************************/
float ComputeUstar(float zref, float z0, float SurfacePotentialTemp,
		   float PotentialTemp, float g,
		   float FirstLevelWindModule)
{
  // Von Karman constant
  float Ka = 0.41;

  // Thermic roughness height
  float z0t = ComputeZ0t(z0);

  // Richardson number
  float Richardson = ComputeRichardson(SurfacePotentialTemp, PotentialTemp,
                                       g, zref, FirstLevelWindModule);

  // Local variables
  float alu = Ka / log( (zref + z0) / z0 );

  // Dynamic stability
  float Stab = ComputeDynamicalStability(zref, z0, z0t, Richardson);

  float ustar = alu * sqrt(Stab) * FirstLevelWindModule;

  return ustar;
}

/*****************************************************/
/*  Compute the aerodynamic resistance               */
/*         (from Louis' formula)                     */
/*                                                   */
/* param:  zref        - height of the first         */
/*                          vertical node            */
/*         z0          - dynamical roughness length  */
/*         Richardson  - richardson number           */
/*         FirstLevelWindModule                      */
/*                     - wind module at the          */
/*                          first vertical node      */
/*                                                   */
/* return: Ra          - aerodynamic resistance      */
/*****************************************************/
float ComputeRaLouis(float zref, float z0, float Richardson,
                     float FirstLevelWindModule)
{
  // Von Karman constant
  float Ka = 0.41;

  // Thermic roughness height
  float z0t = ComputeZ0t(z0);

  // Local variables
  float alu = Ka / log( (zref + z0) / z0 );
  float alt = Ka / log( (zref + z0t) / z0t );

  // Thermic stability
  float Stab = ComputeThermalStability(zref, z0, z0t, Richardson);

  float Ra = 1. / ( alu * alt * Stab * FirstLevelWindModule );

  return Ra;
}

/*****************************************************/
/*  Compute the aerodynamic resistance               */
/*         (from Louis' formula)                     */
/*                                                   */
/* param: zref                 - height of the first */
/*                                 vertical node     */
/*        z0                   - dynamical roughness */
/*                                 length            */
/*        SurfacePotentialTemp - surface potential   */
/*                                 temperature       */
/*        PotentialTemp        - potential           */
/*                                 temperature       */
/*        g                    - gravity             */
/*                                 acceleration      */
/*        FirstLevelWindModule - wind module at the  */
/*                               first vertical node */
/*                                                   */
/* return: Ra                  - aerodynamic         */
/*                                 resistance        */
/*****************************************************/
float ComputeRaLouis(float zref, float z0, float SurfacePotentialTemp,
                     float PotentialTemp, float g,
                     float FirstLevelWindModule)
{
  // Von Karman constant
  float Ka = 0.41;

  // Thermic roughness height
  float z0t = ComputeZ0t(z0);

  // Richardson number
  float Richardson = ComputeRichardson(SurfacePotentialTemp, PotentialTemp,
                                       g, zref, FirstLevelWindModule);

  // Local variables
  float alu = Ka / log( (zref + z0) / z0 );
  float alt = Ka / log( (zref + z0t) / z0t );

  // Thermic stability
  float Stab = ComputeThermalStability(zref, z0, z0t, Richardson);

  float Ra = 1. / ( alu * alt * Stab * FirstLevelWindModule );

  return Ra;
}
