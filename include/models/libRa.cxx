#include "libRa.hxx"

/*********************************************************************************************************/
/*                                   Defining characteristic lengths                                     */
/*********************************************************************************************************/

/*****************************************************/
/*        Compute the displacement height            */
/*               (Macdonald, 1998)                   */
/*                                                   */
/* param:  A        - empirical constant             */
/*         lambda_p - plan area density              */
/*         h        - buildings' mean height         */
/*                                                   */
/* return: d        - displacement height            */
/*****************************************************/
float ComputeD(float A, float lambda_p, float h)
{
  float d = h * ( 1. + pow(A, -lambda_p) * ( lambda_p - 1. ));
  return d;
}

/*****************************************************/
/*          Compute the roughness height             */
/*               (Macdonald, 1998)                   */
/*                                                   */
/* param:  lambda_p - plan area density              */
/*         h        - buildings' mean height         */
/*         d        - displacement height            */
/*         Cd       - drag coefficient (ESDU, 1980)  */
/*                                                   */
/* return: z0       - roughness height               */
/*****************************************************/
/*
float ComputeZ0(float lambda_p, float h, float d, float Cd)
{
  float kappa = 0.41;

  float z0 = (1. - d/h ) * exp( -pow(0.5 * Cd / (kappa * kappa) * ( 1 - d/h ) * lambda_p, -0.5) );

  return z0;
}

float computeZ0(float lambda_p, float h)
{
  float z0 = h * exp( -0.52 * pow(lambda_p, -0.5) );

  return z0;
}
*/
// Macdonald1 = eq. 22 (Macdonald, 1998) - 
float ComputeZ0(float lambda_f, float h, float d, float Cd)
{
  float kappa = 0.40;

  float z0 = h * (1. - d/h ) * exp( - pow(0.5 * Cd * ( 1 - d/h ) * lambda_f / pow(kappa,2), -0.5) );

  return z0;
}

// Macdonald2 = eq 18. (Macdonald, 1998) 
float computeZ0(float lambda_f, float h, float Cd)
{
  float kappa = 0.40;

  float z0 = h * exp( - pow(0.5 * Cd * lambda_f / pow(kappa,2), -0.5) );

  return z0;
}

/*****************************************************/
/*    Compute the canyon's caracteristic length      */
/*           ( harmonic mean of the middle           */
/*             of the length and width of            */
/*             the canyon's region )                 */
/*            1/z_can = 1/(w/2) + 1/(l/2) 					 */
/*                                                   */
/* param:  w        - canyon's region mean width     */
/*         l        - canyon's region mean height    */
/*                                                   */
/* return: z_can    - canyon's caracteristic length  */
/*****************************************************/

float ComputeZcan(float w, float l)
{
  float z_can = w * l / ( 2. * (w + l));
  return z_can;
}

/*****************************************************/
/*  Compute the average mixing length in the canopy  */
/*                                                   */
/* param:  h        - buildings' mean height         */
/*         d        - displacement height            */
/*                                                   */
/* return: lc       - average mixing length in the   */
/*                     canopy                        */
/*****************************************************/
float ComputeLc(float h, float d)
{
  // Von Karman constant
  float Ka = 0.41;

  float lc = Ka * h * ( h - d ) / d;
  return lc;
}

float ComputeLm(float z, float lc)
{
  // Von Karman constant
  float Ka = 0.41;

  float lm = lc * Ka * z / ( lc + Ka * z );
  return lm;
}

float ComputeZlim(float lc, float alpha)
{
  // Von Karman constant
  float Ka = 0.41;

  float zlim = alpha * lc / ( Ka * (1 - alpha) );
  return zlim;
}

/*********************************************************************************************************/
/*                                           Dynamic properties                                          */
/*********************************************************************************************************/

/*****************************************************/
/*  Compute the friction velocity above urban area   */
/*         (neutral condition is supposed)           */
/*                                                   */
/*  \frac{\partial u}{\partial z} =                  */
/*       \frac{u^*}{\kappa \left( z - d \right)}     */
/*  integrated between FirstLevelWindModule (z=zref) */
/*                     and 0 (z=z0+d)                */
/*                                                   */
/* param:  zref        - height of the first         */
/*                          vertical node            */
/*         z0          - dynamical roughness length  */
/*         d           - displacement height         */
/*         FirstLevelWindModule                      */
/*                     - wind module at the          */
/*                          first vertical node      */
/*                                                   */
/* return: ustar       - friction velocity           */
/*****************************************************/
float ComputeUstar_AboveCanopy(float zref, float z0, float d, float FirstLevelWindModule)
{
  // Von Karman constant
  float Ka = 0.41;

  float Ustar = FirstLevelWindModule * Ka / log( (zref - d) / z0 );

  return Ustar;
}
/*
float ComputeUstar_nonurban(float zref, float z0, float Ri, float uw)
{
  return 0;
}
*/
/*****************************************************/
/*  Compute the friction velocity above urban area   */
/*       (stability is taken into account as         */
/*        it is done for the Louis' formula)         */
/*                                                   */
/* param:  zref        - height of the first         */
/*                          vertical node            */
/*         z0          - dynamical roughness length  */
/*         d           - displacement height         */
/*         Richardson  - richardson number           */
/*         FirstLevelWindModule                      */
/*                     - wind module at the          */
/*                          first vertical node      */
/*                                                   */
/* return: ustar       - friction velocity           */
/*****************************************************/
float ComputeUstar_AboveCanopy(float zref, float z0, float d, float Richardson, float FirstLevelWindModule)
{
  // Von Karman constant
  float Ka = 0.41;

  // Thermic roughness height
  float z0t = ComputeZ0t(z0);

  // Dynamic stability
  float Stab = ComputeDynamicalStability(zref, z0, z0t, Richardson);

  float Ustar = FirstLevelWindModule * Ka * sqrt(Stab) / log( (zref - d) / z0 );

  return Ustar;
}

/*****************************************************/
/*  Compute the friction velocity inside urban area  */
/*         (neutral condition is supposed)           */
/*****************************************************/
/*
float ComputeUstar_InsideCanopy(float lc, float h, float Utop, float z0)
{
  // Von Karman constant
  float Ka = 0.41;

  float Ustar = Utop / ( (h - z0) / lc + log(h/z0) / Ka );

  return Ustar;
}

float ComputeUstar_Surface(float z0, float zlim, float h, float WindProfile_param)
{
  // Von Karman constant
  float Ka = 0.41;

  float Ustar = Ka * exp( WindProfile_param * (zlim/h - 1.) ) / log( zlim / z0 );
  return Ustar;
}
*/

float ComputeUstar_Surface(float z0, float zlim, float h, float WindProfile_param, float Utop)
{
  // Von Karman constant
  float Ka = 0.41;

  float Ustar = Utop * Ka * exp( WindProfile_param * (zlim/h - 1.) ) / log( zlim / z0 );

  return Ustar;
}

/*****************************************************/
/* Compute the wind module at the top of the canyon  */
/*                                                   */
/* param:  zref     - Height of the atmospheric      */
/*                     model's first level           */
/*         h        - Mean buildings' height         */
/*         d        - displacement height            */
/*         Ustar    - Friction velocity              */
/*         FirstLevelWindModule                      */
/*                  - Wind module at the first       */
/*                     vertical node                 */
/*                                                   */
/* return: Utop     - Wind module at the top of the  */
/*                      canyon                       */
/*****************************************************/
/*
float ComputeUtop(float zref, float h, float d, float Ustar,
                  float FirstLevelWindModule)
{
  float kappa = 0.41;
//  float Utop = FirstLevelWindModule +
  float Utop = FirstLevelWindModule -
    Ustar * log( (zref - d) / (h - d) ) / kappa;

  return Utop;
}
*/
float ComputeUtop(float zref, float h, float d, float Ustar,
                  float FirstLevelWindModule, float z0)
{
  float Utop = FirstLevelWindModule *
               log( (h - d) / z0 ) / log( (zref - d) / z0 );

  if (Utop < 0.)
    Utop = 0.;

  return Utop;
}

/*
float ComputeUtop(float zref, float z0, float h, float d, float Ustar)
{
  float kappa = 0.41;
  float Utop = Ustar / kappa * log(h - d) / log(z0);

  return Utop;
}
*/

float ComputeUtop_CanyonIntegration(float w, float h, float w_r)
{
  float coeff = 0;

  // Skimming flow regime (narrow canyon)
  if(w < w_r / 2.)
    {
      coeff = 2. / M_PI;
    }
  // Wake interference flow regime
  else if(w > w_r / 2. && w < w_r)
    {
      coeff = 1. + 3. * ( 2. / M_PI - 1.) * ( h / w - 1. / 3.);
    }
  // Isolated roughness flow regime (large canyon)
  else if(w > w_r)
    {
      coeff = 1;
    }

  return coeff;
}

/*********************************************************************************************************/
/*                                                                                                       */
/*                                   Turbulence scheme: Common path                                      */
/*                                                                                                       */
/*********************************************************************************************************/

/*****************************************************/
/*     Compute the resistance in the log region      */
/*                                                   */
/* param:  z        - computed height                */
/*         zref     - height of the atmospheric      */
/*                     model's first level           */
/*         d        - displacement height            */
/*         Ustar    - friction velocity              */
/*                                                   */
/* return: Ra(z)    - aerodynamical resistance       */
/*                     at height z                   */
/*****************************************************/
float ComputeRa_z(float z, float zref, float d, float Ustar)
{
  // Von Karman constant
  float kappa = 0.41;

  float R = log( (zref - d) / (z - d) ) / (kappa * Ustar);

  return R;
}

/*****************************************************/
/*     Compute the resistance above the canyon       */
/*                                                   */
/* param:  zref     - height of the atmospheric      */
/*                     model's first level           */
/*         h        - buildings' mean height         */
/*         d        - displacement height            */
/*         Ustar    - friction velocity              */
/*                                                   */
/* return: RcanSup  - aerodynamical resistance       */
/*                     above the canyon              */
/*****************************************************/
float ComputeRcanSup(float zref, float h, float d, float Ustar)
{
  // Von Karman constant
  float kappa = 0.41;

  float R = log( (zref - d) / (h - d) ) / (kappa * Ustar);

  return R;
}


/*********************************************************************************************************/
/*                                                                                                       */
/*                              Turbulence scheme: Improved Prandtl's model                              */
/*                                                                                                       */
/*********************************************************************************************************/

/*****************************************************/
/*    Compute the resistance of below part of the    */
/*       canyon, in the recirculation region         */
/*                                                   */
/* param:  d        - displacement height            */
/*         lc       - average mixing length in the   */
/*                     canopy                        */
/*         X        - canyon's caracteristic length  */
/*         w        - canyon's mean width            */
/*         h        - buildings' mean height         */
/*         Ustar    - friction velocity              */
/*         FirstLevelWindModule                      */
/*                  - wind module at the first       */
/*                     vertical node                 */
/*         zref     - height of the atmospheric      */
/*                     model's first level           */
/*                                                   */
/* return: RcanInf  - aerodynamical resistance of    */
/*                     below part of the canyon      */
/*****************************************************/
/*
float ComputeRcanInf_recirculation(float d, float lc, float z_can, float w, float h, float zref,
				   float Ustar, float FirstLevelWindModule)
{
  float kappa = 0.41;
  float Utop = ComputeUtop(d, h, zref, Ustar, FirstLevelWindModule);

  float Coeff = 4. * w / ( pow(lc, 2) * pow(kappa, 2) * pow(h, 2) * Utop);
  float Integrale = Trapeze(z_can, h, 1000, func_recirculation, h, w, lc);

  float result = Coeff * Integrale;

  return result;
}

float ComputeRcanInf_recirculation(float lc, float z_can, float w, float h, float U_h)
{
  float kappa = 0.41;

  float Coeff = 4. * w / ( pow(lc, 2) * pow(kappa, 2) * pow(h, 2) * U_h);
  float Integrale = Trapeze(z_can, h, 1000, func_recirculation, h, w, lc);

  float result = Coeff * Integrale;

  return result;
}

float ComputeRcanInf_recirculation_exp_profile(float lc, float z_can, float h,
					       float U_h, float WindProfile_param)
{
  float Integrale = U_h * Trapeze(z_can, h, n, Kt_recirculation, WindProfile_param, h, lc);

  return Integrale;
}
*/

float ComputeRcanInf_recirculation_exp_profile(float lc, float z_can, float h,
					       float U_h, float WindProfile_param)
{
	float Coeff = 1. / (WindProfile_param * U_h / h);

  float Integrale = Coeff * Trapeze(z_can, h, n, Kt_recirculation, WindProfile_param, h, lc);

  return Integrale;
}


/*****************************************************/
/*    Compute the resistance of below part of the    */
/*         canyon, in the ventilated region          */
/*                                                   */
/* param:  d        - displacement height            */
/*         X        - canyon's caracteristic length  */
/*         w        - canyon's mean width            */
/*         h        - buildings' mean height         */
/*         Ustar    - friction velocity              */
/*         FirstLevelWindModule                      */
/*                  - wind module at the first       */
/*                     vertical node                 */
/*         zref     - height of the atmospheric      */
/*                     model's first level           */
/*                                                   */
/* return: RcanInf  - aerodynamical resistance of    */
/*                     below part of the canyon      */
/*****************************************************/
/*
float ComputeRcanInf_ventilated(float d, float z_can, float w, float h, float zref,
				float Ustar, float FirstLevelWindModule)
{
  float kappa = 0.41;
  float Utop = ComputeUtop(d, h, zref, Ustar, FirstLevelWindModule);

  float Coeff = 4. * w / ( pow(kappa, 2) * pow(h, 2) * Utop );
  float Integrale = Trapeze(z_can, h, 1000, func_ventilated, h, w);

  float result = Coeff * Integrale;

  return result;
}

float ComputeRcanInf_ventilated(float z_can, float w, float h, float U_h)
{
  float kappa = 0.41;

  float Coeff = 4. * w / ( pow(kappa, 2) * pow(h, 2) * U_h );
  float Integrale = Trapeze(z_can, h, 1000, func_ventilated, h, w);

  float result = Coeff * Integrale;

  return result;
}

float ComputeRcanInf_ventilated_exp_profile(float z_can, float h, float U_h, float WindProfile_param)
{
  float Integrale = U_h * Trapeze(z_can, h, n, Kt_ventilated, WindProfile_param, h, NULL);

  return Integrale;
}

float ComputeRcanInf_ventilated_exp_profile(float z_can, float h, float U_h, float WindProfile_param)
{
  float Integrale = U_h * Trapeze(z_can, h, n, Kt_ventilated, WindProfile_param, h);

  return Integrale;
}
*/

float ComputeRcanInf_ventilated_exp_profile(float z_can, float h, float U_h, float WindProfile_param)
{
  float Coeff = 1. / (WindProfile_param * U_h / h);
  float Integrale = Coeff * Trapeze(z_can, h, n, Kt_ventilated, WindProfile_param, h);

  return Integrale;
}


/************************************************************/
/* Calcul numérique d'intégrale par la méthode des trapèzes */
/*                                                          */
/* INPUT:                                                   */
/*   borne_inf    - Borne inférieure de l'intervalle        */
/*   borne_sup    - Borne supérieure de l'intervalle        */
/*   n            - Nombre d'itérations                     */
/*   func         - Fonction à intégrer                     */
/*   arg1         - Premier argument de la fonction func    */
/*   arg2         - Second argument de la fonction func     */
/*   arg3         - Troisième argument de la fonction func  */
/* OUTPUT:                                                  */
/*     Valeur de l'intégrale prise entre les deux bornes    */
/************************************************************/
float Trapeze(float borne_inf, float borne_sup, int n, float (&func)(float,float,float,float),
	      float arg1, float arg2, float arg3)
{
  // Valeur de la base des trapèzes
  float k;
  // Valeur de l'aire sous la courbe
  float Int = 0.;

  // Initialisation
  k = (borne_sup - borne_inf) / n;
  Int = k * ( func(borne_inf, arg1, arg2, arg3) + func(borne_sup, arg1, arg2, arg3) ) / 2.;

  // Algorithme
  for(int i=1; i<=n-1; i++)
    Int += k * func(borne_inf+i*k, arg1, arg2, arg3);

  return Int;
}

float Trapeze(float borne_inf, float borne_sup, int n, float (&func)(float,float,float),
	      float arg1, float arg2)
{
  // Valeur de la base des trapèzes
  float k;
  // Valeur de l'aire sous la courbe
  float Int = 0.;

  // Initialisation
  k = (borne_sup - borne_inf) / n;
  Int = k * ( func(borne_inf, arg1, arg2) + func(borne_sup, arg1, arg2) ) / 2.;

  // Algorithme
  for(int i=1; i<=n-1; i++)
    Int += k * func(borne_inf+i*k, arg1, arg2);

  return Int;
}

/************************************************************/
/*          Déclaration de la fonction à intégrer           */
/* INPUT:                                                   */
/*   x            - Point auquel est calculée la fonction   */
/*   h            - Hauteur moyenne des bâtiments           */
/*   w            - Largeur moyenne de la voie              */
/*   lc           - Longueur de mélange de la canopée       */
/* OUTPUT:                                                  */
/*            Valeur de la fonction en un point             */
/************************************************************/
/*
float func_recirculation(float x, float h, float w, float lc)
{
  float kappa = 0.41;
  float result = pow( lc + kappa * x, 2.) * exp ( h * (h - x) / ( 4. * w * x) );

  return result;
}
*/
/************************************************************/
/*          Déclaration de la fonction à intégrer           */
/* INPUT:                                                   */
/*   x            - Point auquel est calculée la fonction   */
/*   h            - Hauteur moyenne des bâtiments           */
/*   w            - Largeur moyenne de la voie              */
/*   empty argument                                         */
/* OUTPUT:                                                  */
/*            Valeur de la fonction en un point             */
/************************************************************/
/*
float func_ventilated(float x, float h, float w, float null_argument)
{
  float result = exp ( h * (h - x) / ( 4. * w * x) );

  return result;
}

float Kt_recirculation(float x, float WindProfile_param, float h, float lc)
{
  // Mixing length
  float lm = ComputeLm(x, lc);

  // Turbulent eddy diffusivity
  float Kt = WindProfile_param * pow(lm, 2.) * exp ( WindProfile_param * ( x / h - 1. ) ) / h;
  return Kt;
}

*/

float Kt_recirculation(float x, float WindProfile_param, float h, float lc)
{
  // Mixing length
  float lm = ComputeLm(x, lc);

  // Turbulent eddy diffusivity
  float Kt = 1./ (pow(lm, 2.) * exp ( WindProfile_param * ( x / h - 1. ) ));
  return Kt;
}

/*
float Kt_ventilated(float x, float WindProfile_param, float h, float null_argument)
{
  // Von Karman constant
  float Ka = 0.41;

  // Mixing length
  float lm = Ka * x;

  // Turbulent eddy diffusivity
  float Kt = WindProfile_param * pow(lm, 2.) * exp ( WindProfile_param * ( x / h - 1. ) ) / h;
  return Kt;
}

float Kt_ventilated(float x, float WindProfile_param, float h)
{
  // Von Karman constant
  float Ka = 0.41;

  // Mixing length
  float lm = Ka * x;

  // Turbulent eddy diffusivity
  float Kt = WindProfile_param * pow(lm, 2.) * exp ( WindProfile_param * ( x / h - 1. ) ) / h;
  return Kt;
}
*/

float Kt_ventilated(float x, float WindProfile_param, float h)
{
  // Von Karman constant
  float Ka = 0.41;

  // Mixing length
  float lm = Ka * x;

  // Turbulent eddy diffusivity
  float Kt = 1./ (pow(lm, 2.) * exp ( WindProfile_param * ( x / h - 1. ) ));

  return Kt;
}

/*
float Kt_log(float x, float U_star)
{
  return 0;
}
*/
/*****************************************************/
/*       Compute the resistance of the surface       */
/*            in the recirculation region            */
/*                                                   */
/* param:  d        - displacement height            */
/*         lc       - average mixing length in the   */
/*                     canopy                        */
/*         z_can    - canyon's caracteristic length  */
/*         z_0      - local roughness height         */
/*         w        - canyon's mean width            */
/*         h        - buildings' mean height         */
/*         zref     - height of the atmospheric      */
/*                     model's first level           */
/*         Ustar    - friction velocity              */
/*         FirstLevelWindModule                      */
/*                  - wind module at the first       */
/*                     vertical node                 */
/*                                                   */
/* return: Rsurface - aerodynamical resistance of    */
/*                     the surface                   */
/*****************************************************/
/*
float ComputeRsurface_recirculation(float d, float lc, float z_can, float z_0, float w,
				    float h, float zref, float Ustar, float FirstLevelWindModule)
{
  float kappa = 0.41;
  float Utop = ComputeUtop(d, h, zref, Ustar, FirstLevelWindModule);

  float Coeff = 4. * w / ( pow(lc, 2) * pow(kappa, 2) * pow(h, 2) * Utop);
  float Integrale = Trapeze(z_0, z_can, 1000, func_recirculation, h, w, lc);

  float result = Coeff * Integrale;

  return result;
}

float ComputeRsurface_recirculation(float lc, float z_can, float z_0, float w, float h, float U_h)
{
  float kappa = 0.41;

  float Coeff = 4. * w / ( pow(lc, 2.) * pow(kappa, 2.) * pow(h, 2.) * U_h);
  float Integrale = Trapeze(z_0, z_can, 1000, func_recirculation, h, w, lc);

  float result = Coeff * Integrale;

  return result;
}

float ComputeRsurface_recirculation_exp_profile(float lc, float z_can, float z_0, float h,
						float U_h, float WindProfile_param, float alpha)
{
  float z_lim = ComputeZlim(lc, alpha);

  float Integrale = 0;

  if(z_lim > z_0)
    {
      float U_star_surface =  ComputeUstar_Surface(z_0, z_lim, h, WindProfile_param);
      Integrale += ComputeRa_z(z_0, z_lim, 0., U_star_surface);
    }

  Integrale += U_h * Trapeze(z_lim, z_can, n, Kt_recirculation, WindProfile_param, h, lc);

  return Integrale;
}

float ComputeRsurface_recirculation_exp_profile(float lc, float z_can, float z_0, float h,
						float U_h, float WindProfile_param, float alpha)
{
  float z_lim = ComputeZlim(lc, alpha);

  float Integrale = 0;

  if(z_lim > z_0)
    {
      float U_star_surface =  ComputeUstar_Surface(z_0, z_lim, h, WindProfile_param, U_h);
      Integrale += ComputeRa_z(z_0, z_lim, 0., U_star_surface);
    }

  Integrale += U_h * Trapeze(z_lim, z_can, n, Kt_recirculation, WindProfile_param, h, lc);

  return Integrale;
}
*/

float ComputeRsurface_recirculation_exp_profile(float lc, float z_can, float z_0, float h, float U_h, float WindProfile_param, float alpha)
{
  float z_lim = ComputeZlim(lc, alpha);

  float Integrale = 0.;
  float Coeff = 1. / (WindProfile_param * U_h / h);

  if(z_lim > z_0)
    {
      float U_star_surface =  ComputeUstar_Surface(z_0, z_lim, h, WindProfile_param, U_h);
      //cout<<"U_star_surface = "<<U_star_surface<<endl;
      Integrale += ComputeRa_z(z_0, z_lim, 0., U_star_surface);
    }

  Integrale += Coeff * Trapeze(z_lim, z_can, n, Kt_recirculation, WindProfile_param, h, lc);

  return Integrale;
}

/*****************************************************/
/*       Compute the resistance of the surface       */
/*             in the ventilated region              */
/*                                                   */
/* param:  d        - displacement height            */
/*         z_can    - canyon's caracteristic length  */
/*         z_0      - local roughness height         */
/*         w        - canyon's mean width            */
/*         h        - buildings' mean height         */
/*         zref     - height of the atmospheric      */
/*                     model's first level           */
/*         Ustar    - friction velocity              */
/*         FirstLevelWindModule                      */
/*                  - wind module at the first       */
/*                     vertical node                 */
/*                                                   */
/* return: Rsurface - aerodynamical resistance of    */
/*                     the surface                   */
/*****************************************************/
/*
float ComputeRsurface_ventilated(float d, float z_can, float z_0, float w, float h,
				 float zref, float Ustar, float FirstLevelWindModule)
{
  float kappa = 0.41;
  float Utop = ComputeUtop(d, h, zref, Ustar, FirstLevelWindModule);

  float Coeff = 1. / ( pow(kappa, 2.) * Utop );
  float Integrale = Trapeze(z_0, z_can, 1000, func_ventilated, h, w);

  float result = Coeff * Integrale;

  return result;
}

float ComputeRsurface_ventilated(float z_can, float z_0, float w, float h, float U_h)
{
  float kappa = 0.41;

  float Coeff = 1. / ( pow(kappa, 2.) * U_h );
  float Integrale = Trapeze(z_0, z_can, 1000, func_ventilated, h, w);

  float result = Coeff * Integrale;

  return result;
}
*/
/*
float ComputeRsurface_ventilated_exp_profile(float lc, float z_can, float z_0, float h,
					     float U_h, float WindProfile_param, float alpha)
{
  float z_lim = ComputeZlim(lc, alpha);

  float Integrale = 0;

  if(z_lim > z_0)
    {
      float U_star_surface =  ComputeUstar_Surface(z_0, z_lim, h, WindProfile_param);
      Integrale += ComputeRa_z(z_0, z_lim, 0., U_star_surface);
      Integrale += U_h * Trapeze(z_lim, z_can, n, Kt_ventilated, WindProfile_param, h, NULL);
    }
  else
    Integrale += U_h * Trapeze(z_0, z_can, n, Kt_ventilated, WindProfile_param, h, NULL);

  return Integrale;
}

float ComputeRsurface_ventilated_exp_profile(float lc, float z_can, float z_0, float h,
					     float U_h, float WindProfile_param, float alpha)
{
  float z_lim = ComputeZlim(lc, alpha);

  float Integrale = 0;

  if(z_lim > z_0)
    {
      float U_star_surface =  ComputeUstar_Surface(z_0, z_lim, h, WindProfile_param);
      Integrale += ComputeRa_z(z_0, z_lim, 0., U_star_surface);
      Integrale += U_h * Trapeze(z_lim, z_can, n, Kt_ventilated, WindProfile_param, h);
    }
  else
    Integrale += U_h * Trapeze(z_0, z_can, n, Kt_ventilated, WindProfile_param, h);

  return Integrale;
}

float ComputeRsurface_ventilated_exp_profile(float lc, float z_can, float z_0, float h,
					     float U_h, float WindProfile_param, float alpha)
{
  float z_lim = ComputeZlim(lc, alpha);

  float Integrale = 0;

  if(z_lim > z_0)
    {
      float U_star_surface =  ComputeUstar_Surface(z_0, z_lim, h, WindProfile_param, U_h);
      Integrale += ComputeRa_z(z_0, z_lim, 0., U_star_surface);
      Integrale += U_h * Trapeze(z_lim, z_can, n, Kt_ventilated, WindProfile_param, h);
    }
  else
    Integrale += U_h * Trapeze(z_0, z_can, n, Kt_ventilated, WindProfile_param, h);

  return Integrale;
}
*/
float ComputeRsurface_ventilated_exp_profile(float lc, float z_can, float z_0, float h,
					     float U_h, float WindProfile_param, float alpha)
{
  float z_lim = ComputeZlim(lc, alpha);

  float Integrale = 0.;
	float Coeff = 1. / (WindProfile_param * U_h / h);

  if(z_lim > z_0)
    {
//      float U_star_surface =  ComputeUstar_Surface(z_0, z_lim, h, WindProfile_param, U_h);
      float U_star_surface =  ComputeUstar_Surface(z_0, z_lim, h, WindProfile_param, U_h);
      Integrale += ComputeRa_z(z_0, z_lim, 0., U_star_surface);
      //Integrale += Coeff * Trapeze(z_lim, z_can, n, Kt_ventilated, WindProfile_param, h);
    }
  //else
    //Integrale += Coeff * Trapeze(z_0, z_can, n, Kt_ventilated, WindProfile_param, h);

  Integrale += Coeff * Trapeze(z_lim, z_can, n, Kt_ventilated, WindProfile_param, h);

  return Integrale;
}

/*********************************************************************************************************/
/*                                                                                                       */
/*                                      Turbulence scheme: k-l                                           */
/*                                                                                                       */
/*********************************************************************************************************/

/*****************************************************/
/*  Compute the resistance of below part of the      */
/*   canyon, in TKE model                            */
/*                                                   */
/* param:  d        - displacement height            */
/*         lc       - average mixing length in the   */
/*                     canopy                        */
/*         X        - canyon's caracteristic length  */
/*         z0       - roughness length of the town   */
/*         a        - empirical constant in          */
/*                                K = a sqrt(k) l    */
/*         TKE      - turbulent kinetic energy       */
/*                                                   */
/* return: RcanInf  - aerodynamical resistance of    */
/*                     below part of the canyon      */
/*                                                   */
/* remarque:        - TKE constante                  */
/*****************************************************/
/*
float ComputeRcanInf_TKE(float lc, float X, float h, float a, float TKE)
{
  // Von Karman constant
  float Ka = 0.41;

  float R = ( ( h - X ) / lc  + log( h / X ) / Ka ) / ( a * sqrt(TKE) );

  return R;
}
*/
/*****************************************************/
/*  Compute the resistance of the surface            */
/*                                                   */
/* param:  lc       - average mixing length in the   */
/*                     canopy                        */
/*         X        - canyon's caracteristic length  */
/*         z0       - roughness length of the town   */
/*         a        - empirical constant in          */
/*                                K = a sqrt(k) l    */
/*         TKE      - turbulent kinetic energy       */
/*                                                   */
/* return: Rstreet  - aerodynamical resistance of    */
/*                     the street/wall in the canyon */
/*****************************************************/
/*
float ComputeRsurface_TKE(float lc, float X, float z0, float a, float TKE)
{
  // Von Karman constant
  float Ka = 0.41;

  float R = ( X - z0 ) / ( a * lc * sqrt(TKE) )
    + log( X / z0 ) / ( Ka * a * sqrt(TKE) );
  return R;
}
*/


/******************************************************/
/* Return the number of urban LUC in either wesely or */
/*    zhang's LUC categories                          */
/* param:  landuse_type     - zhang or wesely         */
/*                                                    */
/* return: UrbanLUC         - number of urban LUC     */
/******************************************************/
/*
int ComputeUrbanLUC(std::string landuse_type)
{
  int UrbanLUC;
  if (landuse_type == "zhang")
    UrbanLUC = 15;
  else if (landuse_type == "wesely")
    UrbanLUC = 1;
  return UrbanLUC;
}
*/
