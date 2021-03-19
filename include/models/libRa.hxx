#ifndef LIBRA_HXX
#define LIBRA_HXX

#include<cmath>
#include<string>
//#include"libRa_Louis.hxx"
#include"libRa_Louis.cxx" // YK

const int n = 1000;

/*********************************************************************/
/*                  Defining characteristic lengths                  */
/*********************************************************************/
float ComputeD(float A, float lambda_p, float h);
// **** correction carole **** //
//float ComputeZ0(float lambda_p, float h, float d, float Cd);
float ComputeZ0(float lambda_f, float h, float d, float Cd);
//float computeZ0(float lambda_p, float h);
float computeZ0(float lambda_f, float h, float Cd);
// **** end correction carole **** //
float ComputeZcan(float w, float l);
float ComputeLc(float h, float d);
float ComputeLm(float z, float lc);
float ComputeZlim(float lc, float alpha);

/*********************************************************************/
/*                         Dynamic properties                        */
/*********************************************************************/
float ComputeUstar_AboveCanopy(float zref, float z0, float d,
			       float FirstLevelWindModule);
float ComputeUstar_AboveCanopy(float zref, float z0, float d,
			       float Richardson, float FirstLevelWindModule);
float ComputeUstar_InsideCanopy(float lc, float h, float Utop, float z0);
float ComputeUstar_Surface(float z0, float zlim, float h, float a);
// **** correction carole **** //
float ComputeUstar_Surface(float z0, float zlim, float h, float a, float Utop);
// **** end correction carole **** //
float ComputeUtop(float zref, float h, float d, float Ustar,
                  float FirstLevelWindModule);
// **** correction carole **** //
float ComputeUtop(float zref, float h, float d, float Ustar,
                  float FirstLevelWindModule, float z0);
// **** end correction carole **** //
float ComputeUtop_CanyonIntegration(float w, float h, float w_r);

/*********************************************************************/
/*                   Turbulence scheme: Common path                  */
/*********************************************************************/
float ComputeRa_z(float z, float zref, float d, float Ustar);
float ComputeRcanSup(float zref, float h, float d, float Ustar);

/*********************************************************************/
/*            Turbulence scheme: Improved Prandtl's model            */
/*********************************************************************/
float ComputeRcanInf_recirculation(float d, float lc, float z_can, float w,
				   float h, float zref, float Ustar,
				   float FirstLevelWindModule);
float ComputeRcanInf_recirculation(float lc, float z_can, float w,
				   float h, float U_h);
float ComputeRcanInf_recirculation_exp_profile(float lc, float z_can, float h,
					       float U_h, float a);
//float ComputeRcanInf_ventilated(float d, float lc, float z_can, float w,
//				float h, float zref, float Ustar,
//				float FirstLevelWindModule);
float ComputeRcanInf_ventilated(float d, float z_can, float w, float h,
				float zref, float Ustar,
				float FirstLevelWindModule);
float ComputeRcanInf_ventilated(float z_can, float w, float h, float U_h);
float ComputeRcanInf_ventilated_exp_profile(float z_can, float h, float U_h,
					    float a);
float Trapeze(float borne_inf, float borne_sup, int n,
	      float (&func)(float,float,float,float),
	      float arg1, float arg2, float arg3=0);
float Trapeze(float borne_inf, float borne_sup, int n,
	      float (&func)(float,float,float),
	      float arg1, float arg2);
float func_recirculation(float x, float h, float w, float lc);
float func_ventilated(float x, float h, float w, float null_argument=0);
float Kt_recirculation(float x, float a, float h, float lc);
//float Kt_ventilated(float x, float a, float h, float null_argument);
float Kt_ventilated(float x, float a, float h);
float Kt_log(float x, float U_star);
float ComputeRsurface_recirculation(float d, float lc, float z_can, float z_0,
				    float w, float h, float zref, float Ustar,
				    float FirstLevelWindModule);
float ComputeRsurface_recirculation(float lc, float z_can, float z_0,
				    float w, float h, float U_h);
float ComputeRsurface_recirculation_exp_profile(float lc, float z_can, float z_0, float h,
						float U_h, float a, float alpha);
float ComputeRsurface_ventilated(float d, float z_can, float z_0, float w,
				 float h, float zref, float Ustar,
				 float FirstLevelWindModule);
float ComputeRsurface_ventilated(float z_can, float z_0, float w,
				 float h, float U_h);
float ComputeRsurface_ventilated_exp_profile(float lc, float z_can, float z_0, float h,
					     float U_h, float a, float alpha);

/*********************************************************************/
/*                    Turbulence scheme: k-l                         */
/*********************************************************************/
float ComputeRcanInf_TKE(float lc, float X, float h, float a, float TKE);
float ComputeRsurface_TKE(float lc, float X, float z0, float a, float TKE);
int ComputeUrbanLUC(std::string landuse_type);


#endif
