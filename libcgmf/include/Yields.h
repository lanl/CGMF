/*------------------------------------------------------------------------------
  CGMF-1.0
  Copyright TRIAD/LANL/DOE - see file COPYRIGHT.md
  For any questions about CGMF, please contact us at cgmf_help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file Yields.h

  \brief Systematics for scission fragment yields

*/

#ifndef __YIELDS_H__
#define __YIELDS_H__

#include <string>

class Yields{
    public:
    	Yields(int id,bool sf_flag_in);
		~Yields();
		void setYieldParameters(void);
		void setTKEParameters(void);
		double ViolaSyst(int id);
		void setFissioningSystem(int id,double ex);
		double yieldA(int A);
		double yieldATKE(int A,double TKE);
		double yieldTKE(int A,double TKE);
		double sampleTKE(int A);
		double get_avTKE(int A);
		void setRescaleTKE(double *,int,int);
		void LoadYAParams(std::string, int idx);
		void LoadTKEParams(std::string, int idx);
		void WriteTKEA(int Amin, int Amax);
    private:
		// Number of Gaussian modes (3)
		int numModes;
		int id_cn0;
		int Acn0,Zcn0;
		int a_cn;

		// For the 5-Gaussian fit (two Gaussians are degenerate)
		double w_e[5],sig_e[5],Abar[5];

		// There are now 4 sets of Gaussian parameters (for 1st, 2nd, 3rd, and 4th chance fission)
		// Each parameter 5 values (one for each Gaussian -- we note that the light and heavy fragment Gaussians are degenerate)
		double w_a[4][5], w_b[4][5], w_c[4][5], sig_a[4][5], sig_b[4][5], sig_c[4][5], mu_a[4][5], mu_b[4][5], mu_c[4][5];
		bool NFsupp_YA[4]; // Lists if there are Gaussian parameters for each chance fission

		// For the <TKE>(En), <TKE>(Ah), and s_TKE(Ah) fits
		double TKE_fit_param[2];
		double TKEA_fit_param[11];
		double sTKEA_fit_param[11];

		// There are 5 parameters for the <TKE>(En) parameterization -- assumes linear relation (4 sets for 1st, 2nd, 3rd, and 4th chance fission)
		double tke_params[4][5];
		bool NFsupp_TKE[4];

		// There are 11 parameters for the power-law fits of <TKE>(Ah) and sigma_TKE(Ah) (9 coefficients, the A-value we expand about, and the maximum A)
		double tkea_coeff[4][11];
		double stkea_coeff[4][11];

		double gaussian(double x, double w, double av, double s);
		double avTKE[300];
		double sTKE[300];
		bool sf_flag;
		void setParametersAverageTKE(void);
		double averageTKE(double e); // Calculates average TKE for the system as a function of incident energy
		double rescaleAverageTKE;
		double sampleGaussian(void);
		double einc_now;
};

#endif //__YIELDS_H__


