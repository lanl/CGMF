/*------------------------------------------------------------------------------
  CGMF-1.1
  Copyright TRIAD/LANL/DOE - see file LICENSE
  For any questions about CGMF, please contact us at cgmf-help@lanl.gov
-------------------------------------------------------------------------------*/

/*! @file yields.cpp

  \brief Fission yields

*/

#include <string.h>
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <stdlib.h>
#include <iterator>
#include <functional>

#include "FissionFragments.h"
#include "Yields.h"
#include "physics.h"
#include "config-ff.h"
#include "config.h"
#include "rngcgm.h"

using namespace std;

/****************************************************************
*   Yields.cpp is responsible for reading in the data files	*
*   containing the parameterization data for the Y(A,TKE). The	*
*   Y(A), <TKE>(En), <TKE>(Ah), and s_TKE(Ah) are read from	*
*   YAModel.dat and TKEModel.dat. These parameterizations can	*
*   be adjusted to fit whatever data the user wants and none	*
*   of this information is hard-coded. This file stores these	*
*   parameterizations (up to 4th chance fission) and samples	*
*   from the Y(A) or <TKE>(Ah) for a particular fission event.	*
*   The <TKE>(Ah) is renormalized such that the weighted sum:	*
*   \Sum [ <TKE>(Ah) * Y(Ah|En) , Ah ] = <TKE>(En)		*
*   We note that the Y(A) have an energy-dependence and the	*
*   <TKE>(Ah) will be scaled based on the energy, but no shape	*
*   change is assumed for <TKE>(Ah) or s_TKE(Ah).		*
****************************************************************/

/*******************************************************************************
 * Yields constructor
 *------------------------------------------------------------------------------
 * Defines and initializes the yields given a compound nucleus id and whether
 * this nucleus is spontaneously fissioning or not
 ******************************************************************************/
Yields::Yields(int id, bool sf_flag_in) {

    id_cn0  = id; // Defines the ID within Yields
    Zcn0    = id/1000;
    Acn0    = id%1000;
    sf_flag = sf_flag_in; // and the sf_flag within Yields

    setYieldParameters();	// Load the Y(A) parameterization
    setTKEParameters();		// Load the <TKE> parameterizations

    return;
}

// Deconstructor
Yields::~Yields(void) { }

/*******************************************************************************
 * LoadYAParams
 *------------------------------------------------------------------------------
 * Reads and assigns the Y(A) parameterizations from the input file.
 ******************************************************************************/
void Yields::LoadYAParams(string line, int idx) {

    int toss;
    if (idx>4) { // Only support up to 4th chance fission
    	return;
    }

    std::istringstream iss(line);
    iss >> toss; // The ZAID

    iss >>   w_a[idx][1]; // 1st peak amplitude
    iss >>   w_b[idx][1];
    iss >>  mu_a[idx][1]; // 1st peak mean
    iss >>  mu_b[idx][1];
    iss >> sig_a[idx][1]; // 1st peak variance
    iss >> sig_b[idx][1];

    iss >>   w_a[idx][2]; // 2nd peak amplitude
    iss >>   w_b[idx][2];
    iss >>  mu_a[idx][2]; // 2nd peak mean
    iss >>  mu_b[idx][2];
    iss >> sig_a[idx][2]; // 2nd peak variance
    iss >> sig_b[idx][2];

    iss >> sig_a[idx][0]; // Symmetric peak variance
    iss >> sig_b[idx][0];

    NFsupp_YA[idx]=true; // Mark that idx-chance fission is supported

/* For debugging
    cout<<"Mode 0"<<endl;
    cout<<"sig(E) = "<<sig_a[idx][0]<<" + "<<sig_b[idx][0]<<"E"<<endl;

    for (unsigned int i=1;i<3;i++) {
	cout<<"Mode "<<i<<endl;
	cout<<"W(E) = 1/(1+exp[("<<w_a[idx][i]<<"+E)/"   <<w_b[idx][i]<<"])"<<endl;
	cout<<"mu(E) = "          <<mu_a[idx][i]<<" + "<< mu_b[idx][i]<<"E" <<endl;
	cout<<"sig(E) = "        <<sig_a[idx][i]<<" + "<<sig_b[idx][i]<<"E" <<endl;
    }
*/
    return;
}

/*******************************************************************************
 * setYieldParameters
 *------------------------------------------------------------------------------
 * Initial read of the Y(A) parameterizations. Based on the compound nucleus id
 * and the sf_flag, we read the appropriate Y(A) parameters. For neutron-induced
 * fission, we store all Y(A) up to 4th-chance fission.
 ******************************************************************************/
void Yields::setYieldParameters(void) {

    numModes = 3; // Three (independent) Gaussians

    string GaussData = datadir;
    GaussData += MASSYIELDFILE; // 3-Gaussian parameterizations
    ifstream gaussfit(GaussData.c_str(),ios::in);
    if (!gaussfit.is_open()) {
    	cerr << "The file " << GaussData << " was not found" << endl; exit(0);
    }

    // Initialize the arrays
    for (unsigned int i=0;i<4;i++) NFsupp_YA[i] = false; // Default is to do a simple Y(A) shift if the N-chance yields aren't available (see below)
    bool not_found = true;
    int this_id;
    char c[300];
    string line, subline;
    
    gaussfit.getline(c,299);
    while (!gaussfit.eof()) { // Loop through the YAModel file and input the necessary info
        if (c[0] == '#') {gaussfit.getline(c,299); continue;} // Ignore comment lines

        line = (string)c;
        std::istringstream iss(line);
        int Nparam = std::distance(std::istream_iterator<string>(iss),std::istream_iterator<string>());
//        if (Nparam!=15) {gaussfit.getline(c,299); continue;} // Ignore any lines without the correct number of parameters

        // Determine the Z and A for this data
        iss.str(line);
        iss.clear();
        iss >> this_id;
        if (!sf_flag) {
            // Accept the ID of the compound nucleus, or its 2nd, 3rd, and 4th - chance daughters
            if (this_id == id_cn0) { LoadYAParams(line,0); not_found=false;}
            if (this_id == id_cn0-1) LoadYAParams(line,1);
            if (this_id == id_cn0-2) LoadYAParams(line,2);
            if (this_id == id_cn0-3) LoadYAParams(line,3);
        } else {
            // For SF, we only accept the compound nucleus (-ZAID means spontaneous fission)
            if (-this_id == id_cn0) { LoadYAParams(line,0); not_found=false; break;}
        }
        gaussfit.getline(c,299);
    }
    // Check if a Y(A) model exists for the initial compound nucleus
    if (not_found) {
        cerr << "No Y(A) model found for ";
        if (sf_flag) {
            cerr << "(sf) ";
        } else {
            cerr << "(n,f) ";
        }
        cerr << id_cn0 << endl;
        exit(0);
    }
    gaussfit.close();
    return;
}

/*******************************************************************************
 * LoadTKEParams
 *------------------------------------------------------------------------------
 * Reads and assigns the TKE parameterizations from the input file.
 ******************************************************************************/
void Yields::LoadTKEParams(string line, int idx) {

    int toss;
        if (idx>4) { // Only support up to 4th chance fission
	return;
    }

    std::istringstream iss(line);
    iss >> toss; // The ZAID
    if (toss<0) {
        toss = -toss;
    }

    /* THE <TKE>(EN) PARAMETERIZATION */
    iss >> tke_params[idx][0]; // The <TKE>(En=0) value
    iss >> tke_params[idx][4]; // The inflection point (i.e. if the <TKE>(En) changes from positive to negative slope as in U235,U238
    iss >> tke_params[idx][1]; // The positive slope
    iss >> tke_params[idx][3]; // The negative slope

    // Calculate the "zero-point" for the negative slope region: f(x) = a + bx and g(x) = c + dx (solve for c)
    tke_params[idx][2] = tke_params[idx][0] + (tke_params[idx][1] - tke_params[idx][3])*tke_params[idx][4];

    /* THE <TKE>(AH) PARAMETERIZATION -- NOTE: NO DEPENDENCE ON INCIDENT NEUTRON ENERGY YET!! */
    iss >> tkea_coeff[idx][0]; // The A value we expand around
    iss >> tkea_coeff[idx][1]; // The maximum A value the power-law fit should sample
    iss >> tkea_coeff[idx][2]; // c_0
    iss >> tkea_coeff[idx][3]; // c_1
    iss >> tkea_coeff[idx][4]; // c_2
    iss >> tkea_coeff[idx][5]; // c_3
    iss >> tkea_coeff[idx][6]; // c_4
    iss >> tkea_coeff[idx][7]; // c_5
    iss >> tkea_coeff[idx][8]; // c_6
    iss >> tkea_coeff[idx][9]; // c_7
    iss >> tkea_coeff[idx][10]; // c_8

    /* THE s_TKE(AH) PARAMETERIZATION -- NOTE: NO DEPENDENCE ON INCIDENT NEUTRON ENERGY YET!! */
    iss >> stkea_coeff[idx][0]; // The A value we expand around
    iss >> stkea_coeff[idx][1]; // The maximum A value the power-law fit should sample
    iss >> stkea_coeff[idx][2]; // c_0
    iss >> stkea_coeff[idx][3]; // c_1
    iss >> stkea_coeff[idx][4]; // c_2
    iss >> stkea_coeff[idx][5]; // c_3
    iss >> stkea_coeff[idx][6]; // c_4
    iss >> stkea_coeff[idx][7]; // c_5
    iss >> stkea_coeff[idx][8]; // c_6
    iss >> stkea_coeff[idx][9]; // c_7
    iss >> stkea_coeff[idx][10]; // c_8

    NFsupp_TKE[idx]=true; // Mark that idx-chance fission is supported

    int Ac = toss % 1000;

/* For debugging
    cout<<"<TKE>(En) = "<<tke_params[idx][0]<<" + En*("<<tke_params[idx][1]<<") for En<="<<tke_params[idx][4]<<endl;
    cout<<"<TKE>(En) = "<<tke_params[idx][2]<<" + En*("<<tke_params[idx][3]<<") for En>="<<tke_params[idx][4]<<endl;

    cout<<"<TKE>(Ah) = "<<tkea_coeff[idx][2];
    for (unsigned int i=1;i<9;i++) {
	cout<<" + ("<<tkea_coeff[idx][i+2]<<")(Ah - "<<tkea_coeff[idx][0]<<")^"<<i;
    }
    cout<<" for "<<Ac/2<<" <= Ah <= "<<tkea_coeff[idx][1]<<endl;

    cout<<"s_TKE(Ah) = "<<stkea_coeff[idx][2];
    for (unsigned int i=1;i<9;i++) {
	cout<<" + ("<<stkea_coeff[idx][i+2]<<")(Ah - "<<stkea_coeff[idx][0]<<")^"<<i;
    }
    cout<<" for "<<Ac/2<<" <= Ah <= "<<stkea_coeff[idx][1]<<endl;
*/
    return;
}

/*******************************************************************************
 * setTKEParameters
 *------------------------------------------------------------------------------
 * Initial read of the TKE parameterizations. Based on the compound nucleus id
 * and the sf_flag, we read the appropriate TKE parameters. For neutron-induced
 * fission, we store all TKE up to 4th-chance fission.
 ******************************************************************************/
void Yields::setTKEParameters(void) {

    string TKEData = datadir;
    TKEData += TKEYIELDFILE; // List of <TKE>(En), <TKE>(Ah), and s_TKE(Ah) parameterizations
    ifstream tkefit(TKEData.c_str(),ios::in);
    if (!tkefit.is_open()) {
    	cerr << "The file " << TKEData << " was not found" << endl; exit(0);
    }

    // Initialize the arrays
    for (unsigned int i=0;i<4;i++) NFsupp_TKE[i] = false; // If the <TKE> parameterizations aren't available for our nucleus, we use the systematics of the initial compound nucleus
    bool not_found = true;
    int this_id;
    char c[400];
    string line, subline;
    
    tkefit.getline(c,399);
    while (!tkefit.eof()) { // Loop through the TKEModel file and input the necessary info
	if (c[0] == '#') {tkefit.getline(c,399); continue;} // Ignore comment lines

	line = (string)c;
	std::istringstream iss(line);
	int Nparam = std::distance(std::istream_iterator<string>(iss),std::istream_iterator<string>());
	if (Nparam!=27) {tkefit.getline(c,399); continue;} // Ignore any lines without the correct number of parameters

	// Determine the Z and A for this data
	iss.str(line);
	iss.clear();
	iss >> this_id;
	if (!sf_flag) {
	    // Accept the ID of the compound nucleus, or its 2nd, 3rd, and 4th - chance daughters
	    if (this_id == id_cn0) {LoadTKEParams(line,0); not_found=false;}
	    if (this_id == id_cn0-1) LoadTKEParams(line,1);
	    if (this_id == id_cn0-2) LoadTKEParams(line,2);
	    if (this_id == id_cn0-3) LoadTKEParams(line,3);

	} else {
	    // For SF, we only accept the compound nucleus
	    if (this_id == -id_cn0) {LoadTKEParams(line,0); not_found=false; break;}
	}
	tkefit.getline(c,399);
    }
    tkefit.close();
    return;
}

/*******************************************************************************
 * ViolaSyst
 *------------------------------------------------------------------------------
 * Viola systematics for <TKE> when the compound nucleus is not supported
 ******************************************************************************/
double Yields::ViolaSyst(int id) {

    // Form the compound Z and A
    int Zcn = id/1000;
    int Acn = id%1000;
    // See <https://journals.aps.org/prc/pdf/10.1103/PhysRevC.31.1550> for background
    double x = Zcn*Zcn / pow((double)Acn,1.0/3.0);

    double Vsyst = 0.1189*x + 7.3;

    // Isotopic scaling -- this is done because we know that the Viola systematics doesn't perform very well along an isotopic chain
    // We multiply the Viola systematics <TKE>(x) = 0.1189*x + 7.3 by a scaling function f(x) = is_a*x + is_b where x = Z^2 / A^1/3
    double is_a = 0.;
    double is_b = 1.; // Defaults to Viola systematics (i.e. no isotopic scaling)

    switch (Zcn) {
	case 92: /* URANIUMS */
	    is_a = -2.43783800e-03;
	    is_b = 4.32954561;
	    break;
	case 93: /* NEPTUNIUMS */
	    is_a = -2.47262338e-03;
	    is_b = 4.46218351;
	    break;
	case 94: /* PLUTONIUMS */
	    is_a = -1.78668219e-03;
	    is_b = 3.55865245;
	    break;
	case 98: /* CALIFORNIUMS */
	    is_a = -3.44312133e-03;
	    is_b = 6.20846742;
	    break;
	default:
	    break;
    }
    double iso_scaling = is_a*x + is_b;
    return Vsyst*iso_scaling;
}

/*******************************************************************************
 * setFissioningSystem
 *------------------------------------------------------------------------------
 * Defines the fissioning system for a given fission event based on the nucleus
 * id and the "effective" incident neutron energy. This is "effective" when we
 * have multi-chance fission.
 ******************************************************************************/
void Yields::setFissioningSystem(int id, double einc) {
   
    // Determine the number of pre-fission neutrons and the effective incident neutron energy
    int emitted_neutrons = id_cn0 - id;
    int idx = emitted_neutrons; // Index for multi-chance fission (0,1,2,3 corresponds to 1st, 2nd, 3rd, 4th chance fission)
    a_cn = id%1000;
    einc_now = einc;
    int a0_cn = a_cn; // Will define the shift between the a_cn and the a0_cn (if the multi-chance fission mode isn't supported)

    /**************** SET THE TKE PARAMETERS FOR THIS FISSION EVENT ****************/

    /*  THE <TKE>(EN) PARAMETERIZATION USES TWO FUNCTIONS F(X) AND G(X) TO CAPTURE THE <TKE>(EN) SHAPE:

	F(X) = A + BX		FOR X<=E0
	G(X) = C + DX		FOR X>=E0

	WHERE E0 IS CALLED THE INFLECTION POINT. THIS MODELS THE CASES WHEN THE <TKE> APPEARS TO RISE
	FOR LOW EN VALUES AND THEN DECREASE AFTER A CERTAIN INCIDENT NEUTRON ENERGY. THIS IS SEEN IN
	U235(N,F) AND U238(N,F) AND NP237(N,F) DATA.
    */

    // Default is to use the Viola systematics and a gradual slope (when we're above 4th-chance fission or when no systematics exist)
    double E0 = 0.;		// No inflection point
    double A  = ViolaSyst(id);	// MeV
    double B  = -0.25;		// MeV/MeV
    double C  = ViolaSyst(id);	// MeV
    double D  = -0.25;		// MeV/MeV

    // If we're looking at 1st, 2nd, 3rd, or 4th-chance fission and a systematics exist for that compound nucleus, use it
    if ((idx<4) && (NFsupp_TKE[idx])) {
    	E0 = tke_params[idx][4]; // The inflection point (E0)
    	A = tke_params[idx][0]; // The zero-energy <TKE> value for En <= E0 region
    	B = tke_params[idx][1]; // The slope for En <= E0 region
    	C = tke_params[idx][2]; // The "zero-energy" <TKE> value for En >= E0 region
    	D = tke_params[idx][3]; // The slope for En >= E0 region
    } else {
    	idx = 0; // Default to the id_cn0 compound nucleus for the <TKE>(A) and s_TKE(A) systematics
    }

    // Set the TKE fit parameters
    if (einc<=E0) {
    	TKE_fit_param[0] = A;
    	TKE_fit_param[1] = B;
    } else {
    	TKE_fit_param[0] = C;
    	TKE_fit_param[1] = D;
    }

    // Set the <TKE> fit (see setParametersAverageTKE for use)
    TKEA_fit_param[0]  = tkea_coeff[idx][0]; // Expansion point
    TKEA_fit_param[1]  = tkea_coeff[idx][1]; // Maximum Ah
    TKEA_fit_param[2]  = tkea_coeff[idx][2];
    TKEA_fit_param[3]  = tkea_coeff[idx][3];
    TKEA_fit_param[4]  = tkea_coeff[idx][4];
    TKEA_fit_param[5]  = tkea_coeff[idx][5];
    TKEA_fit_param[6]  = tkea_coeff[idx][6];
    TKEA_fit_param[7]  = tkea_coeff[idx][7];
    TKEA_fit_param[8]  = tkea_coeff[idx][8];
    TKEA_fit_param[9]  = tkea_coeff[idx][9];
    TKEA_fit_param[10] = tkea_coeff[idx][10];

    // And the s_TKE fit (see setParametersAverageTKE for use)
    sTKEA_fit_param[0]  = stkea_coeff[idx][0]; // Expansion point
    sTKEA_fit_param[1]  = stkea_coeff[idx][1]; // Maximum Ah
    sTKEA_fit_param[2]  = stkea_coeff[idx][2];
    sTKEA_fit_param[3]  = stkea_coeff[idx][3];
    sTKEA_fit_param[4]  = stkea_coeff[idx][4];
    sTKEA_fit_param[5]  = stkea_coeff[idx][5];
    sTKEA_fit_param[6]  = stkea_coeff[idx][6];
    sTKEA_fit_param[7]  = stkea_coeff[idx][7];
    sTKEA_fit_param[8]  = stkea_coeff[idx][8];
    sTKEA_fit_param[9]  = stkea_coeff[idx][9];
    sTKEA_fit_param[10] = stkea_coeff[idx][10];

    /***************** SET THE Y(A) PARAMETERS FOR THIS FISSION EVENT *****************/

    idx = emitted_neutrons; // Reset the idx in case the Y(A) parameterization is supported for this idx-chance daughter nucleus when the <TKE> parameterization was not

    /*  THE Y(A) PARAMETERIZATION USES THREE GAUSSIANS WHOSE MEANS, VARIANCES, AND AMPLITUDES CAN
	HAVE AN ENERGY DEPENDENCE:

	W(E) = 1 / (1 + EXP[(E-W_A)/W_B])
	MU(E) = MU_A + MU_B * (E)
	SIG(E) = SIG_A + SIG_B * (E)
    */

    // If the N-chance value is above 4th chance or it's not supported, do a simple mass shift (see below)
    if (idx>4 || !NFsupp_YA[idx]) {idx = 0; a0_cn = id_cn0%1000;}

    // Calculate the Gaussian parameters for our fissioning nucleus at an "effective" incident neutron energy einc
    for (unsigned int i=1;i<numModes;i++) {
    	// Amplitudes are Fermi functions in incident neutron energy
    	w_e[i] = 1.0/(1.0 + exp((einc - w_a[idx][i])/w_b[idx][i]));

    	// Variances are linear in incident neutron energy
    	sig_e[i] = sig_a[idx][i] + sig_b[idx][i]*einc;

    	// Means are linear in incident neutron energy (with an optional shift if we don't have the appropriate Y(A) for this idx-chance daughter)
    	Abar[i] = (a_cn - a0_cn) + mu_a[idx][i] + mu_b[idx][i]*einc;

    	// We note that if the n-chance fission mode is supported, then a_cn = a0_cn and we will use the Abar
    	// from the 3-Gaussian fit. However, if the n-chance fission mode isn't supported we simply use the
    	// Y(A) distribution from the initial compound nucleus (idx=0) but shift it by the mass difference
    	// (a_cn - a0_cn).

    	// Ex: Let's say we're doing U238(n,f), so we start with U239 as the compound system. We have a Y(A)
    	// for U239, but not for U238. For 2nd chance fission, the a_cn = 238 and idx = 1, but
    	// NFsupp_YA[idx] = false. Therefore, we simply elect to use the idx = 0 yields (U239) with a a0_cn = 239.
    	// Then, when computing the mean Abar, we shift the value by (a_cn - a0_cn) = -1 in order to
    	// arbitrarily shift the Y(A) of U239 to a distribution that conserves mass. Hopefully, this makes
    	// sense but email P. Jaffke <pjaffke@lanl.gov> if not.
    }

    // The symmetric region is a bit different
    w_e[0]   = 2. - 2.*w_e[1] - 2.*w_e[2]; // Using 2*A1 + 2*A2 + A0 = 2
    sig_e[0] = sig_a[idx][0] + sig_b[idx][1]*einc;
    Abar[0]  = a_cn/2.; // The symmetric mean is always a_cn/2 regardless of whether or not the n-chance mode is supported

/* For debugging
    printf("Gaussian Parameters for En = %10.4f MeV:\n",einc);
    printf("w_e   = [%10.7f , %10.7f , %10.7f]\n",w_e[2],w_e[1],w_e[0]);
    printf("sig_e = [%10.4f , %10.4f , %10.4f]\n",sig_e[2],sig_e[1],sig_e[0]);
    printf("Abar  = [%10.4f , %10.4f , %10.4f]\n",Abar[2],Abar[1],Abar[0]);
*/
    return;
}

/*******************************************************************************
 * setRescaleTKE
 *------------------------------------------------------------------------------
 * Determines the rescaling factor to ensure \Sum[ Y(Ah)*<TKE>(Ah) ] = <TKE>
 ******************************************************************************/
void Yields::setRescaleTKE(double *y, int Amin, int Amax) {

    // Set the avTKE array to the current nucleus
    setParametersAverageTKE();
    double sum = 0.;
    for (int a=Amin; a<=Amax; a++) {
    	sum += avTKE[a]*y[a];
    }
    // Set the rescaling value such that \Sum [ Y(Ah)*<TKE>(Ah) ] = <TKE>
    rescaleAverageTKE = averageTKE(einc_now)/sum;
    return;
}

/*******************************************************************************
 * yieldA
 *------------------------------------------------------------------------------
 * Returns the mass yield Y(A) for a fragment mass A
 ******************************************************************************/
double Yields::yieldA(int A) {
    double a;
    if (A>=a_cn/2) {
    	a = (double)A;
    } else {
	   a = (double)(a_cn-A); // Y_pre(A) are symmetric across a_cn/2
    }
   
    double sum=gaussian(a, w_e[0], Abar[0], sig_e[0]); // Calculate the yield from the symmetric Gaussian
    for (int m=1; m<numModes; m++) {
    	sum += gaussian(a,w_e[m],Abar[m],sig_e[m]);
    	sum += gaussian(a,w_e[m],a_cn-Abar[m],sig_e[m]); // Sum the yield from all other Gaussians (light and heavy)
    }

    if (sum>1e-5) return sum;

    // If the yield is less than a cutoff we neglect it
    return 0.;
}

/*******************************************************************************
 * yieldATKE
 *------------------------------------------------------------------------------
 * Returns the yield Y(A,TKE) for a given A and TKE value
 ******************************************************************************/
double Yields::yieldATKE(int A, double TKE) {

    int a;
    if (A>=a_cn/2) {
    	a = A;
    } else {
	   a = a_cn - A;
    }
    double yatke = gaussian(TKE,yieldA(a),rescaleAverageTKE*avTKE[a],sTKE[a]);
    return yatke;
}

/*******************************************************************************
 * yieldTKE
 *------------------------------------------------------------------------------
 * Returns the yield Y(TKE|A) for a given A and TKE value
 ******************************************************************************/
double Yields::yieldTKE(int A, double TKE) {

    int a;
    if (A>=a_cn/2) {
	   a = A;
    } else {
	   a = a_cn - A;
    }
    double ytke = gaussian(TKE,1.0,rescaleAverageTKE*avTKE[a],sTKE[a]);
    if(ytke>1.e-7) return ytke;

    // If the yield is less than a cutoff we neglect it
    return 0.;
}

/*******************************************************************************
 * gaussian
 *------------------------------------------------------------------------------
 * Returns the value of a Gaussian with weight w, mean mu, and width sig at a
 * given x-value
 ******************************************************************************/
double Yields::gaussian(double x, double w, double mu, double sig) {

    double ct = w/(sqrt(PI2)*sig);
    double g = exp(-.5*(x - mu)*(x - mu)/sig/sig);
    return (ct*g);
}

/*******************************************************************************
 * setParametersAverageTKE
 *------------------------------------------------------------------------------
 * Uses the TKE parameters that have been stored by setFissioningSystem to
 * construct the <TKE>(Ah) and s_TKE(Ah). Then ensures that the weighted sum
 * \Sum [ Y(Ah)*<TKE>(Ah) ] = <TKE>.
 ******************************************************************************/
void Yields::setParametersAverageTKE(void) {

    // Default values for the <TKE>(Ah) and s_TKE(Ah)
    fill_n(avTKE, 300, 140.);
    fill_n(sTKE, 300, 7.);

    // Fill the values of <TKE>(Ah) and s_TKE(Ah)

    /*	BOTH THE <TKE>(AH) AND S_TKE(AH) ARE POWER-LAW EXPANSIONS ABOUT A GIVEN MASS

	<TKE>(AH) = C0 + C1*(A-A0) + C2*(A-A0)^2 + C3*(A-A0)^3 + ... + C8*(A-A0)^8
	s_TKE(AH) = C0 + C1*(A-A0) + C2*(A-A0)^2 + C3*(A-A0)^3 + ... + C8*(A-A0)^8

	EACH WITH A SET OF COEFFICIENTS C_I AND A EXPANSION POINT A0
    */

    // The parameterizations are designed to fit data up to a certain mass value: param[1]
    for (int a=(a_cn+1)/2; a<=(int)TKEA_fit_param[1]; a++) {
    	avTKE[a] = 0.;
    	for (int j=0;j<9;j++) avTKE[a] += TKEA_fit_param[j+2]*pow(a-TKEA_fit_param[0],(double)j);
    	avTKE[a_cn-a] = avTKE[a]; // Symmetrize
    }
    for (int a=(a_cn+1)/2; a<=(int)sTKEA_fit_param[1]; a++) {
    	sTKE[a] = 0.;
    	for (int j=0;j<9;j++) sTKE[a] += sTKEA_fit_param[j+2]*pow(a-sTKEA_fit_param[0],(double)j);
    	sTKE[a_cn-a] = sTKE[a]; // Symmetrize
    }
    return;
}

/*******************************************************************************
 * WriteTKEA
 *------------------------------------------------------------------------------
 * Prints the <TKE>(A) to shell
 ******************************************************************************/
void Yields::WriteTKEA(int Amin, int Amax) {

    for (unsigned int i=Amin; i<=Amax; i++) {
    	printf("%3d %7.3f\n",i,avTKE[i]);
    }
    return;
}

/*******************************************************************************
 * averageTKE
 *------------------------------------------------------------------------------
 * Returns the <TKE> given an "effective" incident neutron energy e
 ******************************************************************************/
double Yields::averageTKE(double e) {

    int n_poly = 2;
    double sum = TKE_fit_param[0];
    double e1 = e;
    for (int i=1; i<n_poly; i++) {
    	sum += TKE_fit_param[i]*e1;
    	e1 *= e;
    }
    return (sum);
}

/*******************************************************************************
 * sampleTKE
 *------------------------------------------------------------------------------
 * Samples a TKE value from a Gaussian with mean avTKE[a] and width sTKE[a] for
 * a given mass A.
 ******************************************************************************/
double Yields::sampleTKE(int A) {
   
    int a;
    if (A>=a_cn/2+a_cn%2) {
    	a = A;
    } else {
	   a = a_cn - A;
    }
    return(rescaleAverageTKE*avTKE[a]+sTKE[a]*sampleGaussian());
}

/*******************************************************************************
 * get_avTKE
 *------------------------------------------------------------------------------
 * Returns the <TKE> for a particular mass A
 ******************************************************************************/
double Yields::get_avTKE(int A) {
   
    int a;
    if (A>=a_cn/2+a_cn%2) {
    	a = A;
    } else {
	   a = a_cn - A;
    }
    return(rescaleAverageTKE*avTKE[a]);
}

/*******************************************************************************
 * sampleGaussian
 *------------------------------------------------------------------------------
 * Samples from a Gaussian with mean 1 and width 1
 ******************************************************************************/
double Yields::sampleGaussian(void) {

    double rand1;
    double rand2;
  
    double rsq;
  
    do{
    	rand1 = 2.*rng_cgm() - 1.;
    	rand2 = 2.*rng_cgm() - 1.;
    	rsq = rand1*rand1 + rand2*rand2;
    } while (rsq == 0. || rsq >=1.);
  
    rsq = sqrt(-(double)2.*log(rsq)/rsq);
    rand2 *= rsq;
    return (rand1*rsq);
}
