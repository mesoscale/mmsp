
//------------------------Initialization of global variables/parameters-----------------------//

double R = 8.314 ;
double T = 1023 ;
double T_orig = 1023;
double cr = 0.0;        //Cooling rate (K/s)
double L_orig = 0.139 ;  //Phase field mobility
double c_tot = 0.1 ;
double sigma = 0.150;
double V0 = 1.0e-02/60.0 ;
double epsi_sq = 1.0e-07 ; 
double L = (sigma*V0)/(epsi_sq*c_tot*R*T);
double G_normalize = R*T ;  //Scaling factor for terms having J/mol units
double lold = 9.677*pow(10, -8) ;  //Grid size (m)
double lc = 9.677*pow(10,-8) ;    //Characteristic length scale, or capillary length
double dx_nd = lold/lc ;          // Non-dimensional grid size
double tc = 4.4 ;               //Characteristic time scale (s)
double dt_old = 0.44 ; 			// Dimensional dt (s)
double dt = dt_old/tc;        //Non-dimensional dt
double Dalal = 0.001 ;        // The next four lines are the components of the Onsager mobility matrix
double Dalv = -0.0008 ;
double Dvv = 0.001 ;
double Dval = -0.0008 ;
double Vm = 0.00001;         //Molar volume (m^3/mol)
double w_norm = (R*T)/Vm ;       //Scaling factor for the double well height
double epsi_norm = (R*T*pow(lc,2))/Vm;    //Scaling factor for the gradient energy coeffcient
int variants = 20 ;          // Number of variants : input parameter needed for defining a mmsp-grid
const int dim = 2;     //Spatial dimensions
const int nx = 24;
const int ny = 24;
double W_prefac = 6.095 ;   //Overall scaling for the double well depth              
double W_Al = 0.1 ;       //Relative magnitudes of double well depths for the three components
double W_V =  0.1 ; 
double W_Ti = 0.15 ;
double scaling = 1.0/(2*3.14*3.14) ;        //Scaling factor, needed for Fourier Transformation integrals. 
int steps = 10000 ;    //Number of simulation steps.
double kappa1 = 0.001 ; 
double alpha = 1.0 ;
double gammanuc = 0.5 ;
double kappa_c = 0.076 ;

//------------------------Definition of global containers/arrays-----------------------//

const int node_total = nx*ny; 
double G_Alpha[node_total], G_Beta[node_total], W[node_total] ; 
double hphi[node_total][12], hphiprime[node_total][12], hphidoubleprime[node_total][12];
double gphi[node_total][12], gphiprime[node_total][12], gphidoubleprime[node_total][12];
double strain_energy[node_total][3];
double intenergy[13][nx][ny];
double dGAlpha_dAl,	dGBeta_dAl, dGAlpha_dV, dGBeta_dV, del_dGAlpha_dAl, del_dGBeta_dAl, del_dGAlpha_dV, del_dGBeta_dV, delsq_dGAlpha_dAl, delsq_dGAlpha_dV, delsq_dGBeta_dAl, delsq_dGBeta_dV ;
double G_Al_alpha, G_Ti_alpha, G_V_alpha, G_Al_beta, G_Ti_beta, G_V_beta ; 
double epsi[12][3][3] ;


//-----------------------Controlling the number and type of output files----------------//

// If a particular file is not needed, change the value to 0
bool output_pfsq = 1 ;
bool output_cal = 1 ;
bool output_cv = 1 ;
bool output_chemnuc = 1 ;
bool output_strainintnuc = 1 ;

//---------------------Definition of functions-------------------------//
void thermo_auxillary_terms(MMSP::vector<store_type> gradient, MMSP::vector<store_type> gradientsq, double c_Al, double c_V) ;
void customoutput(int t) ;
void calculate_strains() ;
void fft(CArray& x) ;
void ifft(CArray& x) ;
void intstrain();
double nodesum() ;


//---------------------Definition of grids and grid variables-----------------------//
MMSP::grid<dim, store_type> grid(variants, 0, nx, 0, ny) ;
MMSP::grid<dim, store_type> nuc_grid(variants, 0, nx, 0, ny) ;
MMSP::grid<dim, store_type> gradsqcal_grid(variants, 0, nx, 0, ny) ;
MMSP::grid<dim, store_type> gradsqcv_grid(variants, 0, nx, 0, ny) ;
MMSP::grid<dim, store_type> intenergies(variants, 0, nx, 0, ny) ;
MMSP::grid<dim, store_type> selfenergies(variants, 0, nx, 0, ny) ;

const double Lx = g1(grid,0) - g0(grid,0) ;
const double Ly = g1(grid,1) - g0(grid,1) ;
const double dx = MMSP::dx(grid, 0) ;
const double dy = MMSP::dx(grid, 1) ;
double fs = nx/Lx;
double k[nx], fk[nx] ;

//-------------------Definition of user-defined classes---------------------------------//
class sfts
{
	
	public:
	double e[6] ;
	
	sfts()
	{
		e[0] = 0.0 ;
		e[1] = 0.0 ;
		e[2] = 0.0 ;
		e[3] = 0.0 ;
		e[4] = 0.0 ;
		e[5] = 0.0 ;
		
	}
};

sfts eigen_alpha[12];



//----------Including the various modules-----------//

#include "Array.h"
#include "/home/arun/Documents/mmsp-arun/fftw++-2.05/fftw++.h"

using namespace utils;
using namespace Array;
using namespace fftwpp;



#include "thermo_modules.hpp"
#include "strain_modules.hpp"



//---------User defined functions--------------//

double nodesum(MMSP::grid<dim, store_type> grid, MMSP::vector<int> s)
{
	double phisum = 0 ;
	for(int h = 0 ; h < length(grid(s)) ; h++)
	{
		int hindex = MMSP::index(grid(s),h) ;
		if((hindex <= 12))
		{
			if(grid(s)[hindex] > 0.0001)
			{
				phisum += grid(s)[hindex] ; 
			}	
		}
	}
	return phisum ;
}



void initialize_epsi()
{
	
epsi[0][0][0] = 0.0806 ;
epsi[0][0][1] = 0.0301 ;
epsi[0][0][2] = -0.0642 ;
epsi[0][1][0] = 0.0301 ;
epsi[0][1][1] = 0.0806 ;
epsi[0][1][2] = -0.0277 ;
epsi[0][2][0] = -0.0642 ;
epsi[0][2][1] = -0.0277 ;
epsi[0][2][2] = 0.3834 ;

epsi[1][0][0] = 0.3445 ;
epsi[1][0][1] = -0.1183 ;
epsi[1][0][2] = -0.0486 ;
epsi[1][1][0] = -0.1183 ;
epsi[1][1][1] = 0.1056 ;
epsi[1][1][2] = 0.0011 ;
epsi[1][2][0] = -0.0486 ;
epsi[1][2][1] = 0.0011 ;
epsi[1][2][2] = 0.0999 ;

epsi[2][0][0] = 0.0895 ;
epsi[2][0][1] = -0.0300 ;
epsi[2][0][2] = -0.0635 ;
epsi[2][1][0] = -0.0300 ;
epsi[2][1][1] = 0.0759 ;
epsi[2][1][2] = 0.0214 ;
epsi[2][2][0] = -0.0635 ;
epsi[2][2][1] = 0.0214 ;
epsi[2][2][2] = 0.3847 ;

epsi[3][0][0] = 0.0895 ;
epsi[3][0][1] = 0.0300 ;
epsi[3][0][2] = -0.0635 ;
epsi[3][1][0] = 0.0300 ;
epsi[3][1][1] = 0.0759 ;
epsi[3][1][2] = -0.0214 ;
epsi[3][2][0] = -0.0635 ;
epsi[3][2][1] = -0.0214 ;
epsi[3][2][2] = 0.3847 ;

epsi[4][0][0] = 0.0989 ;
epsi[4][0][1] = -0.0021 ;
epsi[4][0][2] = -0.0226 ;
epsi[4][1][0] = -0.0021 ;
epsi[4][1][1] = 0.1575 ;
epsi[4][1][2] = 0.1592 ;
epsi[4][2][0] = -0.0227 ;
epsi[4][2][1] = 0.1592 ;
epsi[4][2][2] = 0.2936 ;

epsi[5][0][0] = 0.0806 ;
epsi[5][0][1] = -0.0301 ;
epsi[5][0][2] = -0.0642 ;
epsi[5][1][0] = -0.0301 ;
epsi[5][1][1] = 0.0860 ;
epsi[5][1][2] = 0.0277 ;
epsi[5][2][0] = -0.0642 ;
epsi[5][2][1] = 0.0277 ;
epsi[5][2][2] = 0.3834 ;

epsi[6][0][0] = 0.3240 ;
epsi[6][0][1] = -0.1376 ;
epsi[6][0][2] = -0.0411 ;
epsi[6][1][0] = -0.1376 ;
epsi[6][1][1] = 0.1322 ;
epsi[6][1][2] = -0.0016 ;
epsi[6][2][0] = -0.0411 ;
epsi[6][2][1] = -0.0016 ;
epsi[6][2][2] = 0.0938 ;

epsi[7][0][0] = 0.0938 ;
epsi[7][0][1] = 0.0016 ;
epsi[7][0][2] = -0.0411 ;
epsi[7][1][0] = 0.0016 ;
epsi[7][1][1] = 0.1322 ;
epsi[7][1][2] = 0.1376 ;
epsi[7][2][0] = -0.0411 ;
epsi[7][2][1] = 0.1376 ;
epsi[7][2][2] = 0.3240 ;

epsi[8][0][0] = 0.0758 ;
epsi[8][0][1] = -0.0259 ;
epsi[8][0][2] = -0.0278 ;
epsi[8][1][0] = -0.0259 ;
epsi[8][1][1] = 0.0771 ;
epsi[8][1][2] = 0.0101 ;
epsi[8][2][0] = -0.0278 ;
epsi[8][2][1] = 0.0101 ;
epsi[8][2][2] = 0.3971 ;

epsi[9][0][0] = 0.3456 ;
epsi[9][0][1] = -0.1251 ;
epsi[9][0][2] = -0.0102 ;
epsi[9][1][0] = -0.1251 ;
epsi[9][1][1] = 0.1123 ;
epsi[9][1][2] = -0.0155 ;
epsi[9][2][0] = -0.0102 ;
epsi[9][2][1] = -0.0155 ;
epsi[9][2][2] = 0.0121 ;

epsi[10][0][0] = 0.0758 ;
epsi[10][0][1] = 0.0259 ;
epsi[10][0][2] = -0.0278 ;
epsi[10][1][0] = 0.0259 ;
epsi[10][1][1] = 0.0771 ;
epsi[10][1][2] = -0.0101 ;
epsi[10][2][0] = -0.0278 ;
epsi[10][2][1] = -0.0101 ;
epsi[10][2][2] = 0.3971 ;

epsi[11][0][0] = 0.0863 ;
epsi[11][0][1] = 0.0177 ;
epsi[11][0][2] = -0.0254 ;
epsi[11][1][0] = 0.0177 ;
epsi[11][1][1] = 0.0819 ;
epsi[11][1][2] = 0.0729 ;
epsi[11][2][0] = -0.0254 ;
epsi[11][2][1] = 0.0729 ;
epsi[11][2][2] = 0.3818 ;

}


void initialize_alpha_eigen()
{

eigen_alpha[0].e[0] = -0.083  ;
eigen_alpha[0].e[5] = 0.0095  ;
eigen_alpha[0].e[4] = 0.0  ;
eigen_alpha[0].e[1] = 0.123  ;
eigen_alpha[0].e[3] = 0.0  ;
eigen_alpha[0].e[2] = 0.035  ;

eigen_alpha[1].e[0] = -0.083  ;
eigen_alpha[1].e[5] = 0.0  ;
eigen_alpha[1].e[4] = 0.0095  ;
eigen_alpha[1].e[1] = 0.035  ;
eigen_alpha[1].e[3] = 0.0  ;
eigen_alpha[1].e[2] = 0.123  ;

eigen_alpha[2].e[0] = 0.079  ;
eigen_alpha[2].e[5] = -0.0359  ;
eigen_alpha[2].e[4] = -0.0264  ;
eigen_alpha[2].e[1] = 0.0047  ;
eigen_alpha[2].e[3] = 0.0810  ;
eigen_alpha[2].e[2] = -0.0087  ;

eigen_alpha[3].e[0] = -0.079  ;
eigen_alpha[3].e[5] = 0.0359  ;
eigen_alpha[3].e[4] = 0.0264  ;
eigen_alpha[3].e[1] = 0.0047  ;
eigen_alpha[3].e[3] = 0.0810  ;
eigen_alpha[3].e[2] = -0.0087  ;

eigen_alpha[4].e[0] = 0.079  ;
eigen_alpha[4].e[5] = -0.0359  ;
eigen_alpha[4].e[4] = 0.0264  ;
eigen_alpha[4].e[1] = 0.0047  ;
eigen_alpha[4].e[3] = -0.0810  ;
eigen_alpha[4].e[2] = -0.0087  ;

eigen_alpha[5].e[0] = 0.079  ;
eigen_alpha[5].e[5] = 0.0359  ;
eigen_alpha[5].e[4] = -0.0264  ;
eigen_alpha[5].e[1] = 0.0047  ;
eigen_alpha[5].e[3] = -0.0810  ;
eigen_alpha[5].e[2] = -0.0087  ;

eigen_alpha[6].e[0] = -0.083  ;
eigen_alpha[6].e[5] = -0.0095  ;
eigen_alpha[6].e[4] = 0.0  ;
eigen_alpha[6].e[1] = 0.123  ;
eigen_alpha[6].e[3] = 0.0  ;
eigen_alpha[6].e[2] = 0.035  ;

eigen_alpha[7].e[0] = -0.083  ;
eigen_alpha[7].e[5] = 0.0  ;
eigen_alpha[7].e[4] = -0.0095  ;
eigen_alpha[7].e[1] = 0.035  ;
eigen_alpha[7].e[3] = 0.0  ;
eigen_alpha[7].e[2] = 0.123  ;

eigen_alpha[8].e[0] = 0.079  ;
eigen_alpha[8].e[5] = -0.0264  ;
eigen_alpha[8].e[4] = -0.0359  ;
eigen_alpha[8].e[1] = -0.0087  ;
eigen_alpha[8].e[3] = 0.0810  ;
eigen_alpha[8].e[2] = 0.0047  ;

eigen_alpha[9].e[0] = 0.079  ;
eigen_alpha[9].e[5] = 0.0264  ;
eigen_alpha[9].e[4] = 0.0359  ;
eigen_alpha[9].e[1] = -0.0087  ;
eigen_alpha[9].e[3] = 0.0810  ;
eigen_alpha[9].e[2] = 0.0047  ;

eigen_alpha[10].e[0] = 0.079  ;
eigen_alpha[10].e[5] = -0.0264  ;
eigen_alpha[10].e[4] = 0.0359  ;
eigen_alpha[10].e[1] = -0.0087  ;
eigen_alpha[10].e[3] = -0.0810  ;
eigen_alpha[10].e[2] = 0.0047  ;

eigen_alpha[11].e[0] = 0.079  ; 
eigen_alpha[11].e[5] = 0.0264  ; 
eigen_alpha[11].e[4] = -0.0359  ; 
eigen_alpha[11].e[1] = -0.0087  ;
eigen_alpha[11].e[3] = -0.0810  ; 
eigen_alpha[11].e[2] = 0.0047  ; 
}




