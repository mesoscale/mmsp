#include <stdio.h>
#include <math.h>
#include"MMSP.hpp"
using namespace MMSP;

typedef float phi_type; 											
typedef MMSP::sparse<phi_type> store_type; 	

double hphi(double phi)
{
    return pow(phi,3)*(10-15*phi+6*pow(phi,2)) ;
}

double hphiprime(double phi)
{
    return 30*pow(phi,2)*pow((phi-1),2) ;
}
int main(int argc, char* argv[])
{
	Init(argc, argv);

	double R = 8.314 ;
    double T = 1200.0 ;
    double dt = 0.01 ;        
    double Vm = 0.00001;   
    int variants = 20 ;          
    const int dim = 2;  
    const int nx = 100;
    const int ny = 100;
    long steps = 5e+8 ; 
	int iterations = 200;
	const int node_total = nx*ny; 
    double epsisq = 1 ;
    double W[node_total] ;
    double diffusionCoefficient = 1.0;


	grid<1,scalar<double> > oldGrid(1,0,nx);
    grid<2, store_type> newGrid(variants, 0, nx, 0, ny) ;
    grid<2, store_type > update(newGrid);
    
    const double Lx = g1(newGrid,0) - g0(newGrid,0) ;
    const double Ly = g1(newGrid,1) - g0(newGrid,1) ;
    const double dx = MMSP::dx(newGrid, 0) ;
    const double dy = MMSP::dx(newGrid, 1) ;
    
    MMSP::b0(newGrid,0) = MMSP::Dirichlet;
	MMSP::b1(newGrid,0) = MMSP::Dirichlet;

	MMSP::b0(newGrid,1) = MMSP::Dirichlet;
	MMSP::b1(newGrid,1) = MMSP::Dirichlet;
    
    MMSP::b0(update,0) = MMSP::Dirichlet;
    MMSP::b1(update,0) = MMSP::Dirichlet;

    MMSP::b0(update,1) = MMSP::Dirichlet;
    MMSP::b1(update,1) = MMSP::Dirichlet;
    double rstar = 5.0 ;
    double r0 = rstar ;
    double x0 = 0.5*Lx ;
    double y0 = 0.5*Ly ;
    for(int n = 0; n < nodes(newGrid); n++)
	{
		MMSP::vector<int> x = position(newGrid, n);        
        store_type newGrain;
        double x1 = x[0]*dx ;
        double y1 = x[1]*dy ;
        double r = sqrt(pow((x1-x0),2) + pow((y1-y0),2)) ;
        set(newGrain, 1) = 0.5*(1-tanh(r-r0)/sqrt(2)) ;
		newGrid(n) = newGrain;
        update(n) = newGrain ;
	}
    
    
    double F = 0.0 ;
    std::string file_name = "F.txt" ;
    std::ofstream file1;
    file1.open(file_name.c_str());
    
	for (int k=0; k<iterations; k++) 
    {
        if(k%100==0)   
		{
            F = 0 ;
            std::string file_name = "output_" + std::to_string(k) + ".txt" ;
            std::ofstream myfile;
            myfile.open(file_name.c_str());
            for(int n = 0 ; n < nodes(newGrid) ; n++)    
            {
                MMSP::vector<int> s = position(newGrid, n); 
                MMSP::vector<store_type> gradient = grad(update, s) ;
                double maggrad = gradient[0][1] + gradient[1][1] ;  
                F += 0.5*epsisq*pow(maggrad,2) + pow(update(n)[1],2) * pow((1 - update(n)[1]),2) -1.0/(15.0*sqrt(2.0))*hphi(update(n)[1]) ; 
                std::cout<<hphi(update(n)[1])<<" "<<F<<std::endl;
                myfile<<s[0]<<","<<s[1]<<","<<newGrid(s)[1]<<"\n" ; 
            }
            file1<<F<<"\n";
            myfile.close();
        }
        
		for(int n = 0 ; n < nodes(update) ; n++)
		{
			MMSP::vector<int> x = position(update, n); 
			MMSP::vector<store_type> gradsq = gradientsq(update, x) ;
			MMSP::sparse<phi_type> dFdp;
            store_type dpdt;
            phi_type lap_aniso = (gradsq[0][1]*epsisq + gradsq[1][1]*epsisq)  ;
            double interaction_energy = update(n)[1]*pow((1-update(n)[1]), 2) - pow(update(n)[1],2)*(1-update(n)[1]) ;
            set(dFdp, 1) = -(1)*hphiprime(update(n)[1])  + interaction_energy - lap_aniso ; 
            set(dpdt, 1) = -0.1* (dFdp[1]);
            set(update(n), 1) = update(n)[1] + dt * dpdt[1];							
		}
        swap(update, newGrid);
	}
    file1.close();
	Finalize();
	return 0;
}

