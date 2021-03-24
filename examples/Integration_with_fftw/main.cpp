#include<iomanip>
#include<vector>
#include<cmath>
#include<ctime>
#include<cstdio>
#include<cstdlib>
#include<cassert>
#include<iostream>
#include<fstream>
#include<sstream>
#include<string>
#include <time.h>
#include<deque>
#include <stdio.h>
#include <stdlib.h>
#include<vector>
#include <complex>
#include <valarray>
 
#include"MMSP.hpp"


using namespace std;

typedef float phi_type; 											
typedef MMSP::sparse<phi_type> store_type; 	
typedef std::complex<double> Complex;
typedef std::valarray<Complex> CArray;

const double PI = 3.141592653589793238460;

#include "definitions_v2.hpp"

const int variant = 1 ;


//--------Initializing all the grids--------//

void initialize()
{
	
	G_Al_alpha = (Al_alpha_a_3 + Al_alpha_b_3*T+ Al_alpha_c_3*T*log(T) + Al_alpha_d_3*T*T) ;
	G_Al_beta = (Al_beta_a_3 + Al_beta_b_3*T+ Al_beta_c_3*T*log(T) + Al_beta_d_3*T*T) ;
	
	G_Ti_alpha = (Ti_alpha_a_2 + Ti_alpha_b_2*T+ Ti_alpha_c_2*T*log(T) + Ti_alpha_d_2*T*T) ;
	G_Ti_beta = (Ti_beta_a_1 + Ti_beta_b_1*T+ Ti_beta_c_1*T*log(T) + Ti_beta_d_1*T*T) ;
	
	G_V_alpha = (V_alpha_a_2 + V_alpha_b_2*T+ V_alpha_c_2*T*log(T) + V_alpha_d_2*T*T) ; //(21.96 - T/50)
	G_V_beta = (V_beta_a_2 + V_beta_b_2*T+ V_beta_c_2*T*log(T) + V_beta_d_2*T*T) ;	

	for(int n = 0; n < nodes(grid); n++)
	{
		MMSP::vector<int> x = position(grid, n);
        
        store_type newGrain;
        store_type nuc_newGrain ;
        
        set(nuc_newGrain, 1) = 0.0 ;
        set(nuc_newGrain, 3) = 0.0 ;
        nuc_grid(n) = nuc_newGrain ;
        
        set(intenergies(n), 1) = 0.0 ;
        set(intenergies(n), 3) = 0.0 ;
        set(selfenergies(n), 1) = 0.0 ;
        set(selfenergies(n), 3) = 0.0 ;
        
        set(newGrain, 1) = 0.0;
        set(newGrain, 2) = 0.0;
        set(newGrain, 3) = 0.0;
        set(newGrain, 4) = 0.0;
        set(newGrain, 5) = 0.0;
        set(newGrain, 6) = 0.0;
        set(newGrain, 7) = 0.0;
        set(newGrain, 8) = 0.0;
        set(newGrain, 9) = 0.0;
        set(newGrain, 10) = 0.0;
        set(newGrain, 11) = 0.0;
        set(newGrain, 12) = 0.0;
        
        if((pow((x[0]*MMSP::dx(grid,0)-0.2*Lx),2) + pow((x[1]*MMSP::dx(grid,1)-0.8*Ly),2)) < 100*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
			set(newGrain, 1) = 1.0;	
		}
        
        else if((pow((x[0]*MMSP::dx(grid,0)-0.4*Lx),2) + pow((x[1]*MMSP::dx(grid,1)-0.8*Ly),2)) < 100*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
			set(newGrain, 2) = 1.0;	
		}
        
        else if((pow((x[0]*MMSP::dx(grid,0)-0.6*Lx),2) + pow((x[1]*MMSP::dx(grid,1)-0.8*Ly),2)) < 100*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
			set(newGrain, 3) = 1.0;	
		}
        
        else if((pow((x[0]*MMSP::dx(grid,0)-0.8*Lx),2) + pow((x[1]*MMSP::dx(grid,1)-0.8*Ly),2)) < 100*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
			set(newGrain, 4) = 1.0;	
		}
        
        else if((pow((x[0]*MMSP::dx(grid,0)-0.2*Lx),2) + pow((x[1]*MMSP::dx(grid,1)-0.5*Ly),2)) < 100*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
			set(newGrain, 5) = 1.0;	
		}
        
        else if((pow((x[0]*MMSP::dx(grid,0)-0.4*Lx),2) + pow((x[1]*MMSP::dx(grid,1)-0.5*Ly),2)) < 100*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
			set(newGrain, 6) = 1.0;	
		}
        
        else if((pow((x[0]*MMSP::dx(grid,0)-0.6*Lx),2) + pow((x[1]*MMSP::dx(grid,1)-0.5*Ly),2)) < 100*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
			set(newGrain, 7) = 1.0;	
		}
        
        else if((pow((x[0]*MMSP::dx(grid,0)-0.8*Lx),2) + pow((x[1]*MMSP::dx(grid,1)-0.5*Ly),2)) < 100*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
			set(newGrain, 8) = 1.0;	
		}
        
        else if((pow((x[0]*MMSP::dx(grid,0)-0.2*Lx),2) + pow((x[1]*MMSP::dx(grid,1)-0.2*Ly),2)) < 100*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
			set(newGrain, 9) = 1.0;	
		}
        
        else if((pow((x[0]*MMSP::dx(grid,0)-0.4*Lx),2) + pow((x[1]*MMSP::dx(grid,1)-0.2*Ly),2)) < 100*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
			set(newGrain, 10) = 1.0;	
		}
        
        else if((pow((x[0]*MMSP::dx(grid,0)-0.6*Lx),2) + pow((x[1]*MMSP::dx(grid,1)-0.2*Ly),2)) < 100*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
			set(newGrain, 11) = 1.0;	
		}
        
        else if((pow((x[0]*MMSP::dx(grid,0)-0.8*Lx),2) + pow((x[1]*MMSP::dx(grid,1)-0.2*Ly),2)) < 100*MMSP::dx(grid,0)*MMSP::dx(grid,0))
		{
			set(newGrain, 12) = 1.0;	
		}

		
		grid(n) = newGrain;		
		
	}

}



int main()
{	
    
	initialize();
	initialize_epsi();
	initialize_alpha_eigen() ;
	define_c_sigma() ;
	
    fftw::maxthreads=get_max_threads();
     
    unsigned int nyp=ny/2+1;
     
    	

size_t align=sizeof(Complex);

array2<double> f1(nx,ny,align);
array2<Complex> F1(nx,nyp,align);
rcfft2d Forward1(nx,ny,f1,F1);

array2<double> f2(nx,ny,align);
array2<Complex> F2(nx,nyp,align);
rcfft2d Forward2(nx,ny,f2,F2);

array2<double> f3(nx,ny,align);
array2<Complex> F3(nx,nyp,align);
rcfft2d Forward3(nx,ny,f3,F3);

array2<double> f4(nx,ny,align);
array2<Complex> F4(nx,nyp,align);
rcfft2d Forward4(nx,ny,f4,F4);

array2<double> f5(nx,ny,align);
array2<Complex> F5(nx,nyp,align);
rcfft2d Forward5(nx,ny,f5,F5);

array2<double> f6(nx,ny,align);
array2<Complex> F6(nx,nyp,align);
rcfft2d Forward6(nx,ny,f6,F6);

array2<double> f7(nx,ny,align);
array2<Complex> F7(nx,nyp,align);
rcfft2d Forward7(nx,ny,f7,F7);

array2<double> f8(nx,ny,align);
array2<Complex> F8(nx,nyp,align);
rcfft2d Forward8(nx,ny,f8,F8);

array2<double> f9(nx,ny,align);
array2<Complex> F9(nx,nyp,align);
rcfft2d Forward9(nx,ny,f9,F9);

array2<double> f10(nx,ny,align);
array2<Complex> F10(nx,nyp,align);
rcfft2d Forward10(nx,ny,f10,F10);

array2<double> f11(nx,ny,align);
array2<Complex> F11(nx,nyp,align);
rcfft2d Forward11(nx,ny,f11,F11);

array2<double> f12(nx,ny,align);
array2<Complex> F12(nx,nyp,align);
rcfft2d Forward12(nx,ny,f12,F12);

array2<Complex> dfdstr1(nx,nyp,align);
array2<double> dfdstr_real1(nx,ny,align);
crfft2d Backward1(nx,ny,dfdstr1,dfdstr_real1);

array2<Complex> dfdstr2(nx,nyp,align);
array2<double> dfdstr_real2(nx,ny,align);
crfft2d Backward2(nx,ny,dfdstr2,dfdstr_real2);

array2<Complex> dfdstr3(nx,nyp,align);
array2<double> dfdstr_real3(nx,ny,align);
crfft2d Backward3(nx,ny,dfdstr3,dfdstr_real3);

array2<Complex> dfdstr4(nx,nyp,align);
array2<double> dfdstr_real4(nx,ny,align);
crfft2d Backward4(nx,ny,dfdstr4,dfdstr_real4);

array2<Complex> dfdstr5(nx,nyp,align);
array2<double> dfdstr_real5(nx,ny,align);
crfft2d Backward5(nx,ny,dfdstr5,dfdstr_real5);

array2<Complex> dfdstr6(nx,nyp,align);
array2<double> dfdstr_real6(nx,ny,align);
crfft2d Backward6(nx,ny,dfdstr6,dfdstr_real6);

array2<Complex> dfdstr7(nx,nyp,align);
array2<double> dfdstr_real7(nx,ny,align);
crfft2d Backward7(nx,ny,dfdstr7,dfdstr_real7);

array2<Complex> dfdstr8(nx,nyp,align);
array2<double> dfdstr_real8(nx,ny,align);
crfft2d Backward8(nx,ny,dfdstr8,dfdstr_real8);

array2<Complex> dfdstr9(nx,nyp,align);
array2<double> dfdstr_real9(nx,ny,align);
crfft2d Backward9(nx,ny,dfdstr9,dfdstr_real9);

array2<Complex> dfdstr10(nx,nyp,align);
array2<double> dfdstr_real10(nx,ny,align);
crfft2d Backward10(nx,ny,dfdstr10,dfdstr_real10);

array2<Complex> dfdstr11(nx,nyp,align);
array2<double> dfdstr_real11(nx,ny,align);
crfft2d Backward11(nx,ny,dfdstr11,dfdstr_real11);

array2<Complex> dfdstr12(nx,nyp,align);
array2<double> dfdstr_real12(nx,ny,align);
crfft2d Backward12(nx,ny,dfdstr12,dfdstr_real12);

array2<Complex> elint1(nx,nyp,align);
array2<double> elint_real1(nx,ny,align);

array2<Complex> elint3(nx,nyp,align);
array2<double> elint_real3(nx,ny,align);

crfft2d Backward13(nx,ny,elint1,elint_real1);
crfft2d Backward14(nx,ny,elint3,elint_real3);

array2<Complex> self1(nx,nyp,align);
array2<double> self_real1(nx,ny,align);
array2<Complex> self3(nx,nyp,align);
array2<double> self_real3(nx,ny,align);

crfft2d Backward15(nx, ny, self1, self_real1);
crfft2d Backward16(nx, ny, self3, self_real3);

array2<Complex> fstr(nx,nyp,align);
array2<double> fstr_real(nx,ny,align);

crfft2d Backward17(nx,ny,fstr,fstr_real);

	
	//-------Setting the boundary conditions for various grids-------//
	MMSP::b0(grid,0) = MMSP::Dirichlet;
	MMSP::b1(grid,0) = MMSP::Dirichlet;

	MMSP::b0(grid,1) = MMSP::Dirichlet;
	MMSP::b1(grid,1) = MMSP::Dirichlet;

	MMSP::b0(grid,2) = MMSP::Dirichlet;
	MMSP::b1(grid,2) = MMSP::Dirichlet;

	
	//------Defining the lattice in the Fourier space------//
	for(int i = 0 ; i < nx ; i++)
    {
        fk[i] = (2.0*3.14*i)/(double)nx  ;
    }
  
	
	srand(time(NULL));

    for(int t = 0 ; t<steps ; t++)
	{
	
        //cout<<t<<endl;

		if((T-dt*tc*cr) > 100.0)    //Implementing the variation of Temperature based on the given cooling rate
		{ T= T- dt*tc*cr ; }
		
		
		if(cr!=0.0)
		{
			G_Al_alpha = (Al_alpha_a_3 + Al_alpha_b_3*T+ Al_alpha_c_3*T*log(T) + Al_alpha_d_3*T*T) ;
			G_Al_beta = (Al_beta_a_3 + Al_beta_b_3*T+ Al_beta_c_3*T*log(T) + Al_beta_d_3*T*T) ;
	
			G_Ti_alpha = (Ti_alpha_a_2 + Ti_alpha_b_2*T+ Ti_alpha_c_2*T*log(T) + Ti_alpha_d_2*T*T) ;
			G_Ti_beta = (Ti_beta_a_1 + Ti_beta_b_1*T+ Ti_beta_c_1*T*log(T) + Ti_beta_d_1*T*T) ;
	
			G_V_alpha = (V_alpha_a_2 + V_alpha_b_2*T+ V_alpha_c_2*T*log(T) + V_alpha_d_2*T*T) ; //(21.96 - T/50)
			G_V_beta = (V_beta_a_2 + V_beta_b_2*T+ V_beta_c_2*T*log(T) + V_beta_d_2*T*T) ;
		}

		MMSP::grid<dim, store_type > update(grid);
	
		MMSP::b0(update,0) = MMSP::Dirichlet;
		MMSP::b1(update,0) = MMSP::Dirichlet;

		MMSP::b0(update,1) = MMSP::Dirichlet;
		MMSP::b1(update,1) = MMSP::Dirichlet;

		MMSP::b0(update,2) = MMSP::Dirichlet;
		MMSP::b1(update,2) = MMSP::Dirichlet;
	
		
        
        
        
    for(int n = 0 ; n < nodes(grid) ; n++)
    {
        MMSP::vector<int> x = position(grid, n) ;
        f1(x[0], x[1])=grid(n)[1] ;   
        f2(x[0], x[1])=grid(n)[2] ;   
        f3(x[0], x[1])=grid(n)[3] ;   
        f4(x[0], x[1])=grid(n)[4] ;   
        f5(x[0], x[1])=grid(n)[5] ;   
        f6(x[0], x[1])=grid(n)[6] ;   
        f7(x[0], x[1])=grid(n)[7] ;   
        f8(x[0], x[1])=grid(n)[8] ;   
        f9(x[0], x[1])=grid(n)[9] ;   
        f10(x[0], x[1])=grid(n)[10] ;   
        f11(x[0], x[1])=grid(n)[11] ;   
        f12(x[0], x[1])=grid(n)[12] ;   
    }
    
    Forward1.fft0(f1,F1);  
    Forward2.fft0(f2,F2);  
    Forward3.fft0(f3,F3);  
    Forward4.fft0(f4,F4);  
    Forward5.fft0(f5,F5);  
    Forward6.fft0(f6,F6);  
    Forward7.fft0(f7,F7);  
    Forward8.fft0(f8,F8);  
    Forward9.fft0(f9,F9);  
    Forward10.fft0(f10,F10);  
    Forward11.fft0(f11,F11);  
    Forward12.fft0(f12,F12);
    
    for(int i = 0 ; i < nx ; i++)
    {
        for(int j = 0 ; j < nyp ; j++)
        {
                double k1 = fk[i] ;
                double k2 = fk[j] ;
		
                double modk = (sqrt(k1*k1 + k2*k2)) ;
                if(modk==0.0) 
                {
                    elint1(i,j) = 0.0   ;
                    elint3(i,j) = 0.0   ;
                    self1(i,j) = 0.0   ;
                    self3(i,j) = 0.0   ;
                }
                else 
                {
                    self1(i,j) = Bpq(k1, k2, modk, 1, 1) ;
                    self3(i,j) = Bpq(k1, k2, modk, 3, 3) ;
                    elint1(i,j) = Bpq(k1, k2, modk, 1, 1)*F1(i,j) + Bpq(k1, k2, modk, 2, 1)*F2(i,j) + Bpq(k1, k2, modk, 3, 1)*F3(i,j) + Bpq(k1, k2, modk, 4, 1)*F4(i,j) 
                        + Bpq(k1, k2, modk, 5, 1)*F5(i,j) + Bpq(k1, k2, modk, 6, 1)*F6(i,j) + Bpq(k1, k2, modk, 7, 1)*F7(i,j) + Bpq(k1, k2, modk, 8, 1)*F8(i,j) + Bpq(k1, k2, modk, 9, 1)*F9(i,j)
                        + Bpq(k1, k2, modk, 10, 1)*F10(i,j) + Bpq(k1, k2, modk, 11, 1)*F11(i,j) + Bpq(k1, k2, modk, 12, 1)*F12(i,j) ;
                        
                    elint3(i,j) = Bpq(k1, k2, modk, 1, 3)*F1(i,j) + Bpq(k1, k2, modk, 2, 3)*F2(i,j) + Bpq(k1, k2, modk, 3, 3)*F3(i,j) + Bpq(k1, k2, modk, 4, 3)*F4(i,j) 
                        + Bpq(k1, k2, modk, 5, 3)*F5(i,j) + Bpq(k1, k2, modk, 6, 3)*F6(i,j) + Bpq(k1, k2, modk, 7, 3)*F7(i,j) + Bpq(k1, k2, modk, 8, 3)*F8(i,j) + Bpq(k1, k2, modk, 9, 3)*F9(i,j)
                        + Bpq(k1, k2, modk, 10, 3)*F10(i,j) + Bpq(k1, k2, modk, 11, 3)*F11(i,j) + Bpq(k1, k2, modk, 12, 3)*F12(i,j) ;
                }
        }
    }
            
    Backward13.fft0Normalized(elint1, elint_real1); 
    Backward14.fft0Normalized(elint3, elint_real3); 
    Backward15.fft0Normalized(self1, self_real1); 
    Backward16.fft0Normalized(self3, self_real3); 
    
    
    for(int n = 0 ; n < nodes(grid) ; n++)
    {
        MMSP::vector<int> x = position(grid, n) ;
        set(nuc_grid(n), 1) = 1.0 ;
        set(nuc_grid(n), 3) = 1.0 ;
        double phisum = 0.0 ;
        for(int h = 0 ; h < length(grid(x)) ; h++)
        {
            int hindex = MMSP::index(grid(x),h) ;
            if(hindex <= 12)
            {
                phisum += grid(n)[hindex] ;
            }
        }
        if(phisum > 0.1)
        {
            set(intenergies(n), 1) = 0 ;
            set(intenergies(n), 3) = 0 ;
            set(selfenergies(n), 1) = 0 ;
            set(selfenergies(n), 3) = 0 ;
            set(nuc_grid(n), 1) = 0.0 ;
            set(nuc_grid(n), 3) = 0.0 ;            
        }
        else
        {
            set(intenergies(n), 1) = 4*3.14*9*dx*dx*nuc_grid(n)[1]*(elint_real1[x[0]][x[1]]);
            set(intenergies(n), 3) = 4*3.14*9*dx*dx*nuc_grid(n)[3]*(elint_real3[x[0]][x[1]]);
            set(selfenergies(n), 1) = 4*3.14*9*dx*dx*nuc_grid(n)[1]*nuc_grid(n)[1]*(self_real1[x[0]][x[1]]);
            set(selfenergies(n), 3) = 4*3.14*9*dx*dx*nuc_grid(n)[3]*nuc_grid(n)[3]*(self_real3[x[0]][x[1]]);
        }
    }
    
    //cout<<"Finished nuc"<<endl;
      
    for(int i = 0 ; i < nx ; i++)
    {
      for(int j = 0 ; j < nyp ; j++)
      {
          double k1 = fk[i] ;
          double k2 = fk[j] ;
		
          double modk = (sqrt(k1*k1 + k2*k2)) ;
          if(modk==0.0) 
          {
              dfdstr1(i,j) = 0.0   ;
              dfdstr2(i,j) = 0.0   ;
              dfdstr3(i,j) = 0.0   ;
              dfdstr4(i,j) = 0.0   ;
              dfdstr5(i,j) = 0.0   ;
              dfdstr6(i,j) = 0.0   ;
              dfdstr7(i,j) = 0.0   ;
              dfdstr8(i,j) = 0.0   ;
              dfdstr9(i,j) = 0.0   ;
              dfdstr10(i,j) = 0.0   ;
              dfdstr11(i,j) = 0.0   ;
              dfdstr12(i,j) = 0.0   ;
              fstr(i,j) = 0.0 ;
          }
          else 
          {
              dfdstr1(i,j) = Bpq(k1, k2, modk, 1, 1)*F1(i,j) + Bpq(k1, k2, modk, 1, 2)*F2(i,j) + Bpq(k1, k2, modk, 1, 3)*F3(i,j) + Bpq(k1, k2, modk, 1, 4)*F4(i,j) 
                        + Bpq(k1, k2, modk, 1, 5)*F5(i,j) + Bpq(k1, k2, modk, 1, 6)*F6(i,j) + Bpq(k1, k2, modk, 1, 7)*F7(i,j) + Bpq(k1, k2, modk, 1, 8)*F8(i,j) + Bpq(k1, k2, modk, 1, 9)*F9(i,j)
                        + Bpq(k1, k2, modk, 1, 10)*F10(i,j) + Bpq(k1, k2, modk, 1, 11)*F11(i,j) + Bpq(k1, k2, modk, 1, 12)*F12(i,j) ;
              dfdstr2(i,j) = Bpq(k1, k2, modk, 2, 1)*F1(i,j) + Bpq(k1, k2, modk, 2, 2)*F2(i,j) + Bpq(k1, k2, modk, 2, 3)*F3(i,j) + Bpq(k1, k2, modk, 2, 4)*F4(i,j) 
                        + Bpq(k1, k2, modk, 2, 5)*F5(i,j) + Bpq(k1, k2, modk, 2, 6)*F6(i,j) + Bpq(k1, k2, modk, 2, 7)*F7(i,j) + Bpq(k1, k2, modk, 2, 8)*F8(i,j) + Bpq(k1, k2, modk, 2, 9)*F9(i,j)
                        + Bpq(k1, k2, modk, 2, 10)*F10(i,j) + Bpq(k1, k2, modk, 2, 11)*F11(i,j) + Bpq(k1, k2, modk, 2, 12)*F12(i,j) ;
              dfdstr3(i,j) = Bpq(k1, k2, modk, 3, 1)*F1(i,j) + Bpq(k1, k2, modk, 3, 2)*F2(i,j) + Bpq(k1, k2, modk, 3, 3)*F3(i,j) + Bpq(k1, k2, modk, 3, 4)*F4(i,j) 
                        + Bpq(k1, k2, modk, 3, 5)*F5(i,j) + Bpq(k1, k2, modk, 3, 6)*F6(i,j) + Bpq(k1, k2, modk, 3, 7)*F7(i,j) + Bpq(k1, k2, modk, 3, 8)*F8(i,j) + Bpq(k1, k2, modk, 3, 9)*F9(i,j)
                        + Bpq(k1, k2, modk, 3, 10)*F10(i,j) + Bpq(k1, k2, modk, 3, 11)*F11(i,j) + Bpq(k1, k2, modk, 3, 12)*F12(i,j) ;
              dfdstr4(i,j) = Bpq(k1, k2, modk, 4, 1)*F1(i,j) + Bpq(k1, k2, modk, 4, 2)*F2(i,j) + Bpq(k1, k2, modk, 4, 3)*F3(i,j) + Bpq(k1, k2, modk, 4, 4)*F4(i,j) 
                        + Bpq(k1, k2, modk, 4, 5)*F5(i,j) + Bpq(k1, k2, modk, 4, 6)*F6(i,j) + Bpq(k1, k2, modk, 4, 7)*F7(i,j) + Bpq(k1, k2, modk, 4, 8)*F8(i,j) + Bpq(k1, k2, modk, 4, 9)*F9(i,j)
                        + Bpq(k1, k2, modk, 4, 10)*F10(i,j) + Bpq(k1, k2, modk, 4, 11)*F11(i,j) + Bpq(k1, k2, modk, 4, 12)*F12(i,j) ;
              dfdstr5(i,j) = Bpq(k1, k2, modk, 5, 1)*F1(i,j) + Bpq(k1, k2, modk, 5, 2)*F2(i,j) + Bpq(k1, k2, modk, 5, 3)*F3(i,j) + Bpq(k1, k2, modk, 5, 4)*F4(i,j) 
                        + Bpq(k1, k2, modk, 5, 5)*F5(i,j) + Bpq(k1, k2, modk, 5, 6)*F6(i,j) + Bpq(k1, k2, modk, 5, 7)*F7(i,j) + Bpq(k1, k2, modk, 5, 8)*F8(i,j) + Bpq(k1, k2, modk, 5, 9)*F9(i,j)
                        + Bpq(k1, k2, modk, 5, 10)*F10(i,j) + Bpq(k1, k2, modk, 5, 11)*F11(i,j) + Bpq(k1, k2, modk, 5, 12)*F12(i,j) ;
              dfdstr6(i,j) = Bpq(k1, k2, modk, 6, 1)*F1(i,j) + Bpq(k1, k2, modk, 6, 2)*F2(i,j) + Bpq(k1, k2, modk, 6, 3)*F3(i,j) + Bpq(k1, k2, modk, 6, 4)*F4(i,j) 
                        + Bpq(k1, k2, modk, 6, 5)*F5(i,j) + Bpq(k1, k2, modk, 6, 6)*F6(i,j) + Bpq(k1, k2, modk, 6, 7)*F7(i,j) + Bpq(k1, k2, modk, 6, 8)*F8(i,j) + Bpq(k1, k2, modk, 6, 9)*F9(i,j)
                        + Bpq(k1, k2, modk, 6, 10)*F10(i,j) + Bpq(k1, k2, modk, 6, 11)*F11(i,j) + Bpq(k1, k2, modk, 6, 12)*F12(i,j) ;
              dfdstr7(i,j) = Bpq(k1, k2, modk, 7, 1)*F1(i,j) + Bpq(k1, k2, modk, 7, 2)*F2(i,j) + Bpq(k1, k2, modk, 7, 3)*F3(i,j) + Bpq(k1, k2, modk, 7, 4)*F4(i,j) 
                        + Bpq(k1, k2, modk, 7, 5)*F5(i,j) + Bpq(k1, k2, modk, 7, 6)*F6(i,j) + Bpq(k1, k2, modk, 7, 7)*F7(i,j) + Bpq(k1, k2, modk, 7, 8)*F8(i,j) + Bpq(k1, k2, modk, 7, 9)*F9(i,j)
                        + Bpq(k1, k2, modk, 7, 10)*F10(i,j) + Bpq(k1, k2, modk, 7, 11)*F11(i,j) + Bpq(k1, k2, modk, 7, 12)*F12(i,j) ;
              dfdstr8(i,j) = Bpq(k1, k2, modk, 8, 1)*F1(i,j) + Bpq(k1, k2, modk, 8, 2)*F2(i,j) + Bpq(k1, k2, modk, 8, 3)*F3(i,j) + Bpq(k1, k2, modk, 8, 4)*F4(i,j) 
                        + Bpq(k1, k2, modk, 8, 5)*F5(i,j) + Bpq(k1, k2, modk, 8, 6)*F6(i,j) + Bpq(k1, k2, modk, 8, 7)*F7(i,j) + Bpq(k1, k2, modk, 8, 8)*F8(i,j) + Bpq(k1, k2, modk, 8, 9)*F9(i,j)
                        + Bpq(k1, k2, modk, 8, 10)*F10(i,j) + Bpq(k1, k2, modk, 8, 11)*F11(i,j) + Bpq(k1, k2, modk, 8, 12)*F12(i,j) ;
              dfdstr9(i,j) = Bpq(k1, k2, modk, 9, 1)*F1(i,j) + Bpq(k1, k2, modk, 9, 2)*F2(i,j) + Bpq(k1, k2, modk, 9, 3)*F3(i,j) + Bpq(k1, k2, modk, 9, 4)*F4(i,j) 
                        + Bpq(k1, k2, modk, 9, 5)*F5(i,j) + Bpq(k1, k2, modk, 9, 6)*F6(i,j) + Bpq(k1, k2, modk, 9, 7)*F7(i,j) + Bpq(k1, k2, modk, 9, 8)*F8(i,j) + Bpq(k1, k2, modk, 9, 9)*F9(i,j)
                        + Bpq(k1, k2, modk, 9, 10)*F10(i,j) + Bpq(k1, k2, modk, 9, 11)*F11(i,j) + Bpq(k1, k2, modk, 9, 12)*F12(i,j) ;
              dfdstr10(i,j) = Bpq(k1, k2, modk, 10, 1)*F1(i,j) + Bpq(k1, k2, modk, 10, 2)*F2(i,j) + Bpq(k1, k2, modk, 10, 3)*F3(i,j) + Bpq(k1, k2, modk, 10, 4)*F4(i,j) 
                        + Bpq(k1, k2, modk, 10, 5)*F5(i,j) + Bpq(k1, k2, modk, 10, 6)*F6(i,j) + Bpq(k1, k2, modk, 10, 7)*F7(i,j) + Bpq(k1, k2, modk, 10, 8)*F8(i,j) + Bpq(k1, k2, modk, 10, 9)*F9(i,j)
                        + Bpq(k1, k2, modk, 10, 10)*F10(i,j) + Bpq(k1, k2, modk, 10, 11)*F11(i,j) + Bpq(k1, k2, modk, 10, 12)*F12(i,j) ;
              dfdstr11(i,j) = Bpq(k1, k2, modk, 11, 1)*F1(i,j) + Bpq(k1, k2, modk, 11, 2)*F2(i,j) + Bpq(k1, k2, modk, 11, 3)*F3(i,j) + Bpq(k1, k2, modk, 11, 4)*F4(i,j) 
                        + Bpq(k1, k2, modk, 11, 5)*F5(i,j) + Bpq(k1, k2, modk, 11, 6)*F6(i,j) + Bpq(k1, k2, modk, 11, 7)*F7(i,j) + Bpq(k1, k2, modk, 11, 8)*F8(i,j) + Bpq(k1, k2, modk, 11, 9)*F9(i,j)
                        + Bpq(k1, k2, modk, 11, 10)*F10(i,j) + Bpq(k1, k2, modk, 11, 11)*F11(i,j) + Bpq(k1, k2, modk, 11, 12)*F12(i,j) ;
              dfdstr12(i,j) = Bpq(k1, k2, modk, 12, 1)*F1(i,j) + Bpq(k1, k2, modk, 12, 2)*F2(i,j) + Bpq(k1, k2, modk, 12, 3)*F3(i,j) + Bpq(k1, k2, modk, 12, 4)*F4(i,j) 
                        + Bpq(k1, k2, modk, 12, 5)*F5(i,j) + Bpq(k1, k2, modk, 12, 6)*F6(i,j) + Bpq(k1, k2, modk, 12, 7)*F7(i,j) + Bpq(k1, k2, modk, 12, 8)*F8(i,j) + Bpq(k1, k2, modk, 12, 9)*F9(i,j)
                        + Bpq(k1, k2, modk, 12, 10)*F10(i,j) + Bpq(k1, k2, modk, 12, 11)*F11(i,j) + Bpq(k1, k2, modk, 12, 12)*F12(i,j) ;
                
              fstr(i,j) = Bpq(k1, k2, modk, 1, 1)*F1(i,j)*conj(F1(i,j)) ;
          }
      }
	}
    

  
  
  Backward1.fft0Normalized(dfdstr1, dfdstr_real1); 
  Backward2.fft0Normalized(dfdstr2, dfdstr_real2); 
  Backward3.fft0Normalized(dfdstr3, dfdstr_real3); 
  Backward4.fft0Normalized(dfdstr4, dfdstr_real4); 
  Backward5.fft0Normalized(dfdstr5, dfdstr_real5); 
  Backward6.fft0Normalized(dfdstr6, dfdstr_real6); 
  Backward7.fft0Normalized(dfdstr7, dfdstr_real7); 
  Backward8.fft0Normalized(dfdstr8, dfdstr_real8); 
  Backward9.fft0Normalized(dfdstr9, dfdstr_real9); 
  Backward10.fft0Normalized(dfdstr10, dfdstr_real10); 
  Backward11.fft0Normalized(dfdstr11, dfdstr_real11); 
  Backward12.fft0Normalized(dfdstr12, dfdstr_real12); 
  
  Backward17.fft0Normalized(fstr, fstr_real); 
  
    
        if(t%100==0)   //Generating output files. Refer to outputgen.hpp
		{
  
            
            std::string file_name = "output_" + to_string(t) + ".txt" ;
            ofstream myfile;
            myfile.open(file_name.c_str());
            //intstrain();
            for(int n = 0 ; n < nodes(grid) ; n++)    //for n < nodes(grid) ; 
            {
                MMSP::vector<int> s = position(grid,n) ;
                double sum = 0 ;
                for(int h = 0 ; h < length(grid(s)) ; h++)
                {
                    int hindex = MMSP::index(grid(s),h) ;
                    if(hindex <= 12)
                    {
                        sum+=grid(s)[hindex]*grid(s)[hindex] ;
                    }
                }
                
                double sfts_sum = 0 ;
                for(int i = 0 ; i < 6 ; i++)
                {
                    for(int j = 0 ; j < 6 ; j++)
                    {
                        sfts_sum+= eigen_alpha[variant-1].e[i]*c[i][j]*eigen_alpha[variant-1].e[j] ;
                    }
                }
			
                double felastic = 0.5*grid(n)[1]*sfts_sum - 0.5*fstr_real[s[0]][s[1]]/(nx*ny) ;
			
                myfile<<s[0]<<","<<s[1]<<","<<sum<<","<<(dfdstr_real1[s[0]][s[1]])/scaling<<","<<intenergies(s)[1]/(nx*ny)<<","<<intenergies(s)[3]/(nx*ny)<<","<<fstr_real[s[0]][s[1]]/(nx*ny)<<"\n" ; 
			
			
            }
	
            myfile.close();
            
        }
            
            
        for(int n = 0 ; n<nodes(grid) ; n++)
		{
			MMSP::vector<int> x = position(grid, n); 
			MMSP::vector<store_type> gradientsqtemp = gradsq(grid, x) ;
			double temp1 = 0.0 ;
			double temp2 = 0.0 ;
			for(int i = 0 ; i < dim ; i++)
			{
					temp1 += gradientsqtemp[i][20] ;
					temp1 += gradientsqtemp[i][21] ;
			}
			set(gradsqcal_grid(n),1) = temp1 ;
			set(gradsqcv_grid(n), 1) = temp2 ;
		
		}
	
	
	
		for(int n = 0 ; n < nodes(grid) ; n++)
		{
			MMSP::vector<int> x = position(grid, n); 
			
			double strain_energy[13];
			strain_energy[0] = 0.0 ;
			strain_energy[1] = (double)(1.0/(nx*ny))*(dfdstr_real1[x[0]][x[1]]) ; 
            strain_energy[2] = (double)(1.0/(nx*ny))*(dfdstr_real2[x[0]][x[1]]) ; 
            strain_energy[3] = (double)(1.0/(nx*ny))*(dfdstr_real3[x[0]][x[1]]) ; 
            strain_energy[4] = (double)(1.0/(nx*ny))*(dfdstr_real4[x[0]][x[1]]) ; 
            strain_energy[5] = (double)(1.0/(nx*ny))*(dfdstr_real5[x[0]][x[1]]) ; 
            strain_energy[6] = (double)(1.0/(nx*ny))*(dfdstr_real6[x[0]][x[1]]) ; 
            strain_energy[7] = (double)(1.0/(nx*ny))*(dfdstr_real7[x[0]][x[1]]) ; 
            strain_energy[8] = (double)(1.0/(nx*ny))*(dfdstr_real8[x[0]][x[1]]) ;
            strain_energy[9] = (double)(1.0/(nx*ny))*(dfdstr_real9[x[0]][x[1]]) ; 
            strain_energy[10] = (double)(1.0/(nx*ny))*(dfdstr_real10[x[0]][x[1]]) ; 
            strain_energy[11] = (double)(1.0/(nx*ny))*(dfdstr_real11[x[0]][x[1]]) ; 
            strain_energy[12] = (double)(1.0/(nx*ny))*(dfdstr_real12[x[0]][x[1]]) ;
            
			 
			
			
			G_Alpha[n] = ((grid(n)[20]*G_Al_alpha + grid(n)[21]*G_V_alpha + (1-grid(n)[20]-grid(n)[21])*G_Ti_alpha +
					 R*T*(grid(n)[20]*log(grid(n)[20]) + grid(n)[21]*log(grid(n)[21]) + (1-grid(n)[20]-grid(n)[21])*log((1-grid(n)[20]-grid(n)[21]))) +
					 grid(n)[20]*grid(n)[21]*L0_HCP_Al_V + grid(n)[20]*(1-grid(n)[20]-grid(n)[21])*L0_HCP_Al_Ti + grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_HCP_Ti_V))/G_normalize ;
					 
			G_Beta[n] = (grid(n)[20]*G_Al_beta + grid(n)[21]*G_V_beta + (1-grid(n)[20]-grid(n)[21])*G_Ti_beta +
					 R*T*(grid(n)[20]*log(grid(n)[20]) + grid(n)[21]*log(grid(n)[21]) + (1-grid(n)[20]-grid(n)[21])*log((1-grid(n)[20]-grid(n)[21]))) +
					 grid(n)[20]*grid(n)[21]*L0_BCC_Al_V + grid(n)[20]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Al_Ti + grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Ti_V +
					 grid(n)[20]*grid(n)[21]*(1-grid(n)[20]-grid(n)[21])*L0_BCC_Al_Ti_V)/G_normalize ;
					 
			W[n] = W_prefac*(grid(n)[20]*W_Al + grid(n)[21]*W_V + (1-grid(n)[20]-grid(n)[21])*W_Ti) ;
		
		
			MMSP::vector<store_type> gradient = grad(grid, x) ;
			MMSP::vector<store_type> gradientsq = gradsq(grid, x) ;

		
			thermo_auxillary_terms(gradient, gradientsq, grid(n)[20], grid(n)[21]) ;

			MMSP::vector<store_type> gradcalp4 = gradsq(gradsqcal_grid, x) ;	
			MMSP::vector<store_type> gradcvp4 = gradsq(gradsqcv_grid, x) ;
		
			double gradp4_cal = 0 ;
			double gradp4_cv  = 0 ;
		
			for(int i = 0 ; i < dim ; i++)
			{
				gradp4_cal += gradcalp4[i][1]  ;
				gradp4_cv += gradcvp4[i][1] ;
			}
		
			double phi_grad[12] ;
			double phi_gradientsq[12] ;
			for(int h = 0 ; h < length(grid(x)) ; h++)
			{
				int hindex = MMSP::index(grid(x),h) ;
				if(hindex <= 12)
				{
					phi_grad[hindex-1] = 0 ;
					phi_gradientsq[hindex-1] = 0 ;
					for(int i = 0 ; i < dim ; i++)
					{
						phi_grad[hindex-1]+= gradient[i][hindex] ; 
						phi_gradientsq[hindex-1]+= gradientsq[i][hindex] ; 
					}
				}
			}
		
		
		
			for(int h = 0 ; h < length(grid(x)) ; h++)
			{
				int hindex = MMSP::index(grid(x),h) ;
				if(hindex <= 12)
				{
				
					/*hphidoubleprime[n][hindex-1] = 60*grid(n)[hindex]*(2*grid(n)[hindex]-1)*(grid(n)[hindex]-1) ; 
					hphi[n][hindex-1] = pow(grid(n)[hindex],3)*(6*pow(grid(n)[hindex],2) - 15*grid(n)[hindex] + 10) ;
					hphiprime[n][hindex-1] = 30*pow(grid(n)[hindex],2)*pow((grid(n)[hindex]-1),2) ;*/
					
					hphidoubleprime[n][hindex-1] = 2*pow(grid(n)[hindex],1) - 3*pow(grid(n)[hindex],2) ; 
					hphi[n][hindex-1] = pow(grid(n)[hindex],3)/3 - pow(grid(n)[hindex],4)/4 ;
					hphiprime[n][hindex-1] = pow(grid(n)[hindex],2) - pow(grid(n)[hindex],3) ;
		
					gphidoubleprime[n][hindex-1] = 2*(6*pow(grid(n)[hindex],2) - 6*grid(n)[hindex] + 1) ; 
					gphiprime[n][hindex-1] = 2*grid(n)[hindex]*(1-grid(n)[hindex])*(1-2*grid(n)[hindex]) ;
					gphi[n][hindex-1] = pow(grid(n)[hindex],2)*pow((1-grid(n)[hindex]),2) ;
				
				
				}
			}
		
			double hphidoubleprimesum = 0;
			double hphiprimesum1 = 0;
			double hphiprimesum2 = 0;
			double hphisum = 0;
			double gphidoubleprimesum = 0;
			double gphiprimesum = 0;
			double gphisum = 0;
		
			for(int h = 0 ; h < length(grid(x)) ; h++)
			{
				int hindex = MMSP::index(grid(x),h) ;
				if(hindex <= 12)
				{
					hphidoubleprimesum += hphidoubleprime[n][hindex-1]*pow(phi_grad[hindex-1],2) ;
					hphiprimesum1 += hphiprime[n][hindex-1]*phi_gradientsq[hindex-1];
					hphiprimesum2 += hphiprime[n][hindex-1]*phi_grad[hindex-1] ;
					hphisum += hphi[n][hindex-1] ;
				
					gphiprimesum += gphiprime[n][hindex-1]*phi_gradientsq[hindex-1] ;
					gphidoubleprimesum += gphidoubleprime[n][hindex-1]*phi_grad[hindex-1] ;
				}
			}
		
												 
			double c_al_rhs = 2*(del_dGAlpha_dAl - del_dGBeta_dAl)*hphiprimesum2 + delsq_dGAlpha_dAl*hphisum + delsq_dGBeta_dAl*(1-hphisum) + 
						  (dGAlpha_dAl - dGBeta_dAl)*(hphidoubleprimesum + hphiprimesum1) + (W_Al - W_Ti)*(gphiprimesum + gphidoubleprimesum) ;
			double c_v_rhs = 2*(del_dGAlpha_dV - del_dGBeta_dV)*hphiprimesum2 + delsq_dGAlpha_dV*hphisum + delsq_dGBeta_dV*(1-hphisum) + 
						  (dGAlpha_dV - dGBeta_dV)*(hphidoubleprimesum + hphiprimesum1) + (W_V - W_Ti)*(gphiprimesum + gphidoubleprimesum)  ;

			set(update(n), 20) = grid(n)[20] + dt*(Dalal*(c_al_rhs)- kappa_c*gradp4_cal + Dalv*(c_v_rhs)- 0.5*kappa_c*gradp4_cv) ; 
			set(update(n), 21) = grid(n)[21] + dt*(Dvv*(c_v_rhs) - kappa_c*gradp4_cv + Dval*(c_al_rhs)- 0.5*kappa_c*gradp4_cal) ; 
		
			double lap_aniso[13];	
			
			
			
			MMSP::sparse<phi_type> dFdp;
			phi_type dFall = 0.0;		
			phi_type dFall_no = 0 ;	
			phi_type phi_sum = 0.0 ;
			
			for (int j = 0; j < length(grid(n)); j++) 
			{
				int jindex = MMSP::index(grid(n), j);
				if(jindex<=12)
				{		
					phi_sum += grid(n)[jindex] ;
				}
			}
				
			phi_type phi_beta = 1 - phi_sum ; 
			
			for (int j = 0; j < length(grid(n)); j++)   
			{
				int jindex = MMSP::index(grid(n), j);  		
					
				MMSP::vector<int> x = position(grid, n);
				
				W[n] = W_prefac*(grid(n)[20]*W_Al + grid(n)[21]*W_V + (1-grid(n)[20]-grid(n)[21])*W_Ti) ;

				if(jindex <= 12)
				{
					phi_type lap_aniso = 0.0 ;
					lap_aniso+= (gradientsq[0][jindex]*epsi[jindex-1][0][0] + gradientsq[1][jindex]*epsi[jindex-1][1][1])  ;
					
					
				
					
					double interaction_energy = 0 ;
					interaction_energy+= W[n]*grid(n)[jindex]*pow(phi_beta, 2);
					interaction_energy-= W[n]*pow(grid(n)[jindex],2)*(phi_beta);
					
					int check = 0;
					for(int k = 0 ; k < length(grid(n)); k++)
					{
						int tempindex = MMSP::index(grid(n), k); 
						if((tempindex!=jindex)&&(tempindex <=12))
						{
							interaction_energy+= W[n]*pow(grid(n)[tempindex],2)*grid(n)[jindex];
							interaction_energy-= W[n]*pow(grid(n)[tempindex],2)*(phi_beta); 
						}
					}
					
					
					double Gdiff = G_Alpha[n] - G_Beta[n];
					
					set(dFdp, jindex) = -1*hphiprime[n][jindex-1] + strain_energy[jindex] - lap_aniso ;  
					
					
					
					dFall += dFdp[jindex];	
					dFall_no+=1;
				}
										
			}
			L = 0.1;

			for (int h = 0; h < length(grid(n)); h++) 
			{	
					
				int hindex = MMSP::index(grid(n), h);   
				if(hindex<=12)
				{
					store_type dpdt;
					set(dpdt, hindex) = -L * (dFdp[hindex]);
					//set(dpdt, hindex) = -L * ((dFall_no+1)*dFdp[hindex] - dFall);
					set(update(n), hindex) = grid(n)[hindex] + dt * dpdt[hindex];	
				}
			}
			
		}

		swap(grid, update);	

}
  

  
  

    



return 0 ;

}







