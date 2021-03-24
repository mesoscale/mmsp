using namespace std ;

double c[6][6], sigma00[12][6] ;
double c_norm = 1;

double omega_inv_ij(double* cnew, double* k ,  double modk)
{
	double sum = 0.0 ;
	for(int i = 0 ; i < 4 ; i++)
	{
		sum += cnew[i]*k[i] ;
	}
	return (sum)/(modk*modk) ; 
}

void define_c_sigma()
{

	for(int i = 0 ; i < 6 ; i++)
	{
		for(int j = 0 ; j < 6 ; j++)
		{
			c[i][j] = 0.0 ;   
		}
	}
    double scale = 1.0 ;
    c[0][0] = 175.0/scale;
    c[1][1] = 175.0/scale;
    c[2][2] = 220.0/scale;
    c[0][1] = 88.7/scale ;
    c[0][2] = 62.3/scale ;
    c[1][2] = 62.3/scale ;
    c[3][3] = 62.2/scale ;
    c[4][4] = 62.2/scale ;
    c[5][5] = (c[0][0]-c[0][1])/2.0 ;
    c[5][5] /= scale ;



	for(int variant = 0 ; variant < 12 ; variant++)
	{
		for(int i = 0 ; i < 6 ; i++)
		{
			sigma00[variant][i] = 0.0 ;
			for(int j = 0 ; j < 6 ; j++)
			{
				sigma00[variant][i] += c[i][j]*eigen_alpha[variant].e[j] ;
			}
		}			
	} 



}


double B(double k1, double k2, double modk, int variant) {
	double c11[4] = {c[0][0], c[0][5], c[5][0], c[5][5]} ; 
	double c12[4] = {c[0][5], c[0][1], c[5][5], c[5][1]} ; 
	double c21[4] = {c[5][0], c[5][5], c[1][0], c[1][5]} ; 
	double c22[4] = {c[5][5], c[5][1], c[1][5], c[1][1]} ; 
	double k[4] = {k1*k1, k1*k2, k2*k1, k2*k2} ;


	double omega_inv[2][2] = {{omega_inv_ij(c11, k, modk), omega_inv_ij(c12, k, modk)}, {omega_inv_ij(c21, k, modk), omega_inv_ij(c22, k, modk)}} ;


	double det_omega_inv = omega_inv[0][0]*omega_inv[1][1] - omega_inv[1][0]*omega_inv[0][1] + 1e-05  ;
	if(det_omega_inv==0) cout<<"Found zero det"<<endl ;
	
					   
	double omega[3][3] ; 

	for(int i = 0 ; i < dim ; i++)
	{
		for(int j = 0 ; j < dim ; j++)
		{
            if(i==j)
            {
                omega[i][j] =  (1.0/det_omega_inv)*omega[j][i] ;  
            }
            else
            {
                omega[i][j] =  -(1.0/det_omega_inv)*omega[i][j] ;
            }
		}
	}
    
	double k1l = sigma00[variant][0]*k1/modk + sigma00[variant][5]*k2/modk ; 
	double k2l = sigma00[variant][5]*k1/modk + sigma00[variant][1]*k2/modk  ; 
	
	double j1kl = omega[0][0]*k1l + omega[0][1]*k2l  ;
	double j2kl = omega[1][0]*k1l + omega[1][1]*k2l  ;
	
	double i1jkl = (k1/modk)*(sigma00[variant][0]*j1kl + sigma00[variant][5]*j2kl) ;
	double i2jkl = (k2/modk)*(sigma00[variant][5]*j1kl + sigma00[variant][1]*j2kl)  ;
	
	double sfts_sum = 0 ;
	for(int i = 0 ; i < 6 ; i++)
	{
        if(i!=0 or i!=1 or i!=5) continue;
		for(int j = 0 ; j < 6 ; j++)
		{
            if(j!=0 or j!=1 or j!=5) continue;
			sfts_sum+= eigen_alpha[variant-1].e[i]*c[i][j]*eigen_alpha[variant-1].e[j] ;
		}
	}
    
    return (-(i1jkl + i2jkl)) ; }
    
double Bpq(double k1, double k2, double modk, int p, int q) {
	
	double c11[4] = {c[0][0], c[0][5], c[5][0], c[5][5]} ; 
	double c12[4] = {c[0][5], c[0][1], c[5][5], c[5][1]} ; 
	double c21[4] = {c[5][0], c[5][5], c[1][0], c[1][5]} ; 
	double c22[4] = {c[5][5], c[5][1], c[1][5], c[1][1]} ; 
	double k[4] = {k1*k1, k1*k2, k2*k1, k2*k2} ;


	double omega_inv[2][2] = {{omega_inv_ij(c11, k, modk), omega_inv_ij(c12, k, modk)}, {omega_inv_ij(c21, k, modk), omega_inv_ij(c22, k, modk)}} ;
    /*for(int i = 0 ; i < dim ; i++)
	{
		for(int j = 0 ; j < dim ; j++)
		{
            cout<<omega_inv[i][j]<<endl;
        }
    }*/

    
	double det_omega_inv = omega_inv[0][0]*omega_inv[1][1] - omega_inv[1][0]*omega_inv[0][1] + 1e-05  ;
	if(det_omega_inv==0) cout<<"Found zero det"<<endl ;
	
    //cout<<det_omega_inv<<endl;
	double omega[3][3] ; 

	for(int i = 0 ; i < dim ; i++)
	{
		for(int j = 0 ; j < dim ; j++)
		{
            if(i==j)
            {
                omega[i][j] =  (1.0/det_omega_inv)*omega_inv[j][i] ;  
            }
            else
            {
                omega[i][j] =  -(1.0/det_omega_inv)*omega_inv[i][j] ;
            }
		}
	}
    
    /*for(int i = 0 ; i < dim ; i++)
	{
		for(int j = 0 ; j < dim ; j++)
		{
            cout<<omega[i][j]<<endl;
        }
    }*/
    
    double k1l = sigma00[q-1][0]*k1/modk + sigma00[q-1][5]*k2/modk ; 
	double k2l = sigma00[q-1][5]*k1/modk + sigma00[q-1][1]*k2/modk  ; 
	
	double j1kl = omega[0][0]*k1l + omega[0][1]*k2l  ;
	double j2kl = omega[1][0]*k1l + omega[1][1]*k2l  ;
	
	double i1jkl = (k1/modk)*(sigma00[p-1][0]*j1kl + sigma00[p-1][5]*j2kl) ;
	double i2jkl = (k2/modk)*(sigma00[p-1][5]*j1kl + sigma00[p-1][1]*j2kl)  ;
	
	double sfts_sum = 0 ;
	for(int i = 0 ; i < 6 ; i++)
	{
		for(int j = 0 ; j < 6 ; j++)
		{
			sfts_sum+= eigen_alpha[p-1].e[i]*c[i][j]*eigen_alpha[q-1].e[j] ;
		}
	}
			
    
    return (sfts_sum - (i1jkl + i2jkl)) ; } 
    

