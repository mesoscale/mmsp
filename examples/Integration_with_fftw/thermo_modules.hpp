


#define L0_BCC_Al_Ti_V (32045.963)
#define L1_BCC_Al_Ti_V (-113926.0 + 40*T) ;
#define L2_BCC_Al_Ti_V (75972.5 - 150*T) ;
#define L0_HCP_Al_Ti_V 0 
#define L1_HCP_Al_Ti_V (-206074 - 40*T) ;
#define L2_HCP_Al_Ti_V 0 ;
#define L0_BCC_Ti_V (10500.0-1.5*T) 
#define L1_BCC_Ti_V 2025.39 ;
#define L0_HCP_Ti_V 20000.0/(G_normalize)
#define L0_BCC_Al_Ti (-128500.0+39.0*T+2520.0)
#define L1_BCC_Al_Ti 4890.0 ;
#define L2_BCC_Al_Ti 400.0 ;
#define L0_HCP_Al_Ti (-133500.0+39.0*T+2950.0)
#define L1_HCP_Al_Ti -3475.0 + 0.825*T ;
#define L2_HCP_Al_Ti -7756.0 ;
#define L0_HCP_Al_V (-95000.0+20.0*T)
#define L0_BCC_Al_V (-95000.0+20.0*T)
#define L1_BCC_Al_V 6645.0 ;
#define L2_BCC_Al_V -68596.0 ;

#define Al_alpha_a_1 -2495.15 
#define Al_alpha_b_1 135.29 
#define Al_alpha_c_1 -24.37 
#define Al_alpha_d_1 -0.00188 

#define Al_alpha_a_2 -5795.24 
#define Al_alpha_b_2 221.25 
#define Al_alpha_c_2 -38.58 
#define Al_alpha_d_2 -0.01853 

#define Al_alpha_a_3 -5797.36 
#define Al_alpha_b_3 186.88 
#define Al_alpha_c_3 -31.75 
#define Al_alpha_d_3 0.0 

#define Al_beta_a_1 2106.85 ;
#define Al_beta_b_1 132.28 ;
#define Al_beta_c_1 -24.37 ;
#define Al_beta_d_1 -0.00188 ;

#define Al_beta_a_2 -1193.24 ;
#define Al_beta_b_2 218.24 ;
#define Al_beta_c_2 -38.58 ;
#define Al_beta_d_2 -0.01853 ;

#define Al_beta_a_3 -1195.36 
#define Al_beta_b_3 183.87 
#define Al_beta_c_3 -31.75 
#define Al_beta_d_3 0.0 

#define Ti_alpha_a_1 -8059.92 ;
#define Ti_alpha_b_1 133.62 ;
#define Ti_alpha_c_1 -23.99 ;
#define Ti_alpha_d_1 -0.00477 ;

#define Ti_alpha_a_2 -7811.82 
#define Ti_alpha_b_2 132.98 
#define Ti_alpha_c_2 -23.98 
#define Ti_alpha_d_2 -0.0042 

#define Ti_alpha_a_3 908.84 ;  
#define Ti_alpha_b_3 66.98 ;
#define Ti_alpha_c_3 -14.95 ;
#define Ti_alpha_d_3 -0.00815 ;

#define Ti_alpha_a_4 -124526.79 ;  
#define Ti_alpha_b_4 638.81 ;
#define Ti_alpha_c_4 -87.22 ;
#define Ti_alpha_d_4 -0.00821 ;

#define Ti_beta_a_1 -1272.06 
#define Ti_beta_b_1 134.71 
#define Ti_beta_c_1 -25.58 
#define Ti_beta_d_1 -0.00066 

#define Ti_beta_a_2 6667.39 ;
#define Ti_beta_b_2 105.37 ;
#define Ti_beta_c_2 -22.37 ;
#define Ti_beta_d_2 0.00122 ;

#define Ti_beta_a_3 26483.26 ;
#define Ti_beta_b_3 -182.43 ;
#define Ti_beta_c_3 19.09 ;
#define Ti_beta_d_3 -22.01 ;


#define V_alpha_a_1 -3930.43 ;
#define V_alpha_b_1 135.74 ; 
#define V_alpha_c_1 -24.13 ;
#define V_alpha_d_1 -0.0031 ;

#define V_alpha_a_2 -3967.84 
#define V_alpha_b_2 145.69 
#define V_alpha_c_2 -25.90 
#define V_alpha_d_2 0.000063 

#define V_alpha_a_3 -37689.86 ;
#define V_alpha_b_3 323.54 ;
#define V_alpha_c_3 -47.43 ;
#define V_alpha_d_3 0.0 ;

#define V_beta_a_1 -7930.43 ;
#define V_beta_b_1 133.35 ;
#define V_beta_c_1 -24.13 ;
#define V_beta_d_1 -0.0031 ;

#define V_beta_a_2 -7967.84 
#define V_beta_b_2 143.29 
#define V_beta_c_2 -25.90 
#define V_beta_d_2 0.000063 

#define V_beta_a_3 -41689.87 ;
#define V_beta_b_3 321.14 ;
#define V_beta_c_3 -47.43 ;
#define V_beta_d_3 0.0 ;

void thermo_auxillary_terms(MMSP::vector<store_type> gradient, MMSP::vector<store_type> gradientsq, double c_Al, double c_V)
{
	double grad_cal = gradient[0][20] ;// + gradient[1][20] ;// + gradient[2][20] ;
	double grad_cv = gradient[0][21] ;// + gradient[1][21]  ;//+ gradient[2][21] ;
	double gradsq_cal = gradientsq[0][20] ;// + gradientsq[1][20] ;//+ gradientsq[2][20] ;
	double gradsq_cv = gradientsq[0][21] ;//+ gradientsq[1][21] ;//+ gradientsq[2][21] ;
	dGAlpha_dAl = (G_Al_alpha/G_normalize - G_Ti_alpha/G_normalize) + log(c_Al) - log((1-c_Al-c_V)) + c_V*L0_HCP_Al_V/G_normalize + 
							(1-2*c_Al-c_V)*L0_HCP_Al_Ti/G_normalize - c_V*L0_HCP_Ti_V/G_normalize ;
	dGBeta_dAl = (G_Al_beta/G_normalize - G_Ti_beta/G_normalize) + log(c_Al) - log((1-c_Al-c_V)) + c_V*L0_BCC_Al_V/G_normalize +
						    (1-2*c_Al-c_V)*L0_BCC_Al_Ti/G_normalize - c_V*L0_BCC_Ti_V/G_normalize + (c_V*(1-c_Al-c_V)
						     - c_Al*c_V)*L0_BCC_Al_Ti_V/G_normalize ;
	dGAlpha_dV = (G_V_alpha/G_normalize - G_Ti_alpha/G_normalize) + log(c_V) - log((1-c_Al-c_V)) + c_Al*L0_HCP_Al_V/G_normalize +
							(1-c_Al-2*c_V)*L0_HCP_Ti_V/G_normalize - c_Al*L0_HCP_Al_Ti/G_normalize ;	
	dGBeta_dV = (G_V_beta/G_normalize - G_Ti_beta/G_normalize) + log(c_V) - log((1-c_Al-c_V)) + c_Al*L0_BCC_Al_V/G_normalize +
							(1-c_Al-2*c_V)*L0_BCC_Ti_V/G_normalize - c_Al*L0_BCC_Al_Ti/G_normalize + (c_Al*(1-c_Al-c_V) -
							c_Al*c_V)*L0_BCC_Al_Ti_V/G_normalize ;
	del_dGAlpha_dAl = grad_cal*(1.0/c_Al + 1.0/(1.0-c_Al-c_V) - 2*L0_HCP_Al_Ti/G_normalize) + grad_cv*(1.0/(1.0-c_Al-c_V) +
							L0_HCP_Al_V/G_normalize - L0_HCP_Al_Ti/G_normalize -L0_HCP_Ti_V/G_normalize) ;
	del_dGBeta_dAl = grad_cal*(1.0/c_Al + 1/(1-c_Al-c_V) - 2*L0_BCC_Al_Ti/G_normalize - 2*c_V*L0_BCC_Al_Ti_V/G_normalize) + 
								grad_cv*(1/(1-c_Al-c_V) + L0_HCP_Al_V/G_normalize - L0_HCP_Al_Ti/G_normalize -L0_HCP_Ti_V/G_normalize +
								(1-2*c_Al-2*c_V)*L0_BCC_Al_Ti_V/G_normalize) ;	
	del_dGAlpha_dV = grad_cv*(1.0/c_V + 1/(1-c_Al-c_V) - 2*L0_HCP_Ti_V/G_normalize) + grad_cal*(1/(1-c_Al-c_V) + 
						    L0_HCP_Al_V/G_normalize - L0_HCP_Al_Ti/G_normalize - L0_HCP_Ti_V/G_normalize) ;	
	del_dGBeta_dV = 	grad_cv*(1.0/c_V + 1/(1-c_Al-c_V) -2*L0_BCC_Ti_V/G_normalize - 2*c_Al*L0_BCC_Al_Ti_V/G_normalize) + 
							    grad_cal*(1/(1-c_Al-c_V) + L0_BCC_Al_V/G_normalize - L0_BCC_Al_Ti/G_normalize - L0_BCC_Ti_V/G_normalize + 
							    (1-2*c_Al-2*c_V)*L0_BCC_Al_Ti_V/G_normalize) ;										    
	delsq_dGAlpha_dAl = gradsq_cal*(1.0/c_Al + 1/(1-c_Al-c_V) - 2*L0_HCP_Al_Ti/G_normalize) + gradsq_cv*(1/(1-c_Al-c_V) +
								L0_HCP_Al_V/G_normalize - L0_HCP_Al_Ti/G_normalize -L0_HCP_Ti_V/G_normalize) + grad_cal*(-grad_cal/pow(c_Al,2) +
								(grad_cal + grad_cv)/pow((1-c_Al-c_V),2)) + grad_cv*((grad_cal + grad_cv)/pow((1-c_Al-c_V),2)) ;
					
	delsq_dGAlpha_dV = gradsq_cv*(1.0/c_V + 1/(1-c_Al-c_V) - 2*L0_HCP_Ti_V/G_normalize) + gradsq_cal*(1/(1-c_Al-c_V) +
									L0_HCP_Al_V/G_normalize - L0_HCP_Al_Ti/G_normalize - L0_HCP_Ti_V/G_normalize) + grad_cv*(-grad_cv/pow(c_V,2) +
									(grad_cal + grad_cv)/pow((1-c_Al-c_V),2)) + grad_cal*((grad_cal + grad_cv)/pow((1-c_Al-c_V),2)) ;					
							
	delsq_dGBeta_dAl = gradsq_cal*(1.0/c_Al + 1/(1-c_Al-c_V) - 2*L0_BCC_Al_Ti/G_normalize - 2*c_V*L0_BCC_Al_Ti_V/G_normalize) +
								  gradsq_cv*(1/(1-c_Al-c_V) + L0_HCP_Al_V/G_normalize - L0_HCP_Al_Ti/G_normalize -L0_HCP_Ti_V/G_normalize +
								   (1-2*c_Al-2*c_V)*L0_BCC_Al_Ti_V/G_normalize) + grad_cal*(-grad_cal/pow(c_Al,2) + (grad_cal +
								    grad_cv)/pow((1-c_Al-c_V),2) - 2*grad_cv*L0_BCC_Al_Ti_V/G_normalize) + grad_cv*((grad_cal +
								     grad_cv)/pow((1-c_Al-c_V),2) + (-2*grad_cal*grad_cv)*L0_BCC_Al_Ti_V/G_normalize) ;							
	delsq_dGBeta_dV = gradsq_cv*(1.0/c_V + 1/(1-c_Al-c_V) -2*L0_BCC_Ti_V/G_normalize - 2*c_Al*L0_BCC_Al_Ti_V/G_normalize) +
								 gradsq_cal*(1/(1-c_Al-c_V) + L0_BCC_Al_V/G_normalize - L0_BCC_Al_Ti/G_normalize - L0_BCC_Ti_V/G_normalize + 
								(1-2*c_Al-2*c_V)*L0_BCC_Al_Ti_V/G_normalize) + grad_cv*(-grad_cv/pow(c_V,2) + (grad_cal + 
								grad_cv)/pow((1-c_Al-c_V),2) - 2*grad_cal*L0_BCC_Al_Ti_V/G_normalize) + grad_cal*((grad_cal + 
								grad_cv)/pow((1-c_Al-c_V),2) + (-2*grad_cal*grad_cv)*L0_BCC_Al_Ti_V/G_normalize) ;
}

