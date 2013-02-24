// snowflake.hpp
// Anisotropic solidification by 2D phase field method
// Questions/comments to kellet@rpi.edu (Trevor Keller)

#ifndef SNOWFLAKE_UPDATE
#define SNOWFLAKE_UPDATE
#include"MMSP.hpp"
#include<iomanip>
#include<cmath>
#include<ctime>
typedef MMSP::grid<1,MMSP::vector<double> > GRID1D;
typedef MMSP::grid<2,MMSP::vector<double> > GRID2D;
typedef MMSP::grid<3,MMSP::vector<double> > GRID3D;

void print_progress(const int step, const int steps, const int iterations);

namespace MMSP{

void generate(int dim, const char* filename)
{
  const int edge=4096;
  const double deltaX=0.025;
  const double undercooling=-0.5;
  if (dim==2)
  {
    MMSP::grid<2,vector<double> > grid(2,0,edge,0,edge);
    for (int d=0; d<dim; ++d) dx(grid,d)=deltaX;

    // Seed a circle of radius N*dx
    int R=5;
    for (int i=0; i<nodes(grid); ++i)
    {
      grid(i)[1]=undercooling; // Initial undercooling
      vector<int> x = position(grid,i);
      int r=sqrt(pow(x[0]-edge/2,2)+pow(x[1]-edge/2,2));
      if (r<=R) grid(i)[0]=1.;
      else grid(i)[0]=0.;
    }
    #ifdef MPI_VERSION
    output_mpi(grid,filename);
    #else
    output(grid,filename);
    #endif
  }
  else
  {
    std::cerr<<"Snowflake code is intended for 2D, only."<<std::endl;
    exit(1);
  }
}

template <int dim> void update(MMSP::grid<dim,vector<double> >& grid, int steps)
{
  int id=0;
  int np=1;
  #ifdef MPI_VERSION
  id=MPI::COMM_WORLD.Get_rank();
  np=MPI::COMM_WORLD.Get_size();
  #endif

  static int iterations=1;
  static MMSP::grid<2,double> grid_old(1,g0(grid,0),g1(grid,0),g0(grid,1),g1(grid,1));
  if (iterations==1)
    for (int i=0; i<nodes(grid_old); ++i) grid_old(i)=grid(i)[0];

  MMSP::grid<dim,vector<double> > update(grid);
  double      dt=5e-5;    // time-step
  double   theta=0.;      // angle relative to lab frame
  double       c=0.02;   // degree of anisotropy
  double       N=6.;      // symmetry
  double   alpha=0.015;   // gradient-energy coefficient
  double     tau=3e-4;    // time normalization constant
  double      k1=0.9;
  double      k2=20.;
  double   DiffT=2.25;    // thermal diffusivity
  double     CFL=tau/(2*alpha*alpha*(1./pow(dx(grid,0),2)+1./pow(dx(grid,1),2))); // Courant-Friedrich-Lewy condition on dt

  if (dt>0.5*CFL)
  {
    if (id==0) std::cout<<"dt="<<dt<<" is unstable; reduced to ";
    while (dt>0.5*CFL) dt*=3./4;
    if (id==0) std::cout<<dt<<"."<<std::endl;
  }

  std::cout.precision(2);

  int minus=0;
  int plus=0;
  #ifdef DEBUG
  if (id==0 && iterations==1) std::cout<<"dt is "<<dt<<", CFL condition is "<<dt/CFL<<'.'<<std::endl;
  double total_heat=0.;
  double total_phase=0.;
  #endif
  for (int step=0; step<steps; ++step)
  {
    #if (!defined DEBUG) && (!defined SILENT)
    if (id==0) print_progress(step, steps, iterations);
    #endif
    ghostswap(grid);
    MMSP::grid<dim,vector<double> > Dgradphi(grid);

    for (int i=0; i<nodes(grid); ++i)
    {
      vector<int> x=position(grid,i);

      // calculate grad(phi)
      vector<double> gradphi(dim,0.); // (0,0)
      for (int d=0; d<dim; ++d)
      {
        ++x[d];
        double right=grid(x)[0];
        --x[d];
        gradphi[d]=(right-grid(x)[0])/dx(grid,d);
      }
      double psi = theta + atan2(gradphi[1], gradphi[0]);
      double Phi = tan(N*psi/2.);
      double PhiSq = Phi*Phi;
      double beta = (1.-PhiSq)/(1.+PhiSq);
      double dBetadPsi = (-2.*N*Phi)/(1.+PhiSq);
      // Origin of this form for D is uncertain.
      Dgradphi(i)[0]=alpha*alpha*(1.+c*beta)*(   (1.+c*beta)*gradphi[0] - (c*dBetadPsi)*gradphi[1] );
      Dgradphi(i)[1]=alpha*alpha*(1.+c*beta)*( (c*dBetadPsi)*gradphi[0] +   (1.+c*beta)*gradphi[1] );
      #ifdef DEBUG
      if (step==0 && x[1]==g1(grid,1)/2+int(5/sqrt(2)*dx(grid,1)) && x[0]==g1(grid,0)/2-int(5/sqrt(2)*dx(grid,0)))
        std::cout<<"\n\n Iter.:  "<<iterations
                 <<"\n phase:  "<<grid(i)[0]
                 <<"\n alpha:  "<<alpha
                 <<"\n beta:   "<<beta
                 <<"\n psi:    "<<psi
                 <<"\n Phi:    "<<Phi
                 <<"\n dBdpsi: "<<dBetadPsi
                 <<"\n D∇φ:    ["<<std::setw(6)<<std::left<<Dgradphi(i)[0]
                 <<         ","<<std::setw(6)<<std::left<<Dgradphi(i)[1]<<"]";
      #endif
    }
    // Sync parallel grids
    ghostswap(Dgradphi);

    for (int i=0; i<nodes(grid); ++i)
    {
      vector<int> x = position(grid,i);

      // Update phase field
      double divDgradphi = 0.;
      for (int d=0; d<dim; ++d)
      {
        --x[d];
         double left=Dgradphi(x)[d];
         ++x[d];
         divDgradphi+=(Dgradphi(x)[d]-left)/dx(grid,d);
      }
      vector<double> old=grid(i);
      double m_phi=old[0]-0.5-(k1/M_PI)*atan(k2*old[1]);
      //double dphidt=divDgradphi+m_phi*phi_old*(1.-phi_old);
      // Semi-implicit scheme per Warren 2003
      if (m_phi>0)
      {
        ++plus;
        update(x)[0] = ((m_phi+tau/dt)*old[0]+divDgradphi)/(tau/dt+old[0]*m_phi);
      }
      else
      {
        ++minus;
        update(x)[0] = (old[0]*tau/dt+divDgradphi)/(tau/dt-(1.-old[0])*m_phi);
      }
      // Fully explicit forward-Euler discretization
      //update(x)[0] = grid(i)[0] + dt*dphidt/tau;

      // Update undercooling field
      double lapT=0;
      for (int d=0; d<dim; ++d)
      {
        ++x[d];
        double right=grid(x)[1];
        x[d]-=2;
        double left=grid(x)[1];
        ++x[d];
        lapT+=(right-(2*grid(x)[1])+left)/pow(dx(grid,d),2); // Laplacian
      }
      double dTdt = DiffT*lapT+(old[0]-grid_old(i))/dt;
      update(x)[1] = old[1] + dt*dTdt;
      #ifdef DEBUG
      if (step==0)
      {
        total_phase+=grid(i)[0];
        total_heat+=grid(i)[1];
      }
      #endif
    }
    #ifdef DEBUG
    if (step==0)
    {
      #ifdef MPI_VERSION
      double total_total_phase;
      MPI::COMM_WORLD.Allreduce(&total_phase, &total_total_phase, 1, MPI_DOUBLE, MPI_SUM);
      total_phase=total_total_phase;
      double total_total_heat;
      MPI::COMM_WORLD.Allreduce(&total_heat, &total_total_heat, 1, MPI_DOUBLE, MPI_SUM);
      total_heat=total_total_heat;
      #endif
      if (id==0) std::cout<<" Phase:  "<<total_phase
                          <<"\n Heat:   "<<total_heat
                          <<std::endl;
    }
    #endif
    for (int i=0; i<nodes(grid_old); ++i) grid_old(i)=grid(i)[0];
    swap(grid,update);
  }
  ghostswap(grid);
  ++iterations;
  #ifdef MDEBUG
    #ifdef MPI_VERSION
      int total_plus;
      int total_minus;
      MPI::COMM_WORLD.Allreduce(&plus, &total_plus, 1, MPI_INT, MPI_SUM);
      plus=total_plus;
      MPI::COMM_WORLD.Allreduce(&minus, &total_minus, 1, MPI_INT, MPI_SUM);
      minus=total_minus;
    #endif
    #ifndef SILENT
    if (id==0) std::cout<<"m_phi<0 "<<minus<<" times; m_phi≥0 "<<plus<<" times."<<std::endl;
    #endif
  #endif
}

} // namespace MMSP

void print_progress(const int step, const int steps, const int iterations)
{
  char* timestring;
  static unsigned long tstart;
  struct tm* timeinfo;

  if (step==0)
  {
    tstart = time(NULL);
    std::time_t rawtime;
    std::time( &rawtime );
    timeinfo = std::localtime( &rawtime );
    timestring = std::asctime(timeinfo);
    timestring[std::strlen(timestring)-1] = '\0';
    std::cout<<"Pass "<<std::setw(3)<<std::right<<iterations<<": "<<timestring<<" ["<<std::flush;
  }
  else if (step==steps-1)
  {
    unsigned long deltat = time(NULL)-tstart;
    //std::cout << "•] "<<timestring<<" ("<<std::setw(2)<<std::right<<deltat/3600<<"h:"<<std::setw(2)<<std::right<<(deltat%3600)/60<<"m:"<<std::setw(2)<<std::right<<deltat%60<<"s)"<<std::endl;
    std::cout << "•] "
              <<std::setw(2)<<std::right<<deltat/3600<<"h:"
              <<std::setw(2)<<std::right<<(deltat%3600)/60<<"m:"
              <<std::setw(2)<<std::right<<deltat%60<<"s"
              <<" (File "<<std::setw(5)<<std::right<<iterations*steps<<")."<<std::endl;
  }
  else if ((20 * step) % steps == 0) std::cout<<"• "<<std::flush;
}

#endif
