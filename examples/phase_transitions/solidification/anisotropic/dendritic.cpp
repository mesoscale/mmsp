// dendritic.cpp
// Anisotropic solidification by 2D phase field method
// Questions/comments to kellet@rpi.edu (Trevor Keller)

#ifndef DENDRITIC_UPDATE
#define DENDRITIC_UPDATE
#include<iomanip>
#include<cmath>
#include<ctime>
#include"MMSP.hpp"
#include"dendritic.hpp"

namespace MMSP{

void generate(int dim, const char* filename)
{
  const int edge=1024;
  const double deltaX=0.025;
  const double undercooling=-0.5;
  if (dim==2)
  {
    GRID2D initGrid(2,0,edge,0,edge);
    for (int d=0; d<dim; ++d) dx(initGrid,d)=deltaX;

    // Seed a circle of radius N*dx
    int R=5;
    for (int i=0; i<nodes(initGrid); ++i)
    {
      initGrid(i)[1]=undercooling; // Initial undercooling
      vector<int> x = position(initGrid,i);
      int r=sqrt(pow(x[0]-edge/2,2)+pow(x[1]-edge/2,2));
      if (r<=R) initGrid(i)[0]=1.;
      else initGrid(i)[0]=0.;
    }
    output(initGrid,filename);
  }
  else
  {
    std::cerr<<"Anisotropic solidification code is only implemented for 2D."<<std::endl;
    exit(1);
  }
}

template <int dim, typename T> void update(grid<dim,vector<T> >& refGrid, int steps)
{
  int id=0;
  int np=1;
  #ifdef MPI_VERSION
  id=MPI::COMM_WORLD.Get_rank();
  np=MPI::COMM_WORLD.Get_size();
  #endif

  static int iterations=1;
  static grid<2,T> oldGrid(1,g0(refGrid,0),g1(refGrid,0),g0(refGrid,1),g1(refGrid,1));
  if (iterations==1)
    for (int i=0; i<nodes(oldGrid); ++i) oldGrid(i)=refGrid(i)[0];

  grid<dim,vector<T> > newGrid(refGrid);
  double      dt=5e-5;    // time-step
  double   theta=0.;      // angle relative to lab frame
  double       c=0.02;   // degree of anisotropy
  double       N=6.;      // symmetry
  double   alpha=0.015;   // gradient-energy coefficient
  double     tau=3e-4;    // time normalization constant
  double      k1=0.9;
  double      k2=20.;
  double   DiffT=2.25;    // thermal diffusivity
  double     CFL=tau/(2*alpha*alpha*(1./pow(dx(refGrid,0),2)+1./pow(dx(refGrid,1),2))); // Courant-Friedrich-Lewy condition on dt

  if (dt>0.5*CFL)
  {
    if (id==0) std::cout<<"dt="<<dt<<" is unstable; reduced to ";
    while (dt>0.5*CFL) dt*=3./4;
    if (id==0) std::cout<<dt<<"."<<std::endl;
  }

  std::cout.precision(2);

  int minus=0;
  int plus=0;
  for (int step=0; step<steps; ++step)
  {
    if (id==0) print_progress(step, steps, iterations);
    ghostswap(refGrid);
    grid<dim,vector<T> > Dgradphi(refGrid);

    for (int i=0; i<nodes(refGrid); ++i)
    {
      vector<int> x=position(refGrid,i);

      // calculate grad(phi)
      vector<T> gradphi(dim,0.); // (0,0)
      for (int d=0; d<dim; ++d)
      {
        ++x[d];
        T right=refGrid(x)[0];
        --x[d];
        gradphi[d]=(right-refGrid(x)[0])/dx(refGrid,d);
      }
      T psi = theta + atan2(gradphi[1], gradphi[0]);
      T Phi = tan(N*psi/2.);
      T PhiSq = Phi*Phi;
      T beta = (1.-PhiSq)/(1.+PhiSq);
      T dBetadPsi = (-2.*N*Phi)/(1.+PhiSq);
      // Origin of this form for D is uncertain.
      Dgradphi(i)[0]=alpha*alpha*(1.+c*beta)*(   (1.+c*beta)*gradphi[0] - (c*dBetadPsi)*gradphi[1] );
      Dgradphi(i)[1]=alpha*alpha*(1.+c*beta)*( (c*dBetadPsi)*gradphi[0] +   (1.+c*beta)*gradphi[1] );
    }
    // Sync parallel grids
    ghostswap(Dgradphi);

    for (int i=0; i<nodes(refGrid); ++i)
    {
      vector<int> x = position(refGrid,i);

      // Update phase field
      T divDgradphi = 0.;
      for (int d=0; d<dim; ++d)
      {
        --x[d];
         T left=Dgradphi(x)[d];
         ++x[d];
         divDgradphi+=(Dgradphi(x)[d]-left)/dx(refGrid,d);
      }
      vector<T> old=refGrid(i);
      T m_phi=old[0]-0.5-(k1/M_PI)*atan(k2*old[1]);
      // Semi-implicit scheme per Warren 2003
      if (m_phi>0)
      {
        newGrid(x)[0] = ((m_phi+tau/dt)*old[0]+divDgradphi)/(tau/dt+old[0]*m_phi);
      }
      else
      {
        newGrid(x)[0] = (old[0]*tau/dt+divDgradphi)/(tau/dt-(1.-old[0])*m_phi);
      }
      // Fully explicit forward-Euler discretization
      //newGrid(x)[0] = refGrid(i)[0] + dt*dphidt/tau;

      // Update undercooling field
      T lapT=0;
      for (int d=0; d<dim; ++d)
      {
        ++x[d];
        T right=refGrid(x)[1];
        x[d]-=2;
        T left=refGrid(x)[1];
        ++x[d];
        lapT+=(right-(2*refGrid(x)[1])+left)/pow(dx(refGrid,d),2); // Laplacian
      }
      T dTdt = DiffT*lapT+(old[0]-oldGrid(i))/dt;
      newGrid(x)[1] = old[1] + dt*dTdt;
    }
    for (int i=0; i<nodes(oldGrid); ++i) oldGrid(i)=refGrid(i)[0];
    swap(refGrid,newGrid);
  }
  ghostswap(refGrid);
  ++iterations;
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

#include"MMSP.main.hpp"
