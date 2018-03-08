//libraries needed: meep (header meep.hpp)
#include<meep.hpp>
#include<iostream>

using namespace meep;
using namespace std;

class Ceps : public material_function
{
public:
  Ceps(const vec &_center, double _rad, double _e):
    center(_center),rad(_rad),eps(_e){}
  ~Ceps(){};
  virtual bool has_mu(){return true;}
  virtual double chi1p1(field_type ft, const vec &r){
    if(abs(r-center)<rad && ft == E_stuff) return eps;
    return 1.0;
  }
  virtual void sigma_row(component c, double sigrow[3], const vec &r){
    sigrow[0]=0.0; sigrow[1]=0.0; sigrow[2]=0.0;
    if(abs(r-center)<rad){
      sigrow[component_direction(c)]=1.0;
      //###### THE LINES THAT BREAK EVERYTHING: BEGIN
      if(component_direction(c)==X) sigrow[1]+=0.3;
      if(component_direction(c)==Y) sigrow[0]+=0.3;
      //###### END
    }
  }
private:
  vec center;
  double rad, eps;
};


int main(int argc, char **argv){
  initialize mpi(argc,argv);
  const double amicron=50;
  const double fs=0.3;

  double radius_scatterer(0.5);

  double eps_inf(9.0685);
  double omega_d(2.1556/fs);
  double gamma_d(0.01836/fs);

  double pml_thickness=0.30;
  double cell_size=2*pml_thickness + 4.;
  vec center=vec(cell_size/2.,cell_size/2.);
    
  Ceps eps(center,radius_scatterer,eps_inf);
  double courant = 0.125;

  const grid_volume vol=vol2d(2*center.x(),2*center.y(),amicron);
  structure s(vol,eps,pml(pml_thickness));
  lorentzian_susceptibility drude_sus(omega_d,gamma_d,true);
  s.add_susceptibility(eps,E_stuff,drude_sus);  
  fields f(&s);

  const double freq_res=0.6/fs;
  const double tau=2.5*fs;
  const double sigma=tau/(2.0*sqrt(2.0*log(2.0)));
  const component src_comp=Hz;
  gaussian_src_time src(freq_res,sigma,0.0*sigma,10.0*sigma);
  const volume printvol(vec(pml_thickness,pml_thickness),
			center+center-vec(pml_thickness,pml_thickness));
  const volume srcvol(vec(0.,2.*center.y()-pml_thickness),
		      vec(2*center.x(),2.*center.y()-pml_thickness) );

  f.add_volume_source(src_comp,src,srcvol);

  int every_N=floor(0.5/(freq_res*f.dt));
  int counter(0);
  while (f.time() < f.last_source_time()+50*fs){
    f.step();
    if (counter%every_N==0){
      f.output_hdf5(Hz,printvol);
    }
    counter++;
  }
}
