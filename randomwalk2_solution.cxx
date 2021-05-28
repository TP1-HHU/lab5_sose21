#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <random>
#include <iomanip>

using namespace std;
//------------------------------------------------------------------------------
struct particle{
    double x,y;
};

struct statistics{
    double msd,sx,sy;
};
//------------------------------------------------------------------------------
void init(particle* const p, const int N);
void push(particle* const p, const int N, default_random_engine& gen,
          uniform_real_distribution<double>& dist);
statistics stat(const particle* const p, const int N);
void write_to_file(const particle* const p, const int N, const string fname);
void write_to_file_binary(const particle* const p, const int N, const string fname);
string create_filename(const string basename, const int N);
void calc_density(int density[][512], particle* const p, const int Npart,
                  const int Nsteps);
void write_density(int density[][512], const int Nsteps,
                   const string fname);                  
//------------------------------------------------------------------------------
int main(void){
    random_device rd;
    default_random_engine gen(rd());
    uniform_real_distribution<double> dist(0.0,1.0);

    const int Npart  = 50000;  // Number of particles

    particle*  p  = new particle[Npart];
    statistics s;

    int density[512][512] = {0};

    ofstream ofstat("statistics");
    const string basename = "rwalk";

    const int Nsteps = 500;		// total # of steps
    const int Nfiles = 20;		// total # of output files (excl. initial file)
    int Nsubsteps    = Nsteps / Nfiles; // # steps between outputs

    init(p, Npart);

    for(int i = 0; i <= Nfiles; i++){
      s = stat(p,Npart);
      ofstat << i << "\t" << s.msd << "\t" << s.sx << "\t" << s.sy << endl;
      write_to_file_binary(p, Npart, create_filename(basename,i));
	  for(int j = 0; j < Nsubsteps; j++)
	         push(p, Npart, gen, dist);

    }
    ofstat.close();

    calc_density(density, p, Npart, Nsteps);
    write_density(density, Nsteps, "density");

    delete[] p;
    return 0;
}
//------------------------------------------------------------------------------
void calc_density(int density[][512], particle* const p, const int Npart,
                  const int Nsteps)
{

    const double min = -Nsteps;
    const double max =  Nsteps;
    const double h = 2*Nsteps / double(512);

    for(int c=0; c<Npart; c++){
        int j = int( (p[c].x - min) / h);
        int i = int( (max - p[c].y) / h);
        density[i][j]++;
    }

}
//------------------------------------------------------------------------------
void write_density(int density[][512], const int Nsteps,
                   const string fname){

   const int NX = 512;
   const double min = -Nsteps;
   const double h = 2*Nsteps/double(512);
   ofstream out(fname.c_str());

   for(int i=0; i<NX; i++){
     double y = min + (i+0.5)*h;
     for(int j=0; j<NX; j++){
        double x = min + (j+0.5)*h;
        out << x << "\t" << y << "\t" << density[i][j] << endl;
      }
      out << endl;
    }

   out.close();
}
//------------------------------------------------------------------------------
statistics stat(const particle* const p, const int N){
    statistics s;
    s.msd = 0; s.sx = 0; s.sy = 0;
    for(int i = 0; i < N; i++){
        double r2 = p[i].x*p[i].x + p[i].y*p[i].y;
        s.msd += r2;
        s.sx += p[i].x;
        s.sy += p[i].y;
    }
    s.msd /= N;
    return s;
}
//------------------------------------------------------------------------------
void push(particle* const p, const int N, default_random_engine& gen, 
          uniform_real_distribution<double>& dist){
    for(int i = 0; i < N; i++){
        double phi = 2*M_PI*dist(gen);
        double r = 1; //dist(gen);
        p[i].x += r*cos(phi);
	    p[i].y += r*sin(phi);
    }
}
//------------------------------------------------------------------------------
void init(particle* const p, const int N){
    for(int i = 0; i < N; i++){
	      p[i].x = 0;
	      p[i].y = 0;
    }
}
//------------------------------------------------------------------------------
string create_filename(const string basename, const int N){
    stringstream s;
    s.str("");
    s << basename << "_"  << N;
    return s.str();
}
//------------------------------------------------------------------------------
void write_to_file(const particle* const p, const int N, const string fname){
    ofstream out(fname.c_str());
    out << setprecision(17);
    for(int i = 0; i < N; i++)
	     out << p[i].x << "\t" << p[i].y << endl;
    out.close();
}
//------------------------------------------------------------------------------
void write_to_file_binary(const particle* const p, const int N, const string fname){
    ofstream out(fname.c_str(), ios::binary);
    
    // Einzeln
    //for(int i = 0; i < N; i++)
	//     out.write((char *)&p[i], sizeof(particle));
    
    // Oder direkt einfach alle
    out.write((char *)p, sizeof(particle)*N);
    
    out.close();
}
//------------------------------------------------------------------------------
