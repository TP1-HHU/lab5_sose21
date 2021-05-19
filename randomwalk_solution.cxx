#include <fstream>
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cmath>
#include <random>

using namespace std;
//------------------------------------------------------------------------------
struct particle{
    double x,y;
};
//------------------------------------------------------------------------------
struct statistics{
    double msd,sx,sy;
};
//------------------------------------------------------------------------------
void init(particle* const p, const int N);
void push(particle* const p, const int N, default_random_engine& gen, uniform_real_distribution<double>& dist);
statistics stat(const particle* const p, const int N);
void write_to_file(const particle* const p, const int N, const string fname);
string create_filename(const string basename, const int N);
//------------------------------------------------------------------------------
int main(void){
    random_device rd;
    default_random_engine gen(rd());
    uniform_real_distribution<double> dist(0.0,1.0);

    const int Npart  = 50000;  // Number of particles

    particle*  p  = new particle[Npart];
    statistics s;

    ofstream ofstat("statistics");
    const string basename = "rwalk";

    const int Nsteps = 500;		// total # of steps
    const int Nfiles = 20;		// total # of output files (excl. initial file)
    int Nsubsteps    = Nsteps / Nfiles; // # steps between outputs

    init(p, Npart);

    for(int i = 0; i <= Nfiles; i++){
      s = stat(p,Npart);
      ofstat << i << "\t" << s.msd << "\t" << s.sx << "\t" << s.sy << endl;
      write_to_file(p, Npart, create_filename(basename,i));
	  for(int j = 0; j < Nsubsteps; j++)
	         push(p, Npart, gen, dist);

    }

    ofstat.close();

    delete[] p;
    return 0;
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
void push(particle* const p, const int N, default_random_engine& gen, uniform_real_distribution<double>& dist){
    for(int i = 0; i < N; i++){
        double phi = 2*M_PI*dist(gen);
        double r = dist(gen);
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
    for(int i = 0; i < N; i++)
	     out << p[i].x << "\t" << p[i].y << endl;
    out.close();
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
