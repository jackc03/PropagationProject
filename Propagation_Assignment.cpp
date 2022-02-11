#include <cmath>
#include <array>
#include <iostream>
#include <fstream>
#include <iomanip>


int main(){
    double lambda = 8.0e-7;
    double pi     = 3.1415926535;
    double v      = 3.0e8;
    double tau    = 5.0e-15; 
    double xstart = -4.0e-5;
    double xfinal =  4.0e-5;
    const int size= 4000; //Nx is 4000
    double tstart = 0.0;
    double frontratio;
    double deltax, deltat, f, omega, k, xp, t1, tfinal;
    double t[size], x[size], E0[size], E1[size], E2[size];
    int Nx, Nt, n, i;


    deltax = (lambda) / 40;
    deltat = deltax / (v);
    f      = v / lambda;
    omega  = 2 * pi * f;
    k      = (2 * pi) / lambda;
    xp     = 0;
    Nx     = (xfinal - xstart) / deltax;
    Nt     = 1000;
    x[0]   = -4.0e-5;
    t[0]   = 0;
    //v[i] = c
    //vr[i] = c / n0
    //frontratio = v * v * deltat * deltat / (deltax * deltax);
    frontratio = ( ( std::pow(v, 2) * std::pow(deltat, 2) ) / std::pow(deltax, 2) );


    for(i = 1; i < Nx; i++){
        x[i] = x[i-1] + deltax;
        t[i] = t[i-1] + deltat;
    }

    for(i = 0; i < Nx; i++){
        E0[i] = std::exp(-1 * (std::pow( (t[0]-(x[i] / v) + (xp / v)), 2) / (std::pow(tau, 2)))) * cos((k * x[i]) - (omega * t[0]));
        E1[i] = std::exp(-1 * (std::pow( (t[1]-(x[i] / v) + (xp / v)), 2) / (std::pow(tau, 2)))) * cos((k * x[i]) - (omega * t[1]));
    }


    std::ofstream data1, data2;
    data1.open("EiNoRefractionCpp.dat");
    for(i = 1; i < Nx; i++){
        data1 << std::setprecision(20) << "\n" << x[i] << "\t" << E0[i] << "\t"<< E1[i];
    }
    data1.close();


    for( int n = 1; n <= Nt; n++){
        for(i = 1; i < Nx; i++){
            E2[i] = frontratio * (E1[i+1] - 2 * E1[i] + E1[i-1]) + 2 * E1[i] - E0[i];
        }
        for(i = 1; i < Nx; i++){
            E0[i] = E1[i];
            E1[i] = E2[i];
        }
    }


    data2.open("EfNoRefractionCpp.dat");
    for(i = 1; i < Nx; i++){
        data2 << std::setprecision(16) << "\n" << x[i] << std::setprecision(16) << "\t" << E1[i];
    }
    data2.close();

    system("pause");
    return 0;
}