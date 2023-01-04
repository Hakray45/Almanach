#include <stdlib.h>

#include <iostream>
#include <cmath>
using namespace std;

#define PI  3.141592653589793 // Pi number;
#define GM  398600441.8e+6 //- earth gravity
#define a_e  6378136 // the equatorial radius of the earth;
#define J  1082.62575e-6 //vector zonal harmonic decomposition of the geopotential in a series of spherical functions;
#define w_z  7.2921150e-5 // Rad/s is the angular velocity of the earth's rotation;

typedef struct AlmGlon {
    int     N_a; // Calendar number of days within a four-year period, starting from a leap year, to which the amendment applies Taus.
    double t_lam_a; // Moscow maternity time of passage of the j-m NCA of the ascending node of the orbit.
    double lam_a;
    double dI ;
    double w_a;
    double eps_a;
    double dT ;
    double dTT;

}AlmGlon;

AlmGlon constructA(int N_a,double t_lam_a,double del_T_a,double del_Ax_T_a,double lam_a,double w_a,double eps_a,double del_i_a){
    AlmGlon q;
    q.N_a=N_a;
    q.t_lam_a=t_lam_a;
    q.dT=del_T_a;
    q.dTT=del_Ax_T_a;
    q.lam_a=lam_a;
    q.w_a=w_a;
    q.eps_a=eps_a;
    q.dI=del_i_a;
    return q;
}

typedef struct Satellite {

     double x;
     double y;
     double z;

     double x_vel;
     double y_vel;
     double z_vel;

 }Satellite;







double* axis(double eps_a, double w_a, double T_dr, double ci) {


    double a = (pow((T_dr * T_dr) / (4 * PI * PI) * GM, 1. / 3.));
    double p = a * (1 - eps_a * eps_a);
    double Ak = 3. / 2. * J * pow((a_e / (p)), 2);
    double Bk = 2 - 5. / 2. * pow(sin(ci), 2);
    double Ck = pow(1 - eps_a * eps_a, 3. / 2.) / pow(1 + eps_a * cos(w_a * PI), 2);
    double Dk = pow(1 + eps_a * cos(w_a * PI), 3) / (1 - eps_a * eps_a);
    double T_osk = T_dr / (1 - Ak * (Bk * Ck + Dk));
    double del = a;
    double temp;
    do {

        del = a;
        Ak = 3. / 2. * J * pow((a_e / (p)), 2);
        Bk = 2 - 5. / 2. * pow(sin(ci), 2);
        Ck = pow(1 - eps_a * eps_a, 3. / 2.) / pow(1 + eps_a * cos(w_a * PI), 2);
        Dk = pow(1 + eps_a * cos(w_a * PI), 3) / (1 - eps_a * eps_a);
        T_osk = T_dr / (1 - Ak * (Bk * Ck + Dk));
        a = (pow(pow(T_osk / (2 * PI), 2) * GM, 1. / 3.));
        p = a * (1 - eps_a * eps_a);
        temp = a - del;



    } while (fabs(temp) > 10e-2);
    double arr[2];
    arr[0] = a;
    arr[1] = p;
    return arr;
}
double* error(double eps_a, double w, double a, double ci, double L) {
    double h = eps_a * sin(w);
    double l = eps_a * cos(w);
    double B = 1.5 * J * pow(a_e/a,2);
    double del_a = 2 * B * (1 - 1.5 * pow(sin(ci), 2)) * (l * cos(L) + h * sin(L)) + B * pow(sin(ci), 2) * (0.5 * h * sin(L) - 0.5 * l * cos(L) + cos(2 * L) + 3.5 * l * cos(3 * L) + 3.5 * h * sin(3 * L));
    double del_l = B * (1 - 1.5 * pow(sin(ci), 2)) * (cos(L) + 1.5 * l * cos(2 * L) + 1.5 * h * sin(2 * L)) - 0.25 * B * pow(sin(ci), 2) * (-cos(L) - (7. / 3.) * cos(3 * L) - 5 * h * sin(2 * L) - 8.5 * l * cos(4 * L) - 8.5 * h * sin(4 * L) + l * cos(2 * L)) + 0.5 * B * pow(cos(ci), 2) * h * sin(2 * L);
    double del_h = B * (1 - 1.5 * pow(sin(ci), 2)) * (sin(L) + 1.5 * l * sin(2 * L) - 1.5 * h * cos(2 * L)) - 0.25 * B * pow(sin(ci), 2) * (sin(L) - 7. / 3. * sin(3 * L) + 5 * l * sin(2 * L) - 8.5 * l * sin(4 * L) + 8.5 * h * cos(4 * L) + h * cos(2 * L)) + (-0.5 * B * pow(cos(ci), 2) * l * sin(2 * L));
    double del_lam = -B * cos(ci) * (3.5 * l * sin(L) - 2.5 * h * cos(L) - 0.5 * sin(2 * L) - 7. / 6. * l * sin(3 * L) + 7. / 6. * h * cos(3 * L));
    double del_i = 0.5 * B * sin(ci) * cos(ci) * (-l * cos(L) + h * sin(L) + cos(2 * L) + (7. / 3.) * l * cos(3 * L) + (7. / 3.) * h * sin(3 * L));
    double del_L = 2 * B * (1 - 1.5 * pow(sin(ci), 2)) * (1.75 * l * sin(L) - 1.75 * h * cos(L)) + 3 * B * pow(sin(ci), 2) * (-(7. / 24.) * h * cos(L) - (7. / 24.) * l * sin(L) - (49. / 72.) * h * cos(3 * L) + (49. / 72.) * l * sin(3 * L) + 0.25 * sin(2 * L)) + B * pow(cos(ci), 2) * (3.5 * l * sin(L) - 2.5 * h * cos(L) - 0.5 * sin(2 * L) + (7. / 6.) * l * sin(3 * L) + (7. / 6.) * h * cos(3 * L));


    double arr[8];
    arr[0] = h;
    arr[1] = l;
    arr[2] = del_a;
    arr[3] = del_l;
    arr[4] = del_h;
    arr[5] = del_lam;
    arr[6] = del_i;
    arr[7] = del_L;
    return arr;
}
double Kepler_equals(double new_L, double new_w, double new_eps) {

     double del;

    double E = new_L - new_w;
    double temp;


    do {
        temp = E;
        E = new_L - new_w + new_eps * sin(temp);

        del = E - temp;

    } while ( fabs(del) > 10e-9);

    return E;
}

void PredictionAlmanac(int N, int t_i, AlmGlon &Glon, Satellite &G){

    double T_sr = 43200; // ???
    //double N_a=Glon->N_a;
    //
    //double t_lam_a=Glon->t_lam_a;
    //double del_T_a = Glon->del_T_a;
    //double del_Ax_T_a = Glon->del_Ax_T_a;
    //double lam_a = Glon->lam_a;
    //double w_a = Glon->w_a;
    //double eps_a = Glon->eps_a;
    //double del_i_a = Glon->del_i_a;

    //--------------------------------------------
    //1.
    double del_t_pr = ((N - Glon.N_a) * 86400 + (t_i - Glon.t_lam_a));
    //--------------------------------------------
    //2.
    int W = (int)del_t_pr / (T_sr + Glon.dT);
    //--------------------------------------------
    //3.
    double i_sr = 63 * PI / 180;
    double ci = (i_sr / (180 * PI / 180) + Glon.dI) * PI;
    //--------------------------------------------
    //4.
    double T_dr = T_sr + Glon.dT + (2 * W + 1) * Glon.dTT;
    double n = 2 * PI / T_dr;
    //--------------------------------------------
    //5.
    double* test;
    test = axis(Glon.eps_a, Glon.w_a, T_dr, ci);
    double a = test[0];
    double p = test[1];

    //--------------------------------------------
    //6.
    double lam = Glon.lam_a * M_PI - (w_z + 1.5 * J * n * pow(a_e/p,2) * cos(ci)) * del_t_pr;
    double w = Glon.w_a * M_PI - 0.75 * J * n * pow(a_e/p,2) * (1 - 5 * pow(cos(ci), 2)) * del_t_pr;
    //--------------------------------------------
    //7.
    double E0 = -2 * atan(sqrt(1 - Glon.eps_a) / sqrt(1 + Glon.eps_a) * tan(w / 2));
    double L1 = w + E0 - Glon.eps_a * sin(E0);
    //--------------------------------------------
    //8.
    double L2 = L1 + n * (del_t_pr - (T_sr + Glon.dT) * W - Glon.dTT * W * W);

    //--------------------------------------------
    //9.
    double* del_error_L1 = error(Glon.eps_a, w, a, ci, L1);
    double temp[8], temp2[8];
    for (int i = 0;i < 8;i++) {
        temp[i] = del_error_L1[i];
    }

    double* del_error_L2 = error(Glon.eps_a, w, a, ci, L2);
    for (int i = 0;i < 8;i++) {
        temp2[i] = del_error_L2[i];
    }


    double new_h = temp[0] + temp2[4] - temp[4];
    double new_l = temp[1] + del_error_L2[3] - temp[3];
    double new_L = L2 + temp2[7] - temp[7];

    double new_w = atan(new_h / new_l);
    double new_lam = lam + temp2[5] - temp[5];
    double new_i = ci + temp2[6] - temp[6];
    double new_eps = sqrt(new_h * new_h + new_l * new_l);
    double new_a = a + temp2[2] * 10e+2 - temp[2] * 10e+2;
    //--------------------------------------------
   //10.

    double E = Kepler_equals(new_L, new_w, new_eps);
    //--------------------------------------------
   //11.
    double v = 2 * atan(sqrt((1 + new_eps) / (1 - new_eps)) * tan(E / 2));
    double u = v + new_w;
    //--------------------------------------------
   //12.
    p = new_a * (1 - new_eps * new_eps);

    double r = p / (1 + new_eps * cos(v));
    double x = r * (cos(new_lam) * cos(u) - sin(new_lam) * sin(u) * cos(new_i));
    double y = r * (sin(new_lam) * cos(u) + cos(new_lam) * sin(u) * cos(new_i));
    double z = r * (sin(u) * sin(ci));
    //--------------------------------------------
  //13.

    double v_r = sqrt(GM / p) * (new_eps * sin(v));
    double v_u = sqrt(GM / p) * ((1 + new_eps * cos(v)));
    double x_vel = v_r * (cos(new_lam) * cos(u) - sin(new_lam) * sin(u) * cos(new_i)) - v_u * (cos(new_lam) * sin(u) + sin(new_lam) * cos(u) * cos(new_i)) + w_z * y;
    double y_vel = v_r * (sin(new_lam) * cos(u) + cos(new_lam) * sin(u) * cos(new_i)) - v_u * (sin(new_lam) * sin(u) - cos(new_lam) * cos(u) * cos(new_i)) - w_z * x;
    double z_vel = v_r * sin(u) * sin(new_i) + v_u * cos(u) * sin(new_i);
    G.x=x;
    G.y=y;
    G.z=z;
    G.x_vel=x_vel;
    G.y_vel=y_vel;
    G.z_vel=z_vel;



}



void printA(Satellite G){
    cout<<"x:"<<G.x<<" "<<"y:"<<G.y<<" "<<"z:"<<G.z<<" "<<"Vx:"<<G.x_vel<<" "<<"Vy:"<<G.y_vel<<" "<<"Vz:"<<G.z_vel<<" "<<"\n";
   // cout <<"del  x:" << fabs(10697.116424527978-G.x)<<"\n"<<"del  y:" << fabs(21058.292414091822-G.y)<<"\n"
             //<<"del  z:" << fabs(-9635.6794316575106-G.z)<<"\n" <<"del  Vx:" << fabs(-0.68610081793104882
    //-G.x_vel)<<"\n" <<"del Vy:" << fabs(-1.1365486509759850-G.y_vel)<<"\n" <<"del Vz:" << fabs(-3.2499858708515017-G.z_vel)<<"\n";
}
