from math import  pow, sin, cos, atan2, tan, atan, sqrt
import numpy as np

PI = 3.141592653589793
GM = 398600441.8e+6         #- earth gravity
a_e = 6378136        # the equatorial radius of the earth;
J = 1082.62575e-6    #vector zonal harmonic decomposition of the geopotential in a series of spherical functions;
w_z = 7.2921150e-5   #Rad/s is the angular velocity of the earth's rotation;

def axis(eps_a:float, w_a:float, T_dr:float, ci:float)->int:

    a = (pow((T_dr * T_dr) / (4 * PI * PI) * GM, 1. / 3.))
    p = a * (1 - eps_a * eps_a)
    Ak = 3. / 2. * J * pow((a_e / (p)), 2)
    Bk = 2 - 5. / 2. * pow(sin(ci), 2)
    Ck = pow(1 - eps_a * eps_a, 3. / 2.) / pow(1 + eps_a * cos(w_a * PI), 2)
    Dk = pow(1 + eps_a * cos(w_a * PI), 3) / (1 - eps_a * eps_a)
    T_osk = T_dr / (1 - Ak * (Bk * Ck + Dk))
    i=0
    while True:
        delta = a
        Ak = 3. / 2. * J * pow((a_e / (p)), 2)
        T_osk = T_dr / (1 - Ak * (Bk * Ck + Dk))
        a = (pow(pow(T_osk / (2 * PI), 2) * GM, 1. / 3.))
        p = a * (1 - eps_a * eps_a)
        temp = a - delta
        if (abs(temp) > 10e-2) or i==100: break
        i+=1
    return a,p

def error(eps_a:float, w:float, a:float, ci:float, L:float)->list:
    h = eps_a * sin(w)
    l = eps_a * cos(w)
    B = 1.5 * J * pow(a_e / a, 2)

    del_a = 2 * B * (1 - 1.5 * pow(sin(ci), 2)) * (l * cos(L) + h * sin(L)) + B * pow(sin(ci), 2) * (
                0.5 * h * sin(L) - 0.5 * l * cos(L) + cos(2 * L) + 3.5 * l * cos(3 * L) + 3.5 * h * sin(3 * L))

    del_l = B * (1 - 1.5 * pow(sin(ci), 2)) * (cos(L) + 1.5 * l * cos(2 * L) + 1.5 * h * sin(2 * L)) - 0.25 * B * pow(
        sin(ci), 2) * (-cos(L) - (7. / 3.) * cos(3 * L) - 5 * h * sin(2 * L) - 8.5 * l * cos(4 * L) - 8.5 * h * sin(
        4 * L) + l * cos(2 * L)) + 0.5 * B * pow(cos(ci), 2) * h * sin(2 * L)

    del_h = B * (1 - 1.5 * pow(sin(ci), 2)) * (sin(L) + 1.5 * l * sin(2 * L) - 1.5 * h * cos(2 * L)) - 0.25 * B * pow(
        sin(ci), 2) * (sin(L) - 7. / 3. * sin(3 * L) + 5 * l * sin(2 * L) - 8.5 * l * sin(4 * L) + 8.5 * h * cos(
        4 * L) + h * cos(2 * L)) + (-0.5 * B * pow(cos(ci), 2) * l * sin(2 * L))

    del_lam = -B * cos(ci) * (
                3.5 * l * sin(L) - 2.5 * h * cos(L) - 0.5 * sin(2 * L) - 7. / 6. * l * sin(3 * L) + 7. / 6. * h * cos(
            3 * L))

    del_i = 0.5 * B * sin(ci) * cos(ci) * (
                -l * cos(L) + h * sin(L) + cos(2 * L) + (7. / 3.) * l * cos(3 * L) + (7. / 3.) * h * sin(3 * L))

    del_L = 2 * B * (1 - 1.5 * pow(sin(ci), 2)) * (1.75 * l * sin(L) - 1.75 * h * cos(L)) + 3 * B * pow(sin(ci), 2) * (
                -(7. / 24.) * h * cos(L) - (7. / 24.) * l * sin(L) - (49. / 72.) * h * cos(3 * L) + (
                    49. / 72.) * l * sin(3 * L) + 0.25 * sin(2 * L)) + B * pow(cos(ci), 2) * (
                        3.5 * l * sin(L) - 2.5 * h * cos(L) - 0.5 * sin(2 * L) + (7. / 6.) * l * sin(3 * L) + (
                            7. / 6.) * h * cos(3 * L))

    val=[h,l,del_a,del_l,del_h,del_lam,del_i,del_L]
    return val

def Kepler_equals(new_L:float, new_w:float, new_eps:float)->float:
    E = new_L - new_w
    i = 0

    while True:
        temp = E
        E = new_L - new_w + new_eps *sin(temp)
        delta = E - temp
        if (abs(delta)>10e-9) or i==100: break
        i+=1
    return E


def PredictionAlmanac(N:int, t_i:int,parameters:list)->list:
    T_sr = 43200
    #--------------------------------------------
    #1.
    del_t_pr = ((N - parameters[0]) * 86400 + (t_i - parameters[1]))
    # --------------------------------------------
    # 2
    W = int(del_t_pr / (T_sr + parameters[2]))
    # --------------------------------------------
    # 3.
    i_sr = 63 * PI / 180
    ci = (i_sr / (180 * PI / 180) + parameters[7]) * PI
    #--------------------------------------------
    #4.
    T_dr = T_sr + parameters[2] + (2 * W + 1) * parameters[3]
    n = 2 * PI / T_dr
    #--------------------------------------------
    #5
    a,p = axis(parameters[6],parameters[5],T_dr,ci)
    #--------------------------------------------
    #6.
    lam = parameters[4] * PI - (w_z + 1.5 * J * n * pow(a_e / p, 2) * cos(ci)) * del_t_pr
    w = parameters[5]* PI - 0.75 * J * n * pow(a_e / p, 2) * (1 - 5 * pow(cos(ci), 2)) * del_t_pr
    #--------------------------------------------
    #7.
    E0 = -2 * atan(sqrt(1 - parameters[6]) / sqrt(1 + parameters[6]) * tan(w / 2))
    L1 = w + E0 - parameters[6] * sin(E0)
    #--------------------------------------------
    #8.
    L2 = L1 + n * (del_t_pr - (T_sr + parameters[2]) * W - parameters[3] * W * W)
    #--------------------------------------------
    #9.
    del_error_L1 = error(parameters[6], w, a, ci, L1)
    del_error_L2 = error(parameters[6], w, a, ci, L2)

    new_h = del_error_L1[0] + del_error_L2[4] - del_error_L1[4]
    new_l = del_error_L1[1] + del_error_L2[3] - del_error_L1[3]
    new_L = L2 + del_error_L2[7] - del_error_L1[7]

    new_w = atan2(new_h,new_l)
    new_lam = lam + del_error_L2[5] - del_error_L1[5]
    new_i = ci + del_error_L2[6] - del_error_L1[6]
    new_eps = sqrt(new_h * new_h + new_l * new_l)
    new_a = a + del_error_L2[2] * 10e+2 - del_error_L1[2] * 10e+2
    #--------------------------------------------
    #10.
    E = Kepler_equals(new_L, new_w, new_eps)
    #--------------------------------------------
    #11.
    v = 2 * atan(sqrt((1 + new_eps) / (1 - new_eps)) * tan(E / 2))
    u = v + new_w
    #--------------------------------------------
    #12.
    p = new_a * (1 - new_eps * new_eps)

    r = p / (1 + new_eps * cos(v))
    x = (r * (cos(new_lam) * cos(u) - sin(new_lam) * sin(u) * cos(new_i)))/1e+3
    y = (r * (sin(new_lam) * cos(u) + cos(new_lam) * sin(u) * cos(new_i)))/1e+3
    z = (r * (sin(u) * sin(ci)))/1e+3
    #--------------------------------------------
    #13.
    v_r = sqrt(GM / p) * (new_eps * sin(v))
    v_u = sqrt(GM / p) * ((1 + new_eps * cos(v)))
    x_vel = (v_r * (cos(new_lam) * cos(u) - sin(new_lam) * sin(u) * cos(new_i)) - v_u * (cos(new_lam) * sin(u) + sin(new_lam) * cos(u) * cos(new_i)) + w_z * y)/1e+3
    y_vel = (v_r * (sin(new_lam) * cos(u) + cos(new_lam) * sin(u) * cos(new_i)) - v_u * (sin(new_lam) * sin(u) - cos(new_lam) * cos(u) * cos(new_i)) - w_z * x)/1e+3
    z_vel = (v_r * sin(u) * sin(new_i) + v_u * cos(u) * sin(new_i))/1e+3

    return [x,y,z,x_vel,y_vel,z_vel]


if __name__ == '__main__':

    N = 1453
    t_i = 51300

    N_a = 1452
    t_lam_a = 33571.625
    del_T_a = -2655.98046875
    del_Ax_T_a = 6.103515625e-05
    lam_a = -0.293967247009277
    w_a = 0.57867431640625
    eps_a = 0.000432968139648438
    del_i_a = 0.00987052917480469

    parameters = [N_a,t_lam_a,del_T_a,del_Ax_T_a,lam_a,w_a,eps_a,del_i_a]

    coord = PredictionAlmanac(N, t_i, parameters)
    print(coord)
