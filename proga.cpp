#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#define N 300   // ДОЛЖНО ДЕЛИТСЯ НА 2
#define M N   // ДОЛЖНО ДЕЛИТСЯ НА 2
#define x_min -1.0
#define x_max 1.0
#define y_min -1.0
#define y_max 1.0
#define R 0.099
#define F 1.0
#define Al 1.5

#define K (N * M)
#define dx ((x_max - x_min)/(N-1))
#define dy ((y_max - y_min)/(M-1))
#define kv(x) ( (x) * (x) )

#define k1 (0.5 * a * ny/dy)
#define k2 (a * nx/dx)
#define k3 (b * nx/dx)
#define k4 (0.5 * a * ny/dx)
#define k5 (b * nx/dy)
#define k6 (a * nx/dx)
#define k7 (b * nx/dx)
#define k8 (0.5 * a * ny/dx)
#define k9 (0.5 * a * ny/dy)
#define k10 (b * nx/dy)

#define l1 (0.5 * a * nx/dy)
#define l2 (b * ny/dx)
#define l3 (a * ny/dy)
#define l4 (0.5 * a * nx/dx)
#define l5 (b * ny/dy)
#define l6 (b * ny/dx)
#define l7 (0.5 * a * nx/dx)
#define l8 (0.5 * a * nx/dy)
#define l9 (a * ny/dy)
#define l10 (b * ny/dy)

#define s1 (0.5 * a * ny/dy)
#define s2 (0.5 * a * nx/dx)
#define s3 (0.5 * a * nx/dx)
#define s4 (0.5 * a * ny/dy)

#define PI 3.1415926535

const double c49 = (4.0 / 9.0);
const double c13 = (1.0 / 3.0);
const double c23 = (2.0 / 3.0);
const double c29 = (2.0 / 9.0);
const double an = sqrt(2.0) / 2.0;
const double BA = (5.0 / 3.0);
const double CA = 0.3;
double AAk = 0.01;

const double a = c23;
const double b = -c29 + 1.0 / BA;

int mayak = 0;

int sign(double x)
{
    if (x > 0)
    {
        return 1;
    }
    else if (x < 0)
    {
        return -1;
    }
    else
    {
        return 0;
    }
}

double max(double x, double y)
{
    if (x >= y)
    {
        return x;
    }
    else
    {
        return y;
    }
}


double polar_angle(double x, double y)
{
    if (x < 0)
    {
        return atan(y / x) + 1.0 * PI;
    }
    else if (x > 0 && y >= 0)
    {
        return atan(y / x);
    }
    else if (x > 0 && y < 0)
    {
        return atan(y / x) + 2.0 * PI;
    }
    else if (y > 0 && x >= 0 && x <= 0)
    {
        return PI / 2.0;
    }
    else if (y < 0 && x >= 0 && x <= 0)
    {
        return  3.0 * PI / 2.0;
    }
    return 0.0;
}


double FF1(double gam)
{
    return gam + BA;
}

double dFF1(double gam)
{
    return 1.0;
}

double FF2(double gam)
{
    return 2.0 / gam;
}

double dFF2(double gam)
{
    return -2.0 / kv(gam);
}

double f1(const double& u1, const double& u2, const double& u3, const double& dxu1, const double& dxu2, const double& dxu3, //
    const double& dyu1, const double& dyu2, const double& dyu3, const double& ddxxu1, const double& ddxxu2, const double& ddxxu3, //
    const double& ddyyu1, const double& ddyyu2, const double& ddyyu3, const double& ddxyu1, const double& ddxyu2, const double& ddxyu3, bool dev_mode = false)
{
    double E = dxu1 + dyu2;
    double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
    if (c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)) < 0.001)
    {
        return 0.0;
    }
    double gam = E / E0;
    double dxE0 = 0.5 * (c49 * (2.0 * dxu1 * ddxxu1 + 2.0 * dyu2 * ddxyu2 - ddxxu1 * dyu2 - dxu1 * ddxyu2) + //
        c13 * (2.0 * dyu1 * ddxyu1 + 2.0 * dxu2 * ddxxu2 + 2.0 * ddxyu1 * dxu2 + 2.0 * dyu1 * ddxxu2 + 2.0 * dxu3 * ddxxu3 + //
            2.0 * dyu3 * ddxyu3)) / E0;
    double dyE0 = 0.5 * (c49 * (2.0 * dxu1 * ddxyu1 + 2.0 * dyu2 * ddyyu2 - ddxyu1 * dyu2 - dxu1 * ddyyu2) + //
        c13 * (2.0 * dyu1 * ddyyu1 + 2.0 * dxu2 * ddxyu2 + 2.0 * ddyyu1 * dxu2 + 2.0 * dyu1 * ddxyu2 + //
            2.0 * dxu3 * ddxyu3 + 2.0 * dyu3 * ddyyu3)) / E0;
    double dxE = ddxxu1 + ddxyu2;
    double dyE = ddxyu1 + ddyyu2;
    double dxgam = (dxE * E0 - dxE0 * E) / (kv(E0));
    double dygam = (dyE * E0 - dyE0 * E) / (kv(E0));
    double Exx = dxu1;
    double Eyy = dyu2;
    double Exy = 0.5 * (dyu1 + dxu2);
    double Exz = 0.5 * dxu3;
    double Eyz = 0.5 * dyu3;
    double dxExx = ddxxu1;
    double dyExx = ddxyu1;
    double dxEyy = ddxyu2;
    double dyEyy = ddyyu2;
    double dxExy = 0.5 * (ddxyu1 + ddxxu2);
    double dyExy = 0.5 * (ddyyu1 + ddxyu2);
    double dxExz = 0.5 * ddxxu3;
    double dyExz = 0.5 * ddxyu3;
    double dxEyz = 0.5 * ddxyu3;
    double dyEyz = 0.5 * ddyyu3;
    
    //cout << gam << endl;
    if (dev_mode)
    {
        cout << "ddyyu3  " << ddyyu3 << endl;
        cout << "ddxyu3  " << ddxyu3 << endl;
        cout << "ddxxu3  " << ddxxu3 << endl;
        cout << "ddyyu2  " << ddyyu2 << endl;
        cout << "dyu3  " << dyu3 << endl;
        cout << "dxu3  " << dxu3 << endl;
        cout << "dyE0  " << dyE0 << endl;
        cout << "dxE0  " << dxE0 << endl;
        cout << "dxgam   " << dxgam << endl;
        cout << "Exx   " << Exx  << endl;
        cout << "E   " << E << endl;
        cout << "1 = " << -1.5 * (CA / (BA - CA * CA)) * dxgam * (Exx - E / 3.0) << endl;
        cout << "2 = " << 1.5 * (CA / (BA - CA * CA)) * (CA - gam) * (dxExx - dxE/3.0) << endl;
        cout << "3 = " << 1.5 * (CA / (BA - CA * CA)) * dygam * Exy << endl;
        cout << "4 = " << 1.5 * (CA / (BA - CA * CA)) * (CA - gam) * dyExy << endl;
        cout << "5 = " << (CA / (BA - CA * CA)) * (CA * dxE - BA * dxE0) / BA << endl;
    }


    
    if (mayak == 1)
    {
        //для функций F1 и F2
            return -c23 * FF1(gam) * AAk * Al * pow(E0 , Al - 2) * (dxExx + dyExy) + c29 * FF1(gam) * AAk * Al * pow(E0 , Al - 2) * dxE - //
                c23 * dFF1(gam) * AAk * Al * pow(E0, Al - 2) * (Exx * dxgam + Exy * dygam - c13 * E * dxgam) - c23 * FF1(gam) * AAk * Al * //
                (Al -2) * pow(E0 , Al - 3) * (Exx * dxE0 + Exy * dyE0 - c13 * E * dxE0) - FF2(gam) * AAk * Al * pow(E0, Al - 2) * dxE - //
                dFF2(gam) * AAk * Al * pow(E0, Al - 2) * E * dxgam - FF2(gam) * AAk * Al * E * (Al - 2) * pow(E0, Al - 3)* dxE0;
    }
    else
    {
        //Для функций омега
        return -c23 * (CA / (BA - CA * CA)) * dxgam * (Exx - E / 3.0) + c23 * (CA / (BA - CA * CA)) * (CA - gam) * (dxExx - dxE / 3.0) - //
            c23 * (CA / (BA - CA * CA)) * dygam * Exy + c23 * (CA / (BA - CA * CA)) * (CA - gam) * dyExy + (CA / (BA - CA * CA)) * //
            (CA * dxE - BA * dxE0) / BA;
    }


}

double f2(const double& u1, const double& u2, const double& u3, const double& dxu1, const double& dxu2, const double& dxu3, //
    const double& dyu1, const double& dyu2, const double& dyu3, const double& ddxxu1, const double& ddxxu2, const double& ddxxu3, //
    const double& ddyyu1, const double& ddyyu2, const double& ddyyu3, const double& ddxyu1, const double& ddxyu2, const double& ddxyu3, bool dev_mode = false)
{
    double E = dxu1 + dyu2;
    double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
    if (c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)) < 0.001)
    {
        return 0.0;
    }
    double gam = E / E0;
    double dxE0 = 0.5 * (c49 * (2.0 * dxu1 * ddxxu1 + 2.0 * dyu2 * ddxyu2 - ddxxu1 * dyu2 - dxu1 * ddxyu2) + //
        c13 * (2.0 * dyu1 * ddxyu1 + 2.0 * dxu2 * ddxxu2 + 2.0 * ddxyu1 * dxu2 + 2.0 * dyu1 * ddxxu2 + 2.0 * dxu3 * ddxxu3 + //
            2.0 * dyu3 * ddxyu3)) / E0;
    double dyE0 = 0.5 * (c49 * (2.0 * dxu1 * ddxyu1 + 2.0 * dyu2 * ddyyu2 - ddxyu1 * dyu2 - dxu1 * ddyyu2) + //
        c13 * (2.0 * dyu1 * ddyyu1 + 2.0 * dxu2 * ddxyu2 + 2.0 * ddyyu1 * dxu2 + 2.0 * dyu1 * ddxyu2 + //
            2.0 * dxu3 * ddxyu3 + 2.0 * dyu3 * ddyyu3)) / E0;
    double dxE = ddxxu1 + ddxyu2;
    double dyE = ddxyu1 + ddyyu2;
    double dxgam = (dxE * E0 - dxE0 * E) / (kv(E0));
    double dygam = (dyE * E0 - dyE0 * E) / (kv(E0));
    double Exx = dxu1;
    double Eyy = dyu2;
    double Exy = 0.5 * (dyu1 + dxu2);
    double Exz = 0.5 * dxu3;
    double Eyz = 0.5 * dyu3;
    double dxExx = ddxxu1;
    double dyExx = ddxyu1;
    double dxEyy = ddxyu2;
    double dyEyy = ddyyu2;
    double dxExy = 0.5 * (ddxyu1 + ddxxu2);
    double dyExy = 0.5 * (ddyyu1 + ddxyu2);
    double dxExz = 0.5 * ddxxu3;
    double dyExz = 0.5 * ddxyu3;
    double dxEyz = 0.5 * ddxyu3;
    double dyEyz = 0.5 * ddyyu3;

    if (dev_mode)
    {
        cout << "2^^^^^^^^^^^^      " << endl;
        cout << "ddyyu3  " << ddyyu3 << endl;
        cout << "ddxyu3  " << ddxyu3 << endl;
        cout << "ddxxu3  " << ddxxu3 << endl;
        cout << "ddyyu2  " << ddyyu2 << endl;
        cout << "dyu3  " << dyu3 << endl;
        cout << "dxu3  " << dxu3 << endl;
        cout << "dyE  " << dyE << endl;
        cout << "dyE0  " << dyE0 << endl;
        cout << "dxgam   " << dxgam << endl;
        cout << "Exx   " << Exx << endl;
        cout << "E   " << E << endl;
        cout << "1 = " << -1.5 * (CA / (BA - CA * CA)) * dxgam * Exy << endl;
        cout << "2 = " << 1.5 * (CA / (BA - CA * CA)) * (CA - gam) * dxExy << endl;
        cout << "3 = " << 1.5 * (CA / (BA - CA * CA)) * dygam * (Eyy - E / 3.0) << endl;
        cout << "4 = " << 1.5 * (CA / (BA - CA * CA)) * (CA - gam) * (dyEyy - dyE / 3.0) << endl;
        cout << "5 = " << (CA / (BA - CA * CA)) *  (CA * dyE - BA * dyE0) / BA << endl;
    }

    if (mayak == 1)
    {
        // для F1 b и F2
        return -c23 * FF1(gam) * AAk * Al * pow(E0, Al - 2) * (dxExy + dyEyy) + c29 * FF1(gam) * AAk * Al * pow(E0, Al - 2) * dyE - //
            c23 * dFF1(gam) * AAk * Al * pow(E0, Al - 2) * (Exy * dxgam + Eyy * dygam - c13 * E * dygam) - c23 * FF1(gam) * AAk * Al * //
            (Al - 2) * pow(E0, Al - 3) * (Exy * dxE0 + Eyy * dyE0 - c13 * E * dyE0) - FF2(gam) * AAk * Al * pow(E0, Al - 2) * dyE - //
            dFF2(gam) * AAk * Al * pow(E0, Al - 2) * E * dygam - FF2(gam) * AAk * Al * E * (Al - 2) * pow(E0, Al - 3) * dyE0;
    }
    else
    {
        //Для функций омега
        return -c23 * (CA / (BA - CA * CA)) * dxgam * Exy + c23 * (CA / (BA - CA * CA)) * (CA - gam) * dxExy - //
            c23 * (CA / (BA - CA * CA)) * dygam * (Eyy - E / 3.0) + c23 * (CA / (BA - CA * CA)) * (CA - gam) * (dyEyy - dyE / 3.0) + //
            (CA / (BA - CA * CA)) * (CA * dyE - BA * dyE0) / BA;
    }

}

double f3(const double& u1, const double& u2, const double& u3, const double& dxu1, const double& dxu2, const double& dxu3, //
    const double& dyu1, const double& dyu2, const double& dyu3, const double& ddxxu1, const double& ddxxu2, const double& ddxxu3, //
    const double& ddyyu1, const double& ddyyu2, const double& ddyyu3, const double& ddxyu1, const double& ddxyu2, const double& ddxyu3, bool dev_mode = false)
{
    double E = dxu1 + dyu2;
    double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
    if (c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)) < 0.001)
    {
        return 0.0;
    }
    double gam = E / E0;
    double dxE0 = 0.5 * (c49 * (2.0 * dxu1 * ddxxu1 + 2.0 * dyu2 * ddxyu2 - ddxxu1 * dyu2 - dxu1 * ddxyu2) + //
        c13 * (2.0 * dyu1 * ddxyu1 + 2.0 * dxu2 * ddxxu2 + 2.0 * ddxyu1 * dxu2 + 2.0 * dyu1 * ddxxu2 + 2.0 * dxu3 * ddxxu3 + //
            2.0 * dyu3 * ddxyu3)) / E0;
    double dyE0 = 0.5 * (c49 * (2.0 * dxu1 * ddxyu1 + 2.0 * dyu2 * ddyyu2 - ddxyu1 * dyu2 - dxu1 * ddyyu2) + //
        c13 * (2.0 * dyu1 * ddyyu1 + 2.0 * dxu2 * ddxyu2 + 2.0 * ddyyu1 * dxu2 + 2.0 * dyu1 * ddxyu2 + //
            2.0 * dxu3 * ddxyu3 + 2.0 * dyu3 * ddyyu3)) / E0;
    double dxE = ddxxu1 + ddxyu2;
    double dyE = ddxyu1 + ddyyu2;
    double dxgam = (dxE * E0 - dxE0 * E) / (kv(E0));
    double dygam = (dyE * E0 - dyE0 * E) / (kv(E0));
    double Exx = dxu1;
    double Eyy = dyu2;
    double Exy = 0.5 * (dyu1 + dxu2);
    double Exz = 0.5 * dxu3;
    double Eyz = 0.5 * dyu3;
    double dxExx = ddxxu1;
    double dyExx = ddxyu1;
    double dxEyy = ddxyu2;
    double dyEyy = ddyyu2;
    double dxExy = 0.5 * (ddxyu1 + ddxxu2);
    double dyExy = 0.5 * (ddyyu1 + ddxyu2);
    double dxExz = 0.5 * ddxxu3;
    double dyExz = 0.5 * ddxyu3;
    double dxEyz = 0.5 * ddxyu3;
    double dyEyz = 0.5 * ddyyu3;

    if (mayak == 1)
    {
        // для F1 b и F2
        return -c23 * FF1(gam) * AAk * Al * pow(E0, Al - 2) * (dxExz + dyEyz) - c23 * dFF1(gam) * AAk * Al * pow(E0, Al - 2) * //
            (Exz * dxgam + Eyz * dygam) - c23 * FF1(gam) * AAk * Al * (Al - 2) * pow(E0, Al - 3) * (Exz * dxE0 + Eyz * dyE0);
    }
    else
    {
        //Для функций омега
        return -c23 * (CA / (BA - CA * CA)) * dxgam * Exz + c23 * (CA / (BA - CA * CA)) * (CA - gam) * dxExz - //
            c23 * (CA / (BA - CA * CA)) * dygam * Eyz + c23 * (CA / (BA - CA * CA)) * (CA - gam) * dyEyz;
    }

}

int main()
{

    // на нижнем слое
    double** u1 = new double* [N];
    double** u2 = new double* [N];
    double** u3 = new double* [N];

    // на верхнем слое
    double** uu1 = new double* [N];
    double** uu2 = new double* [N];
    double** uu3 = new double* [N];



    bool** IF = new bool* [N];
    for (int i = 0; i < N; i++)
    {
        u1[i] = new double[M];
        u2[i] = new double[M];
        u3[i] = new double[M];
        uu1[i] = new double[M];
        uu2[i] = new double[M];
        uu3[i] = new double[M];
        IF[i] = new bool[M];
    }

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < M; j++)
        {
            double x = x_min + i * dx;
            double y = y_min + j * dy;
            u1[i][j] = 0.0;
            u2[i][j] = 0.0;
            u3[i][j] = 0.0;
            uu1[i][j] = 0.0;
            uu2[i][j] = 0.0;
            uu3[i][j] = 0.0;
            if (sqrt( x*x + y*y ) <= R)
            {
                IF[i][j] = true;      // Если внутри круга
            }
            else
            {
                IF[i][j] = false;
            }
        }
    }

   /* for (int i = 0; i < N; i++)
    {
        u3[i][0] = 0.0;
        u3[i][M - 1] = 0.0;
        uu3[i][0] = 0.0;
        uu3[i][M - 1] = 0.0;
    }*/


    double A1, B1, C1, D1, E1, F1, G1, H1, I1;
    double A2, B2, C2, D2, E2, F2, G2, H2, I2;
    double A3, B3, C3, D3, E3, F3, G3, H3, I3;
    
    A1 = 2.0 * (a + b) / (dx * dx) + a / (dy * dy);
    B1 = F1 = (a + b) / (dx * dx);
    C1 = G1 = -(b + 0.5 * a) / (4.0 * dy * dx);
    D1 = H1 = 0.5 * a / (dy * dy);
    E1 = I1 = (b + 0.5 * a) / (4.0 * dy * dx);

    A2 = (a + b) * 2.0 / (dy * dy) + a / (dx * dx);
    B2 = F2 = 0.5 * a / (dx * dx);
    C2 = G2 = -(0.5 * a + b) / (4.0 * dy * dx);
    D2 = H2 = (a + b) / (dy * dy);
    E2 = I2 = (0.5 * a + b) / (4.0 * dy * dx);

    A3 = a / (dx * dx) + a / (dy * dy);
    B3 = F3 = 0.5 * a / (dx * dx);
    C3 = 0.0;
    D3 = H3 = 0.5 * a / (dy * dy);
    E3 = 0.0;
    G3 = 0.0;
    I3 = 0.0;



    double dxu1, dxu2, dxu3, dyu1, dyu2, dyu3, ddxxu1, ddxxu2, ddxxu3, ddyyu1, ddyyu2, ddyyu3, ddxyu1, ddxyu2, ddxyu3;
    double delta = 1.0;
    double** c;
    int step = 0;

    /*ifstream fin;
    fin.open("file_1.txt");
    double m1, m2, m3;
    for (int i = 0; i < N ; i++)
    {
        for (int j = 0; j < M ; j++)
        {
            fin >> m1 >> m2 >> m3;
            u1[i][j] = m1;
            u2[i][j] = m2;
            u3[i][j] = m3;
            uu1[i][j] = m1;
            uu2[i][j] = m2;
            uu3[i][j] = m3;
        }
    }
    fin.close();*/
    
    int St = 100000; //включение ненулевых граничных условий
    int St2 = 1000;

    ofstream foutf;
    foutf.open("xx_.txt");

    mayak = 0;
    cout << "START" << endl;
    while (step < 1100)//( (delta > 0.000001)||(step < 1500) )
    {
        step++;

        if (step > 5000000000)
        {
            mayak = 1; //влючаем финкции F1 и F2
        }

        if (step % 30 == 0)
        {
            cout << "delta = " << delta << "      step = " << step << "    may = " << mayak << "  AAK = " << AAk << endl;
        }

        delta = 0.0;
        // Задаём граничные условия (снос) на внешние грани

        // 1 (верх лево)
        double nx = 0.0;
        double ny = 1.0;
        for (int i = 1; i < N / 2; i++)
        {
            int j = M - 1;
            double dxu1 = (u1[i + 1][j] - u1[i - 1][j]) / 2.0 / dx;
            double dxu2 = (u2[i + 1][j] - u2[i - 1][j]) / 2.0 / dx;
            double dxu3 = (u3[i + 1][j] - u3[i - 1][j]) / 2.0 / dx;
            /*double dxu1 = (u1[i + 1][j] - u1[i][j]) / dx;
            double dxu2 = (u2[i + 1][j] - u2[i][j]) / dx;
            double dxu3 = (u3[i + 1][j] - u3[i][j]) / dx;*/
            double dyu1 = (u1[i][j] - u1[i][j - 1]) / dy;
            double dyu2 = (u2[i][j] - u2[i][j - 1]) / dy;
            double dyu3 = (u3[i][j] - u3[i][j - 1]) / dy;

            double E = dxu1 + dyu2;
            double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
            double gam = E / E0;
            double exx = dxu1 - E / 3.0;
            double eyy = dyu2 - E / 3.0;
            double ezz = -E / 3.0;
            double exy = 0.5 * (dyu1 + dxu2);
            double exz = 0.5 * dxu3;
            double eyz = 0.5 * dyu3;
            double t1, t2, t3;

            double Exx = dxu1;
            double Eyy = dyu2;
            double Exy = 0.5 * (dyu1 + dxu2);
            double Exz = 0.5 * dxu3;
            double Eyz = 0.5 * dyu3;

            if (step <= St)
            {
                t1 = t2 = t3 = 0.0;
            }
            else
            {
                if (mayak == 0)
                {
                    t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                }
                else
                {
                    t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                        FF2(gam) * AAk * Al * E * nx);
                    t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                        FF2(gam) * AAk * Al * E * ny);
                    t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                }
            }

            uu1[i][j] = (-t1 + a * u1[i][j - 1] / 2.0 / dy - a * (u2[i + 1][j] - u2[i - 1][j]) / 4.0 / dx) / (a / 2.0 / dy);
            uu2[i][j] = (-t2 + (a + b) * u2[i][j - 1] / dy - b * (u1[i + 1][j] - u1[i - 1][j]) / 2.0 / dx) / ((a + b) / dy);

            //double a1 = + k1 + k2 + k3;
            //double b1 = + k4 + k5;
            //double ff1 = -(( - k6 - k7) * u1[i - 1][j] - k8 * u2[i - 1][j] - k9 * u1[i][j - 1] - k10 * u2[i][j - 1]) - t1 - 0.1 * F;

            //double a2 = +l1 + l2;
            //double b2 = +l3 + l4 + l5;
            //double ff2 = -(-l6 * u1[i - 1][j] - l7 * u2[i - 1][j] - l8 * u1[i][j - 1] + (-l9 - l10) * u2[i][j - 1]) - t2;

            ///*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
            //uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
            //uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
            //uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;

            uu3[i][j] = (F - t3 + a * u3[i][j - 1] / dy / 2.0) / (a / 2.0 / dy);
        }


        // 8 (лево верх)
        nx = -1.0;
        ny = 0.0;
        for (int j = M - 2; j > M / 2; j--)
        {
            int i = 0;
            double dxu1 = (u1[i + 1][j] - u1[i][j]) / dx;
            double dxu2 = (u2[i + 1][j] - u2[i][j]) / dx;
            double dxu3 = (u3[i + 1][j] - u3[i][j]) / dx;
            double dyu1 = (u1[i][j + 1] - u1[i][j - 1]) / 2.0 / dy;
            double dyu2 = (u2[i][j + 1] - u2[i][j - 1]) / 2.0 / dy;
            double dyu3 = (u3[i][j + 1] - u3[i][j - 1]) / 2.0 / dy;

            double E = dxu1 + dyu2;
            double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
            double gam = E / E0;
            double exx = dxu1 - E / 3.0;
            double eyy = dyu2 - E / 3.0;
            double ezz = -E / 3.0;
            double exy = 0.5 * (dyu1 + dxu2);
            double exz = 0.5 * dxu3;
            double eyz = 0.5 * dyu3;
            double t1, t2, t3;
            if (step <= St)
            {
                t1 = t2 = t3 = 0.0;
            }
            else
            {
                if (mayak == 0)
                {
                    t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                }
                else
                {
                    t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                        FF2(gam) * AAk * Al * E * nx);
                    t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                        FF2(gam) * AAk * Al * E * ny);
                    t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                }
            }


            uu1[i][j] = (-t1 + (a + b) * u1[i + 1][j] / dx + b * (u2[i][j + 1] - u2[i][j - 1]) / (2.0 * dy)) / (a + b) * dx;
            uu2[i][j] = (-t2 + a * (u1[i][j + 1] - u1[i][j - 1]) / 4.0 / dy + a * u2[i + 1][j] / 2.0 / dx) / a * 2.0 * dx;



            //double a1 = +k1 - k2 - k3;
            //double b1 = -k4 + k5;
            //double ff1 = -((+k6 + k7) * u1[i + 1][j] + k8 * u2[i + 1][j] - k9 * u1[i][j - 1] - k10 * u2[i][j - 1]) - t1;

            //double a2 = +l1 - l2;
            //double b2 = +l3 - l4 + l5;
            //double ff2 = -(+l6 * u1[i + 1][j] + l7 * u2[i + 1][j] - l8 * u1[i][j - 1] + (-l9 - l10) * u2[i][j - 1]) - t2;

            ///*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
            //uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
            //uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
            //uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;

            uu3[i][j] = (-t3 + a * u3[i + 1][j] / 2 / dx) / a * 2.0 * dx;
        }


        // 2 (верх право)
        nx = 0.0;
        ny = 1.0;
        for (int i = N / 2; i < N - 1; i++)
        {
            int j = M - 1;
            double dxu1 = (u1[i + 1][j] - u1[i - 1][j]) / 2.0 / dx;
            double dxu2 = (u2[i + 1][j] - u2[i - 1][j]) / 2.0 / dx;
            double dxu3 = (u3[i + 1][j] - u3[i - 1][j]) / 2.0 / dx;
            /*double dxu1 = (u1[i][j] - u1[i - 1][j]) / dx;
            double dxu2 = (u2[i][j] - u2[i - 1][j])  / dx;
            double dxu3 = (u3[i][j] - u3[i - 1][j])  / dx;*/
            double dyu1 = (u1[i][j] - u1[i][j - 1]) / dy;
            double dyu2 = (u2[i][j] - u2[i][j - 1]) / dy;
            double dyu3 = (u3[i][j] - u3[i][j - 1]) / dy;

            double E = dxu1 + dyu2;
            double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
            double gam = E / E0;
            double exx = dxu1 - E / 3.0;
            double eyy = dyu2 - E / 3.0;
            double ezz = -E / 3.0;
            double exy = 0.5 * (dyu1 + dxu2);
            double exz = 0.5 * dxu3;
            double eyz = 0.5 * dyu3;
            double t1, t2, t3;
            if (step <= St)
            {
                t1 = t2 = t3 = 0.0;
            }
            else
            {
                if (mayak == 0)
                {
                    t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                }
                else
                {
                    t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                        FF2(gam) * AAk * Al * E * nx);
                    t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                        FF2(gam) * AAk * Al * E * ny);
                    t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                }
            }

            uu1[i][j] = (-t1 + a * u1[i][j - 1] / 2.0 / dy - a * (u2[i + 1][j] - u2[i - 1][j]) / 4.0 / dx) / (a / 2.0 / dy);
            uu2[i][j] = (-t2 + (a + b) * u2[i][j - 1] / dy - b * (u1[i + 1][j] - u1[i - 1][j]) / 2.0 / dx) / ((a + b) / dy);

            //double a1 = + k1 + k2 + k3;
            //double b1 = + k4 + k5;
            //double ff1 = -(( - k6 - k7) * u1[i - 1][j] - k8 * u2[i - 1][j] - k9 * u1[i][j - 1] - k10 * u2[i][j - 1]) - t1 - 0.1 * F;

            //double a2 = +l1 + l2;
            //double b2 = +l3 + l4 + l5;
            //double ff2 = -(-l6 * u1[i - 1][j] - l7 * u2[i - 1][j] - l8 * u1[i][j - 1] + (-l9 - l10) * u2[i][j - 1]) - t2;

            ///*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
            //uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
            //uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
            //uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;

            uu3[i][j] = (F - t3 + a * u3[i][j - 1] / dy / 2.0) / (a / 2.0 / dy);
        }


        //3 (право верх)
        nx = 1.0;
        ny = 0.0;
        for (int j = M - 2; j > M / 2; j--)
        {
            int i = N - 1;
            double dxu1 = (u1[i][j] - u1[i - 1][j]) / dx;
            double dxu2 = (u2[i][j] - u2[i - 1][j]) / dx;
            double dxu3 = (u3[i][j] - u3[i - 1][j]) / dx;
            double dyu1 = (u1[i][j + 1] - u1[i][j - 1]) / 2.0 / dy;
            double dyu2 = (u2[i][j + 1] - u2[i][j - 1]) / 2.0 / dy;
            double dyu3 = (u3[i][j + 1] - u3[i][j - 1]) / 2.0 / dy;

            double E = dxu1 + dyu2;
            double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
            double gam = E / E0;
            double exx = dxu1 - E / 3.0;
            double eyy = dyu2 - E / 3.0;
            double ezz = -E / 3.0;
            double exy = 0.5 * (dyu1 + dxu2);
            double exz = 0.5 * dxu3;
            double eyz = 0.5 * dyu3;
            double t1, t2, t3;
            if (step <= St)
            {
                t1 = t2 = t3 = 0.0;
            }
            else
            {
                if (mayak == 0)
                {
                    t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                }
                else
                {
                    t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                        FF2(gam) * AAk * Al * E * nx);
                    t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                        FF2(gam) * AAk * Al * E * ny);
                    t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                }
            }

            uu1[i][j] = (-t1 + (a + b) * u1[i - 1][j] / dx - b * (u2[i][j + 1] - u2[i][j - 1]) / 2.0 / dy) / (a + b) * dx;
            uu2[i][j] = (-t2 - a * (u1[i][j + 1] - u1[i][j - 1]) / 4.0 / dy + a * u2[i - 1][j] / 2.0 / dx) / a * 2.0 * dx;


            //double a1 = +k1 + k2 + k3;
            //double b1 = +k4 + k5;
            //double ff1 = -((-k6 - k7) * u1[i - 1][j] - k8 * u2[i - 1][j] - k9 * u1[i][j - 1] - k10 * u2[i][j - 1]) - t1;

            //double a2 = +l1 + l2;
            //double b2 = +l3 + l4 + l5;
            //double ff2 = -(-l6 * u1[i - 1][j] - l7 * u2[i - 1][j] - l8 * u1[i][j - 1] + (-l9 - l10) * u2[i][j - 1]) - t2;

            ///*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
            //uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
            //uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
            //uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;

            uu3[i][j] = (-t3 + a * u3[i - 1][j] / 2 / dx) / a * 2.0 * dx;
        }


        //4 (право низ)
        nx = 1.0;
        ny = 0.0;
        for (int j = M / 2; j > 0; j--)
        {
            int i = N - 1;
            double dxu1 = (u1[i][j] - u1[i - 1][j]) / dx;
            double dxu2 = (u2[i][j] - u2[i - 1][j]) / dx;
            double dxu3 = (u3[i][j] - u3[i - 1][j]) / dx;
            double dyu1 = (u1[i][j + 1] - u1[i][j - 1]) / 2.0 / dy;
            double dyu2 = (u2[i][j + 1] - u2[i][j - 1]) / 2.0 / dy;
            double dyu3 = (u3[i][j + 1] - u3[i][j - 1]) / 2.0 / dy;

            double E = dxu1 + dyu2;
            double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
            double gam = E / E0;
            double exx = dxu1 - E / 3.0;
            double eyy = dyu2 - E / 3.0;
            double ezz = -E / 3.0;
            double exy = 0.5 * (dyu1 + dxu2);
            double exz = 0.5 * dxu3;
            double eyz = 0.5 * dyu3;
            double t1, t2, t3;
            if (step <= St)
            {
                t1 = t2 = t3 = 0.0;
            }
            else
            {
                if (mayak == 0)
                {
                    t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                }
                else
                {
                    t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                        FF2(gam) * AAk * Al * E * nx);
                    t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                        FF2(gam) * AAk * Al * E * ny);
                    t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                }
            }


            uu1[i][j] = (-t1 + (a + b) * u1[i - 1][j] / dx - b * (u2[i][j + 1] - u2[i][j - 1]) / 2.0 / dy) / (a + b) * dx;
            uu2[i][j] = (-t2 - a * (u1[i][j + 1] - u1[i][j - 1]) / 4.0 / dy + a * u2[i - 1][j] / 2.0 / dx) / a * 2.0 * dx;

            //double a1 = -k1 + k2 + k3;
            //double b1 = +k4 - k5;
            //double ff1 = -((-k6 - k7) * u1[i - 1][j] - k8 * u2[i - 1][j] + k9 * u1[i][j + 1] + k10 * u2[i][j + 1]) - t1;

            //double a2 = -l1 + l2;
            //double b2 = -l3 + l4 - l5;
            //double ff2 = -(-l6 * u1[i - 1][j] - l7 * u2[i - 1][j] + l8 * u1[i][j + 1] + (+l9 + l10) * u2[i][j + 1]) - t2;

            ///*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
            //uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
            //uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
            //uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;

            uu3[i][j] = (-t3 + a * u3[i - 1][j] / 2 / dx) / a * 2.0 * dx;
        }


        // 5 (низ право)
        nx = 0.0;
        ny = -1.0;
        for (int i = N / 2; i < N - 1; i++)
        {
            int j = 0;
            double dxu1 = (u1[i + 1][j] - u1[i - 1][j]) / 2.0 / dx;
            double dxu2 = (u2[i + 1][j] - u2[i - 1][j]) / 2.0 / dx;
            double dxu3 = (u3[i + 1][j] - u3[i - 1][j]) / 2.0 / dx;
            /*double dxu1 = (u1[i][j] - u1[i - 1][j]) / dx;
            double dxu2 = (u2[i][j] - u2[i - 1][j])  / dx;
            double dxu3 = (u3[i][j] - u3[i - 1][j]) / dx;*/
            double dyu1 = (u1[i][j + 1] - u1[i][j]) / dy;
            double dyu2 = (u2[i][j + 1] - u2[i][j]) / dy;
            double dyu3 = (u3[i][j + 1] - u3[i][j]) / dy;

            double E = dxu1 + dyu2;
            double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
            double gam = E / E0;
            double exx = dxu1 - E / 3.0;
            double eyy = dyu2 - E / 3.0;
            double ezz = -E / 3.0;
            double exy = 0.5 * (dyu1 + dxu2);
            double exz = 0.5 * dxu3;
            double eyz = 0.5 * dyu3;
            double t1, t2, t3;
            if (step <= St)
            {
                t1 = t2 = t3 = 0.0;
            }
            else
            {
                if (mayak == 0)
                {
                    t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                }
                else
                {
                    t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                        FF2(gam) * AAk * Al * E * nx);
                    t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                        FF2(gam) * AAk * Al * E * ny);
                    t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                }
            }

            uu1[i][j] = (-t1 + a * u1[i][j + 1] / 2.0 / dy + a * (u2[i + 1][j] - u2[i - 1][j]) / 4.0 / dx) / (a / 2.0 / dy);
            uu2[i][j] = (-t2 + (a + b) * u2[i][j + 1] / dy + b * (u1[i + 1][j] - u1[i - 1][j]) / 2.0 / dx) / ((a + b) / dy);

            //double a1 = -k1 + k2 + k3;
            //double b1 = +k4 - k5;
            //double ff1 = -((-k6 - k7) * u1[i - 1][j] - k8 * u2[i - 1][j] + k9 * u1[i][j + 1] + k10 * u2[i][j + 1]) - t1;

            //double a2 = -l1 + l2;
            //double b2 = -l3 + l4 - l5;
            //double ff2 = -(-l6 * u1[i - 1][j] - l7 * u2[i - 1][j] + l8 * u1[i][j + 1] + (+l9 + l10) * u2[i][j + 1]) - t2;

            ///*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
            //uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
            //uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
            //uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;

            uu3[i][j] = (-F - t3 + a *  u3[i][j + 1] / dy / 2.0) / (a / 2.0 / dy);
        }


        // 6 (низ лево)
        nx = 0.0;
        ny = -1.0;
        for (int i = 1; i < N / 2; i++)
        {
            int j = 0;
            double dxu1 = (u1[i + 1][j] - u1[i - 1][j]) / 2.0 / dx;
            double dxu2 = (u2[i + 1][j] - u2[i - 1][j]) / 2.0 / dx;
            double dxu3 = (u3[i + 1][j] - u3[i - 1][j]) / 2.0 / dx;
            /*double dxu1 = (u1[i + 1][j] - u1[i][j]) / dx;
            double dxu2 = (u2[i + 1][j] - u2[i][j])  / dx;
            double dxu3 = (u3[i + 1][j] - u3[i][j])  / dx;*/
            double dyu1 = (u1[i][j + 1] - u1[i][j]) / dy;
            double dyu2 = (u2[i][j + 1] - u2[i][j]) / dy;
            double dyu3 = (u3[i][j + 1] - u3[i][j]) / dy;

            double E = dxu1 + dyu2;
            double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
            double gam = E / E0;
            double exx = dxu1 - E / 3.0;
            double eyy = dyu2 - E / 3.0;
            double ezz = -E / 3.0;
            double exy = 0.5 * (dyu1 + dxu2);
            double exz = 0.5 * dxu3;
            double eyz = 0.5 * dyu3;
            double t1, t2, t3;
            if (step <= St)
            {
                t1 = t2 = t3 = 0.0;
            }
            else
            {
                if (mayak == 0)
                {
                    t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                }
                else
                {
                    t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                        FF2(gam) * AAk * Al * E * nx);
                    t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                        FF2(gam) * AAk * Al * E * ny);
                    t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                }
            }

            uu1[i][j] = (-t1 + a * u1[i][j + 1] / 2.0 / dy + a * (u2[i + 1][j] - u2[i - 1][j]) / 4.0 / dx) / (a / 2.0 / dy);
            uu2[i][j] = (-t2 + (a + b) * u2[i][j + 1] / dy + b * (u1[i + 1][j] - u1[i - 1][j]) / 2.0 / dx) / ((a + b) / dy);
            //double a1 = -k1 - k2 - k3;
            //double b1 = -k4 - k5;
            //double ff1 = -((+k6 + k7) * u1[i + 1][j] + k8 * u2[i + 1][j] + k9 * u1[i][j + 1] + k10 * u2[i][j + 1]) - t1;

            //double a2 = -l1 - l2;
            //double b2 = -l3 - l4 - l5;
            //double ff2 = -(+l6 * u1[i + 1][j] + l7 * u2[i + 1][j] + l8 * u1[i][j + 1] + (+l9 + l10) * u2[i][j + 1]) - t2;

            ///*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
            //uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
            //uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
            //uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;

            uu3[i][j] = (-F - t3 + a * u3[i][j + 1] / dy / 2.0) / (a / 2.0 / dy);
        }


        // 7 (лево низ)
        nx = -1.0;
        ny = 0.0;
        for (int j = M / 2; j > 0; j--)
        {
            int i = 0;
            double dxu1 = (u1[i + 1][j] - u1[i][j]) / dx;
            double dxu2 = (u2[i + 1][j] - u2[i][j]) / dx;
            double dxu3 = (u3[i + 1][j] - u3[i][j]) / dx;
            double dyu1 = (u1[i][j + 1] - u1[i][j - 1]) / 2.0 / dy;
            double dyu2 = (u2[i][j + 1] - u2[i][j - 1]) / 2.0 / dy;
            double dyu3 = (u3[i][j + 1] - u3[i][j - 1]) / 2.0 / dy;

            double E = dxu1 + dyu2;
            double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
            double gam = E / E0;
            double exx = dxu1 - E / 3.0;
            double eyy = dyu2 - E / 3.0;
            double ezz = -E / 3.0;
            double exy = 0.5 * (dyu1 + dxu2);
            double exz = 0.5 * dxu3;
            double eyz = 0.5 * dyu3;
            double t1, t2, t3;
            if (step <= St)
            {
                t1 = t2 = t3 = 0.0;
            }
            else
            {
                if (mayak == 0)
                {
                    t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                }
                else
                {
                    t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                        FF2(gam) * AAk * Al * E * nx);
                    t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                        FF2(gam) * AAk * Al * E * ny);
                    t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                }
            }


            uu1[i][j] = (-t1 + (a + b) * u1[i + 1][j] / dx + b * (u2[i][j + 1] - u2[i][j - 1]) / (2.0 * dy)) / (a + b) * dx;
            uu2[i][j] = (-t2 + a * (u1[i][j + 1] - u1[i][j - 1]) / 4.0 / dy + a * u2[i + 1][j] / 2.0 / dx) / a * 2.0 * dx;

            //double a1 = -k1 - k2 - k3;
            //double b1 = -k4 - k5;
            //double ff1 = -((+k6 + k7) * u1[i + 1][j] + k8 * u2[i + 1][j] + k9 * u1[i][j + 1] + k10 * u2[i][j + 1]) - t1;

            //double a2 = -l1 - l2;
            //double b2 = -l3 - l4 - l5;
            //double ff2 = -(+l6 * u1[i + 1][j] + l7 * u2[i + 1][j] + l8 * u1[i][j + 1] + (+l9 + l10) * u2[i][j + 1]) - t2;

            ///*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
            //uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
            //uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
            //uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;

            uu3[i][j] = (-t3 + a * u3[i + 1][j] / 2 / dx) / a * 2.0 * dx;
        }


        /*uu1[0][0] = (u1[1][0] + u1[0][1])/2.0;
        uu2[0][0] = (u2[1][0] + u2[0][1]) / 2.0;
        uu3[0][0] = (u3[1][0] + u3[0][1]) / 2.0;

        uu1[N-1][0] = (u1[N-2][0] + u1[N-1][1]) / 2.0;
        uu2[N-1][0] = (u2[N-2][0] + u2[N-1][1]) / 2.0;
        uu3[N-1][0] = (u3[N-2][0] + u3[N-1][1]) / 2.0;

        uu1[0][M-1] = (u1[1][M-1] + u1[0][M-2]) / 2.0;
        uu2[0][M-1] = (u2[1][M-1] + u2[0][M-2]) / 2.0;
        uu3[0][M-1] = (u3[1][M-1] + u3[0][M-2]) / 2.0;

        uu1[N - 1][M - 1] = (u1[N - 2][M - 1] + u1[N - 1][M - 2]) / 2.0;
        uu2[N - 1][M - 1] = (u2[N - 2][M - 1] + u2[N - 1][M - 2]) / 2.0;
        uu3[N - 1][M - 1] = (u3[N - 2][M - 1] + u3[N - 1][M - 2]) / 2.0;*/
        /*uu1[0][0] = 0.0;
        uu2[0][0] = 0.0;
        uu3[0][0] = 0.0;

        uu1[N - 1][0] = 0.0;
        uu2[N - 1][0] = 0.0;
        uu3[N - 1][0] = 0.0;

        uu1[0][M - 1] = 0.0;
        uu2[0][M - 1] = 0.0;
        uu3[0][M - 1] = 0.0;

        uu1[N - 1][M - 1] = 0.0;
        uu2[N - 1][M - 1] = 0.0;
        uu3[N - 1][M - 1] = 0.0;*/

        if(true)// Углы
        {
        // правый верхний угол
        /*nx = an;
        ny = an;*/
        nx = 0.0;
        ny = 1.0;
        if (true)
        {
            int i = N - 1;
            int j = M - 1;
            double dxu1 = (u1[i][j] - u1[i - 1][j]) / dx;
            double dxu2 = (u2[i][j] - u2[i - 1][j]) / dx;
            double dxu3 = (u3[i][j] - u3[i - 1][j]) / dx;
            double dyu1 = (u1[i][j] - u1[i][j - 1]) / dy;
            double dyu2 = (u2[i][j] - u2[i][j - 1]) / dy;
            double dyu3 = (u3[i][j] - u3[i][j - 1]) / dy;

            double E = dxu1 + dyu2;
            double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
            double gam = E / E0;
            double exx = dxu1 - E / 3.0;
            double eyy = dyu2 - E / 3.0;
            double ezz = -E / 3.0;
            double exy = 0.5 * (dyu1 + dxu2);
            double exz = 0.5 * dxu3;
            double eyz = 0.5 * dyu3;
            double t1, t2, t3;
            if (step <= St)
            {
                t1 = t2 = t3 = 0.0;
            }
            else
            {
                if (mayak == 0)
                {
                    t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                }
                else
                {
                    t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                        FF2(gam) * AAk * Al * E * nx);
                    t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                        FF2(gam) * AAk * Al * E * ny);
                    t3 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny));
                }
            }


            /*uu1[i][j] = (-t1 + a * u1[i][j - 1] / 2.0 / dy - a * (u2[i][j] - u2[i - 1][j]) / 2.0 / dx) / (a / 2.0 / dy);
            uu2[i][j] = (-t2 + (a + b) * u2[i][j - 1] / dy - b * (u1[i][j] - u1[i - 1][j]) / dx) / ((a + b) / dy);
            uu3[i][j] = (F - t3 + a * u3[i][j - 1] / dy / 2.0) / (a / 2.0 / dy);*/

            //double E = dxu1 + dyu2;
            //double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
            //double gam = E / E0;
            //double exx = dxu1 - E / 3.0;
            //double eyy = dyu2 - E / 3.0;
            //double ezz = -E / 3.0;
            //double exy = 0.5 * (dyu1 + dxu2);
            //double exz = 0.5 * dxu3;
            //double eyz = 0.5 * dyu3;
            //double t1, t2, t3;
            //if (step <= St)
            //{
            //    t1 = t2 = t3 = 0.0;
            //}
            //else
            //{
            //    if (mayak == 0)
            //    {
            //        t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
            //        t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
            //        t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
            //    }
            //    else
            //    {
            //        t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
            //            FF2(gam) * AAk * Al * E * nx);
            //        t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
            //            FF2(gam) * AAk * Al * E * ny);
            //        t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
            //    }
            //}



            double a1 = +k1 + k2 + k3;
            double b1 = +k4 + k5;
            double ff1 = -((-k6 - k7) * u1[i - 1][j] - k8 * u2[i - 1][j] - k9 * u1[i][j - 1] - k10 * u2[i][j - 1]) - t1;

            double a2 = +l1 + l2;
            double b2 = +l3 + l4 + l5;
            double ff2 = -(-l6 * u1[i - 1][j] - l7 * u2[i - 1][j] - l8 * u1[i][j - 1] + (-l9 - l10) * u2[i][j - 1]) - t2;

            ///*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
            //uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
            uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
            uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;

            uu3[i][j] = (-t3 - (-s3 * u3[i - 1][j] - s4 * u3[i][j - 1])) / (+s1 + s2);
        }


        // левый верхний угол 
        nx = -an;
        ny = an;
        /*nx = 0.0;
        ny = 1.0;*/
        if (true)
        {
            int i = 0;
            int j = M - 1;
            double dxu1 = (u1[i + 1][j] - u1[i][j]) / dx;
            double dxu2 = (u2[i + 1][j] - u2[i][j]) / dx;
            double dxu3 = (u3[i + 1][j] - u3[i][j]) / dx;
            double dyu1 = (u1[i][j] - u1[i][j - 1]) / dy;
            double dyu2 = (u2[i][j] - u2[i][j - 1]) / dy;
            double dyu3 = (u3[i][j] - u3[i][j - 1]) / dy;

            double E = dxu1 + dyu2;
            double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
            double gam = E / E0;
            double exx = dxu1 - E / 3.0;
            double eyy = dyu2 - E / 3.0;
            double ezz = -E / 3.0;
            double exy = 0.5 * (dyu1 + dxu2);
            double exz = 0.5 * dxu3;
            double eyz = 0.5 * dyu3;
            double t1, t2, t3;
            if (step <= St)
            {
                t1 = t2 = t3 = 0.0;
            }
            else
            {
                if (mayak == 0)
                {
                    t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                }
                else
                {
                    t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                        FF2(gam) * AAk * Al * E * nx);
                    t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                        FF2(gam) * AAk * Al * E * ny);
                    t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                }
            }


            double a1 = +k1 - k2 - k3;
            double b1 = -k4 + k5;
            double ff1 = -((+k6 + k7) * u1[i + 1][j] + k8 * u2[i + 1][j] - k9 * u1[i][j - 1] - k10 * u2[i][j - 1]) - t1;

            double a2 = +l1 - l2;
            double b2 = +l3 - l4 + l5;
            double ff2 = -(+l6 * u1[i + 1][j] + l7 * u2[i + 1][j] - l8 * u1[i][j - 1] + (-l9 - l10) * u2[i][j - 1]) - t2;

            ///*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
            //uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
            uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
            uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;
            uu3[i][j] = (-t3 - (+s3 * u3[i + 1][j] - s4 * u3[i][j - 1])) / (+s1 - s2);
        }

        // правый нижний угол
        nx = an;
        ny = -an;
        /*nx = 0.0;
        ny = -1.0;*/
        if (true)
        {
            int i = N - 1;
            int j = 0;
            double dxu1 = (u1[i][j] - u1[i - 1][j]) / dx;
            double dxu2 = (u2[i][j] - u2[i - 1][j]) / dx;
            double dxu3 = (u3[i][j] - u3[i - 1][j]) / dx;
            double dyu1 = (u1[i][j + 1] - u1[i][j]) / dy;
            double dyu2 = (u2[i][j + 1] - u2[i][j]) / dy;
            double dyu3 = (u3[i][j + 1] - u3[i][j]) / dy;

            double E = dxu1 + dyu2;
            double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
            double gam = E / E0;
            double exx = dxu1 - E / 3.0;
            double eyy = dyu2 - E / 3.0;
            double ezz = -E / 3.0;
            double exy = 0.5 * (dyu1 + dxu2);
            double exz = 0.5 * dxu3;
            double eyz = 0.5 * dyu3;
            double t1, t2, t3;
            if (step <= St)
            {
                t1 = t2 = t3 = 0.0;
            }
            else
            {
                if (mayak == 0)
                {
                    t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                }
                else
                {
                    t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                        FF2(gam) * AAk * Al * E * nx);
                    t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                        FF2(gam) * AAk * Al * E * ny);
                    t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                }
            }


            double a1 = -k1 + k2 + k3;
            double b1 = +k4 - k5;
            double ff1 = -((-k6 - k7) * u1[i - 1][j] - k8 * u2[i - 1][j] + k9 * u1[i][j + 1] + k10 * u2[i][j + 1]) - t1;

            double a2 = -l1 + l2;
            double b2 = -l3 + l4 - l5;
            double ff2 = -(-l6 * u1[i - 1][j] - l7 * u2[i - 1][j] + l8 * u1[i][j + 1] + (+l9 + l10) * u2[i][j + 1]) - t2;

            /*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
            uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
            uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
            uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;

            uu3[i][j] = (-t3 - (-s3 * u3[i - 1][j] + s4 * u3[i][j + 1])) / (-s1 + s2);
        }

        // левый нижний угол
        nx = -an;
        ny = -an;
        /*nx = 0.0;
        ny = -1.0;*/
        if (true)
        {
            int i = 0;
            int j = 0;
            double dxu1 = (u1[i + 1][j] - u1[i][j]) / dx;
            double dxu2 = (u2[i + 1][j] - u2[i][j]) / dx;
            double dxu3 = (u3[i + 1][j] - u3[i][j]) / dx;
            double dyu1 = (u1[i][j + 1] - u1[i][j]) / dy;
            double dyu2 = (u2[i][j + 1] - u2[i][j]) / dy;
            double dyu3 = (u3[i][j + 1] - u3[i][j]) / dy;

            double E = dxu1 + dyu2;
            double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
            double gam = E / E0;
            double exx = dxu1 - E / 3.0;
            double eyy = dyu2 - E / 3.0;
            double ezz = -E / 3.0;
            double exy = 0.5 * (dyu1 + dxu2);
            double exz = 0.5 * dxu3;
            double eyz = 0.5 * dyu3;
            double t1, t2, t3;
            if (step <= St)
            {
                t1 = t2 = t3 = 0.0;
            }
            else
            {
                if (mayak == 0)
                {
                    t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                }
                else
                {
                    t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                        FF2(gam) * AAk * Al * E * nx);
                    t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                        FF2(gam) * AAk * Al * E * ny);
                    t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                }
            }


            double a1 = -k1 - k2 - k3;
            double b1 = -k4 - k5;
            double ff1 = -((+k6 + k7) * u1[i + 1][j] + k8 * u2[i + 1][j] + k9 * u1[i][j + 1] + k10 * u2[i][j + 1]) - t1;

            double a2 = -l1 - l2;
            double b2 = -l3 - l4 - l5;
            double ff2 = -(+l6 * u1[i + 1][j] + l7 * u2[i + 1][j] + l8 * u1[i][j + 1] + (+l9 + l10) * u2[i][j + 1]) - t2;

            /*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
            uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
            uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
            uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;

            uu3[i][j] = (-t3 - (+s3 * u3[i + 1][j] + s4 * u3[i][j + 1])) / (-s1 - s2);
        }
    }



        // Задаём граничные условия в кругу
        double p1, p2, p3, p4, p5, p6, p7, p8;
        double S;
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < M; j++)
            {
                if (IF[i][j] == true)
                {
                    double x = x_min + i * dx;
                    double y = y_min + j * dy;
                    double r = sqrt(x *x + y * y);
                    if (r < 0.3 * R)
                    {
                        continue;
                    }
                    nx = x / r;
                    ny = y / r;
                    int kk = (int)(M / 2);
                    int kk2 = (int)(N / 2);
                    //if ((j == kk + 1 || j == kk + 2 || j == kk || j == kk + 3 || j == kk - 1) &&(x < 0) )
                    //{
                    //    double dxu1 = (u1[i][j] - u1[i - 1][j]) / dx;
                    //    double dxu2 = (u2[i][j] - u2[i - 1][j]) / dx;
                    //    double dxu3 = (u3[i][j] - u3[i - 1][j]) / dx;
                    //    double dyu1 = (u1[i][j + 1] - u1[i][j - 1]) / 2.0 / dy;
                    //    double dyu2 = (u2[i][j + 1] - u2[i][j - 1]) / 2.0 / dy;
                    //    double dyu3 = (u3[i][j + 1] - u3[i][j - 1]) / 2.0 / dy;

                    //    double E = dxu1 + dyu2;
                    //    double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
                    //    double gam = E / E0;
                    //    double exx = dxu1 - E / 3.0;
                    //    double eyy = dyu2 - E / 3.0;
                    //    double ezz = -E / 3.0;
                    //    double exy = 0.5 * (dyu1 + dxu2);
                    //    double exz = 0.5 * dxu3;
                    //    double eyz = 0.5 * dyu3;
                    //    double t1, t2, t3;
                    //    if (step <= St)
                    //    {
                    //        t1 = t2 = t3 = 0.0;
                    //    }
                    //    else
                    //    {
                    //        if (mayak == 0)
                    //        {
                    //            t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    //            t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    //            t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                    //        }
                    //        else
                    //        {
                    //            t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                    //                FF2(gam) * AAk * Al * E * nx);
                    //            t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                    //                FF2(gam) * AAk * Al * E * ny);
                    //            t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                    //        }
                    //    }

                    //    uu1[i][j] = (-t1 + (a + b) * u1[i - 1][j] / dx - b * (u2[i][j + 1] - u2[i][j - 1]) / 2.0 / dy) / (a + b) * dx;
                    //    uu2[i][j] = (-t2 - a * (u1[i][j + 1] - u1[i][j - 1]) / 4.0 / dy + a * u2[i - 1][j] / 2.0 / dx) / a * 2.0 * dx;

                    //    uu3[i][j] = (-t3 - (-s3 * u3[i - 1][j] - s4 * u3[i][j - 1])) / (+s1 + s2);
                    //}

                    //else if ((j == kk + 1 || j == kk + 2 || j == kk || j == kk + 3 || j == kk - 1) && (x > 0))
                    //{
                    //    double dxu1 = (u1[i + 1][j] - u1[i][j]) / dx;
                    //    double dxu2 = (u2[i + 1][j] - u2[i][j]) / dx;
                    //    double dxu3 = (u3[i + 1][j] - u3[i][j]) / dx;
                    //    double dyu1 = (u1[i][j + 1] - u1[i][j - 1]) / 2.0 / dy;
                    //    double dyu2 = (u2[i][j + 1] - u2[i][j - 1]) / 2.0 / dy;
                    //    double dyu3 = (u3[i][j + 1] - u3[i][j - 1]) / 2.0 / dy;

                    //    double E = dxu1 + dyu2;
                    //    double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
                    //    double gam = E / E0;
                    //    double exx = dxu1 - E / 3.0;
                    //    double eyy = dyu2 - E / 3.0;
                    //    double ezz = -E / 3.0;
                    //    double exy = 0.5 * (dyu1 + dxu2);
                    //    double exz = 0.5 * dxu3;
                    //    double eyz = 0.5 * dyu3;
                    //    double t1, t2, t3;
                    //    if (step <= St)
                    //    {
                    //        t1 = t2 = t3 = 0.0;
                    //    }
                    //    else
                    //    {
                    //        if (mayak == 0)
                    //        {
                    //            t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    //            t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    //            t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                    //        }
                    //        else
                    //        {
                    //            t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                    //                FF2(gam) * AAk * Al * E * nx);
                    //            t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                    //                FF2(gam) * AAk * Al * E * ny);
                    //            t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                    //        }
                    //    }


                    //    uu1[i][j] = (-t1 + (a + b) * u1[i + 1][j] / dx + b * (u2[i][j + 1] - u2[i][j - 1]) / (2.0 * dy)) / (a + b) * dx;
                    //    uu2[i][j] = (-t2 + a * (u1[i][j + 1] - u1[i][j - 1]) / 4.0 / dy + a * u2[i + 1][j] / 2.0 / dx) / a * 2.0 * dx;

                    //    uu3[i][j] = (-t3 - (+s3 * u3[i + 1][j] - s4 * u3[i][j - 1])) / (+s1 - s2);
                    //}

                    //else if ((i == kk2 + 1 || i == kk2 + 2 || i == kk2  || i == kk2 + 3 || i == kk2 - 1) && (y > 0))
                    //{
                    //    double dxu1 = (u1[i + 1][j] - u1[i - 1][j]) / 2.0 / dx;
                    //    double dxu2 = (u2[i + 1][j] - u2[i - 1][j]) / 2.0 / dx;
                    //    double dxu3 = (u3[i + 1][j] - u3[i - 1][j]) / 2.0 / dx;
                    //    double dyu1 = (u1[i][j + 1] - u1[i][j]) / dy;
                    //    double dyu2 = (u2[i][j + 1] - u2[i][j]) / dy;
                    //    double dyu3 = (u3[i][j + 1] - u3[i][j]) / dy;

                    //    double E = dxu1 + dyu2;
                    //    double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
                    //    double gam = E / E0;
                    //    double exx = dxu1 - E / 3.0;
                    //    double eyy = dyu2 - E / 3.0;
                    //    double ezz = -E / 3.0;
                    //    double exy = 0.5 * (dyu1 + dxu2);
                    //    double exz = 0.5 * dxu3;
                    //    double eyz = 0.5 * dyu3;
                    //    double t1, t2, t3;
                    //    if (step <= St)
                    //    {
                    //        t1 = t2 = t3 = 0.0;
                    //    }
                    //    else
                    //    {
                    //        if (mayak == 0)
                    //        {
                    //            t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    //            t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    //            t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                    //        }
                    //        else
                    //        {
                    //            t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                    //                FF2(gam) * AAk * Al * E * nx);
                    //            t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                    //                FF2(gam) * AAk * Al * E * ny);
                    //            t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                    //        }
                    //    }

                    //    uu1[i][j] = (-t1 + a * u1[i][j + 1] / 2.0 / dy + a * (u2[i + 1][j] - u2[i - 1][j]) / 4.0 / dx) / (a / 2.0 / dy);
                    //    uu2[i][j] = (-t2 + (a + b) * u2[i][j + 1] / dy + b * (u1[i + 1][j] - u1[i - 1][j]) / 2.0 / dx) / ((a + b) / dy);

                    //    uu3[i][j] = (-F - t3 - (-s3 * u3[i - 1][j] + s4 * u3[i][j + 1])) / (-s1 + s2);
                    //}

                    //else if ((i == kk2 + 1 ||  i == kk2 + 1 || i == kk2 + 2 || i == kk2 || i == kk2 + 3 || i == kk2 - 1) && (y < 0))
                    //{
                    //double dxu1 = (u1[i + 1][j] - u1[i - 1][j]) / 2.0 / dx;
                    //double dxu2 = (u2[i + 1][j] - u2[i - 1][j]) / 2.0 / dx;
                    //double dxu3 = (u3[i + 1][j] - u3[i - 1][j]) / 2.0 / dx;
                    //double dyu1 = (u1[i][j] - u1[i][j - 1]) / dy;
                    //double dyu2 = (u2[i][j] - u2[i][j - 1]) / dy;
                    //double dyu3 = (u3[i][j] - u3[i][j - 1]) / dy;

                    //double E = dxu1 + dyu2;
                    //double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
                    //double gam = E / E0;
                    //double exx = dxu1 - E / 3.0;
                    //double eyy = dyu2 - E / 3.0;
                    //double ezz = -E / 3.0;
                    //double exy = 0.5 * (dyu1 + dxu2);
                    //double exz = 0.5 * dxu3;
                    //double eyz = 0.5 * dyu3;
                    //double t1, t2, t3;
                    //if (step <= St)
                    //{
                    //    t1 = t2 = t3 = 0.0;
                    //}
                    //else
                    //{
                    //    if (mayak == 0)
                    //    {
                    //        t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                    //        t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                    //        t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                    //    }
                    //    else
                    //    {
                    //        t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                    //            FF2(gam) * AAk * Al * E * nx);
                    //        t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                    //            FF2(gam) * AAk * Al * E * ny);
                    //        t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                    //    }
                    //}

                    //uu1[i][j] = (-t1 + a * u1[i][j - 1] / 2.0 / dy - a * (u2[i + 1][j] - u2[i - 1][j]) / 4.0 / dx) / (a / 2.0 / dy);
                    //uu2[i][j] = (-t2 + (a + b) * u2[i][j - 1] / dy - b * (u1[i + 1][j] - u1[i - 1][j]) / 2.0 / dx) / ((a + b) / dy);

                    //uu3[i][j] = (F - t3 - (-s3 * u3[i - 1][j] - s4 * u3[i][j - 1])) / (+s1 + s2);
                    //}

                    // 1
                    if (x <= 0.0 && y >= 0.0)  // 1
                    {
                        double dxu1 = (u1[i][j] - u1[i - 1][j]) / dx;
                        double dxu2 = (u2[i][j] - u2[i - 1][j]) / dx;
                        double dxu3 = (u3[i][j] - u3[i - 1][j]) / dx;
                        double dyu1 = (u1[i][j + 1] - u1[i][j]) / dy;
                        double dyu2 = (u2[i][j + 1] - u2[i][j]) / dy;
                        double dyu3 = (u3[i][j + 1] - u3[i][j]) / dy;

                        double E = dxu1 + dyu2;
                        double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
                        double gam = E / E0;
                        double exx = dxu1 - E / 3.0;
                        double eyy = dyu2 - E / 3.0;
                        double ezz = -E / 3.0;
                        double exy = 0.5 * (dyu1 + dxu2);
                        double exz = 0.5 * dxu3;
                        double eyz = 0.5 * dyu3;
                        double t1, t2, t3;
                        if (step <= St)
                        {
                            t1 = t2 = t3 = 0.0;
                        }
                        else
                        {
                            if (mayak == 0)
                            {
                                t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                                t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                                t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                            }
                            else
                            {
                                t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                                    FF2(gam) * AAk * Al * E * nx);
                                t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                                    FF2(gam) * AAk * Al * E * ny);
                                t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                            }
                        }

                        double a1 = -k1 + k2 + k3;
                        double b1 = +k4 - k5;
                        double ff1 = -((-k6 - k7) * u1[i - 1][j] - k8 * u2[i - 1][j] + k9 * u1[i][j + 1] + k10 * u2[i][j + 1]) - t1;

                        double a2 = -l1 + l2;
                        double b2 = -l3 + l4 - l5;
                        double ff2 = -(-l6 * u1[i - 1][j] - l7 * u2[i - 1][j] + l8 * u1[i][j + 1] + (+l9 + l10) * u2[i][j + 1]) - t2;

                        /*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
                        uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
                        uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
                        uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;

                        uu3[i][j] = (-t3 - (-s3 * u3[i - 1][j] + s4 * u3[i][j + 1])) / (-s1 + s2);

                    }

                    // 2
                    else if (x > 0.0 && y >= 0.0)  // 2
                    {
                        double dxu1 = (u1[i + 1][j] - u1[i][j]) / dx;
                        double dxu2 = (u2[i + 1][j] - u2[i][j]) / dx;
                        double dxu3 = (u3[i + 1][j] - u3[i][j]) / dx;
                        double dyu1 = (u1[i][j + 1] - u1[i][j]) / dy;
                        double dyu2 = (u2[i][j + 1] - u2[i][j]) / dy;
                        double dyu3 = (u3[i][j + 1] - u3[i][j]) / dy;

                        double E = dxu1 + dyu2;
                        double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
                        double gam = E / E0;
                        double exx = dxu1 - E / 3.0;
                        double eyy = dyu2 - E / 3.0;
                        double ezz = -E / 3.0;
                        double exy = 0.5 * (dyu1 + dxu2);
                        double exz = 0.5 * dxu3;
                        double eyz = 0.5 * dyu3;
                        double t1, t2, t3;
                        if (step <= St)
                        {
                            t1 = t2 = t3 = 0.0;
                        }
                        else
                        {
                            if (mayak == 0)
                            {
                                t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                                t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                                t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                            }
                            else
                            {
                                t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                                    FF2(gam) * AAk * Al * E * nx);
                                t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                                    FF2(gam) * AAk * Al * E * ny);
                                t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                            }
                        }

                        double a1 = -k1 - k2 - k3;
                        double b1 = -k4 - k5;
                        double ff1 = -((+k6 + k7) * u1[i + 1][j] + k8 * u2[i + 1][j] + k9 * u1[i][j + 1] + k10 * u2[i][j + 1]) - t1;

                        double a2 = -l1 - l2;
                        double b2 = -l3 - l4 - l5;
                        double ff2 = -(+l6 * u1[i + 1][j] + l7 * u2[i + 1][j] + l8 * u1[i][j + 1] + (+l9 + l10) * u2[i][j + 1]) - t2;

                        /*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
                        uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
                        uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
                        uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;

                        uu3[i][j] = (-t3 - (+s3 * u3[i + 1][j] + s4 * u3[i][j + 1])) / (-s1 - s2);
                    }

                    // 3
                    else if (x > 0.0 && y < 0.0)  //3
                    {
                        double dxu1 = (u1[i + 1][j] - u1[i][j]) / dx;
                        double dxu2 = (u2[i + 1][j] - u2[i][j]) / dx;
                        double dxu3 = (u3[i + 1][j] - u3[i][j]) / dx;
                        double dyu1 = (u1[i][j] - u1[i][j - 1]) / dy;
                        double dyu2 = (u2[i][j] - u2[i][j - 1]) / dy;
                        double dyu3 = (u3[i][j] - u3[i][j - 1]) / dy;

                        double E = dxu1 + dyu2;
                        double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
                        double gam = E / E0;
                        double exx = dxu1 - E / 3.0;
                        double eyy = dyu2 - E / 3.0;
                        double ezz = -E / 3.0;
                        double exy = 0.5 * (dyu1 + dxu2);
                        double exz = 0.5 * dxu3;
                        double eyz = 0.5 * dyu3;
                        double t1, t2, t3;
                        if (step <= St)
                        {
                            t1 = t2 = t3 = 0.0;
                        }
                        else
                        {
                            if (mayak == 0)
                            {
                                t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                                t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                                t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                            }
                            else
                            {
                                t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                                    FF2(gam) * AAk * Al * E * nx);
                                t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                                    FF2(gam) * AAk * Al * E * ny);
                                t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                            }
                        }

                        double a1 = +k1 - k2 - k3;
                        double b1 = -k4 + k5;
                        double ff1 = -((+k6 + k7) * u1[i + 1][j] + k8 * u2[i + 1][j] - k9 * u1[i][j - 1] - k10 * u2[i][j - 1]) - t1;

                        double a2 = +l1 - l2;
                        double b2 = +l3 - l4 + l5;
                        double ff2 = -(+l6 * u1[i + 1][j] + l7 * u2[i + 1][j] - l8 * u1[i][j - 1] + (-l9 - l10) * u2[i][j - 1]) - t2;

                        ///*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
                        //uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
                        uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
                        uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;
                        uu3[i][j] = (-t3 - (+s3 * u3[i + 1][j] - s4 * u3[i][j - 1])) / (+s1 - s2);
                    }

                    // 4
                    else // (x <= 0.0 && y < 0.0)  //4
                    {
                        double dxu1 = (u1[i][j] - u1[i - 1][j]) / dx;
                        double dxu2 = (u2[i][j] - u2[i - 1][j]) / dx;
                        double dxu3 = (u3[i][j] - u3[i - 1][j]) / dx;
                        double dyu1 = (u1[i][j] - u1[i][j - 1]) / dy;
                        double dyu2 = (u2[i][j] - u2[i][j - 1]) / dy;
                        double dyu3 = (u3[i][j] - u3[i][j - 1]) / dy;

                        double E = dxu1 + dyu2;
                        double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
                        double gam = E / E0;
                        double exx = dxu1 - E / 3.0;
                        double eyy = dyu2 - E / 3.0;
                        double ezz = -E / 3.0;
                        double exy = 0.5 * (dyu1 + dxu2);
                        double exz = 0.5 * dxu3;
                        double eyz = 0.5 * dyu3;
                        double t1, t2, t3;
                        if (step <= St)
                        {
                            t1 = t2 = t3 = 0.0;
                        }
                        else
                        {
                            if (mayak == 0)
                            {
                                t1 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exx * nx + exy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * nx - BA * E0 * nx);
                                t2 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exy * nx + eyy * ny) + CA / (BA - kv(CA)) / BA * (CA * E * ny - BA * E0 * ny);
                                t3 = c23 * CA / (BA - kv(CA)) * (CA - gam) * (exz * nx + eyz * ny);
                            }
                            else
                            {
                                t1 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (dxu1 * nx + 0.5 * (dyu1 + dxu2) * ny - E * nx / 3.0) - //
                                    FF2(gam) * AAk * Al * E * nx);
                                t2 = pow(E0, Al - 2.0) * (-c23 * FF1(gam) * AAk * Al * (0.5 * (dyu1 + dxu2) * nx + dyu2 * ny - E * ny / 3.0) - //
                                    FF2(gam) * AAk * Al * E * ny);
                                t3 = pow(E0, Al - 2.0) * (-c23) * FF1(gam) * AAk * Al * 0.5 * (dxu3 * nx + dyu3 * ny);
                            }
                        }

                        double a1 = +k1 + k2 + k3;
                        double b1 = +k4 + k5;
                        double ff1 = -((-k6 - k7) * u1[i - 1][j] - k8 * u2[i - 1][j] - k9 * u1[i][j - 1] - k10 * u2[i][j - 1]) - t1;

                        double a2 = +l1 + l2;
                        double b2 = +l3 + l4 + l5;
                        double ff2 = -(-l6 * u1[i - 1][j] - l7 * u2[i - 1][j] - l8 * u1[i][j - 1] + (-l9 - l10) * u2[i][j - 1]) - t2;

                        /*uu1[i][j] = (b2 * ff1 - b1 * ff2) / (a1 * b2 - b1 * a2);
                        uu2[i][j] = (a2 * ff1 - a1 * ff2) / (b1 * a2 - b2 * a1);*/
                        uu1[i][j] = (ff1 - b1 * u2[i][j]) / a1;
                        uu2[i][j] = (ff2 - a2 * u1[i][j]) / b2;

                        uu3[i][j] = (-t3 - (-s3 * u3[i - 1][j] - s4 * u3[i][j - 1])) / (+s1 + s2);
                    }

                }
            }
        }


       

        // Пишем уравнения
        if (true)
        {
            for (int j = 1; j < M - 1; j++)
            {
                for (int i = 1; i < N - 1; i++)
                {
                    if (IF[i][j] == false)
                    {
                        dxu1 = (u1[i + 1][j] - u1[i - 1][j]) / (2.0 * dx);
                        dxu2 = (u2[i + 1][j] - u2[i - 1][j]) / (2.0 * dx);
                        dxu3 = (u3[i + 1][j] - u3[i - 1][j]) / (2.0 * dx);
                        dyu1 = (u1[i][j + 1] - u1[i][j - 1]) / (2.0 * dy);
                        dyu2 = (u2[i][j + 1] - u2[i][j - 1]) / (2.0 * dy);
                        dyu3 = (u3[i][j + 1] - u3[i][j - 1]) / (2.0 * dy);
                        ddxxu1 = (u1[i + 1][j] - 2 * u1[i][j] + u1[i - 1][j]) / (kv(dx));
                        ddxxu2 = (u2[i + 1][j] - 2 * u2[i][j] + u2[i - 1][j]) / (kv(dx));
                        ddxxu3 = (u3[i + 1][j] - 2 * u3[i][j] + u3[i - 1][j]) / (kv(dx));
                        ddyyu1 = (u1[i][j + 1] - 2 * u1[i][j] + u1[i][j - 1]) / (kv(dy));
                        ddyyu2 = (u2[i][j + 1] - 2 * u2[i][j] + u2[i][j - 1]) / (kv(dy));
                        ddyyu3 = (u3[i][j + 1] - 2 * u3[i][j] + u3[i][j - 1]) / (kv(dy));
                        ddxyu1 = (u1[i + 1][j + 1] - u1[i + 1][j - 1] - u1[i - 1][j + 1] + u1[i - 1][j - 1]) / (4.0 * dx * dy);
                        ddxyu2 = (u2[i + 1][j + 1] - u2[i + 1][j - 1] - u2[i - 1][j + 1] + u2[i - 1][j - 1]) / (4.0 * dx * dy);
                        ddxyu3 = (u3[i + 1][j + 1] - u3[i + 1][j - 1] - u3[i - 1][j + 1] + u3[i - 1][j - 1]) / (4.0 * dx * dy);

                        //if (i == 30 && j == 30 && step > 300)
                        //{
                        //    cout << "step = " << step << endl;
                        //    cout << "F1 = " << f1(u1[i][j], u2[i][j], u3[i][j], dxu1, dxu2, dxu3, dyu1, dyu2, dyu3, ddxxu1, ddxxu2, ddxxu3, //
                        //        ddyyu1, ddyyu2, ddyyu3, ddxyu1, ddxyu2, ddxyu3) << endl;
                        //    cout << "F2 = " << f2(u1[i][j], u2[i][j], u3[i][j], dxu1, dxu2, dxu3, dyu1, dyu2, dyu3, ddxxu1, ddxxu2, ddxxu3, //
                        //        ddyyu1, ddyyu2, ddyyu3, ddxyu1, ddxyu2, ddxyu3) << endl;
                        //    cout << "F3 = " << f3(u1[i][j], u2[i][j], u3[i][j], dxu1, dxu2, dxu3, dyu1, dyu2, dyu3, ddxxu1, ddxxu2, ddxxu3, //
                        //        ddyyu1, ddyyu2, ddyyu3, ddxyu1, ddxyu2, ddxyu3) << endl;
                        //}

                        if (step > St2)
                        {
                            uu1[i][j] = (f1(u1[i][j], u2[i][j], u3[i][j], dxu1, dxu2, dxu3, dyu1, dyu2, dyu3, ddxxu1, ddxxu2, ddxxu3, //
                                ddyyu1, ddyyu2, ddyyu3, ddxyu1, ddxyu2, ddxyu3) + B1 * u1[i + 1][j] + C1 * u2[i + 1][j - 1] + D1 * u1[i][j - 1] + //
                                E1 * u2[i - 1][j - 1] + F1 * u1[i - 1][j] + G1 * u2[i - 1][j + 1] + H1 * u1[i][j + 1] + I1 * u2[i + 1][j + 1]) / A1;
                            delta = max(delta, fabs(uu1[i][j] - u1[i][j]));

                            uu2[i][j] = (f2(u1[i][j], u2[i][j], u3[i][j], dxu1, dxu2, dxu3, dyu1, dyu2, dyu3, ddxxu1, ddxxu2, ddxxu3, //
                                ddyyu1, ddyyu2, ddyyu3, ddxyu1, ddxyu2, ddxyu3) + B2 * u2[i + 1][j] + C2 * u1[i + 1][j - 1] + D2 * u2[i][j - 1] + //
                                E2 * u1[i - 1][j - 1] + F2 * u2[i - 1][j] + G2 * u1[i - 1][j + 1] + H2 * u2[i][j + 1] + I2 * u1[i + 1][j + 1]) / A2;
                            delta = max(delta, fabs(uu2[i][j] - u2[i][j]));

                            uu3[i][j] = (f3(u1[i][j], u2[i][j], u3[i][j], dxu1, dxu2, dxu3, dyu1, dyu2, dyu3, ddxxu1, ddxxu2, ddxxu3, //
                                ddyyu1, ddyyu2, ddyyu3, ddxyu1, ddxyu2, ddxyu3) + B3 * u3[i + 1][j] + C3 * u3[i + 1][j - 1] + D3 * u3[i][j - 1] + //
                                E3 * u3[i - 1][j - 1] + F3 * u3[i - 1][j] + G3 * u3[i - 1][j + 1] + H3 * u3[i][j + 1] + I3 * u3[i + 1][j + 1]) / A3;
                            delta = max(delta, fabs(uu3[i][j] - u3[i][j]));
                        }
                        else
                        {
                            uu1[i][j] = (B1 * u1[i + 1][j] + C1 * u2[i + 1][j - 1] + D1 * u1[i][j - 1] + //
                                E1 * u2[i - 1][j - 1] + F1 * u1[i - 1][j] + G1 * u2[i - 1][j + 1] + H1 * u1[i][j + 1] + I1 * u2[i + 1][j + 1]) / A1;

                            delta = max(delta, fabs(uu1[i][j] - u1[i][j]));

                            uu2[i][j] = (B2 * u2[i + 1][j] + C2 * u1[i + 1][j - 1] + D2 * u2[i][j - 1] + //
                                E2 * u1[i - 1][j - 1] + F2 * u2[i - 1][j] + G2 * u1[i - 1][j + 1] + H2 * u2[i][j + 1] + I2 * u1[i + 1][j + 1]) / A2;
                            delta = max(delta, fabs(uu2[i][j] - u2[i][j]));

                            uu3[i][j] = (B3 * u3[i + 1][j] + C3 * u3[i + 1][j - 1] + D3 * u3[i][j - 1] + //
                                E3 * u3[i - 1][j - 1] + F3 * u3[i - 1][j] + G3 * u3[i - 1][j + 1] + H3 * u3[i][j + 1] + I3 * u3[i + 1][j + 1]) / A3;
                            delta = max(delta, fabs(uu3[i][j] - u3[i][j]));
                        }
                    }
                }
            }

            
            
        }

        c = u1;
        u1 = uu1;
        uu1 = c;

        c = u2;
        u2 = uu2;
        uu2 = c;

        c = u3;
        u3 = uu3;
        uu3 = c;
        
        if ( (step % 5000000 == 0 && step < 30000) || (step >= 30000 && step%2500000000 == 0) )
        {
            for (int i = 1; i < N - 1; i++)
            {
                for (int j = 1; j < M - 1; j++)
                {
                    if (IF[i][j] == false && j == M / 2)
                    {
                        double x = x_min + i * dx;
                        double y = y_min + j * dy;
                        dxu1 = (u1[i + 1][j] - u1[i - 1][j]) / (2.0 * dx);
                        dxu2 = (u2[i + 1][j] - u2[i - 1][j]) / (2.0 * dx);
                        dxu3 = (u3[i + 1][j] - u3[i - 1][j]) / (2.0 * dx);
                        dyu1 = (u1[i][j + 1] - u1[i][j - 1]) / (2.0 * dy);
                        dyu2 = (u2[i][j + 1] - u2[i][j - 1]) / (2.0 * dy);
                        dyu3 = (u3[i][j + 1] - u3[i][j - 1]) / (2.0 * dy);
                        double E = dxu1 + dyu2;
                        double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
                        double gam = E / E0;
                        double Exx = dxu1;
                        double Eyy = dyu2;
                        double Exy = 0.5 * (dyu1 + dxu2);
                        double Exz = 0.5 * dxu3;
                        double Eyz = 0.5 * dyu3;
                        double exx = dxu1 - E / 3.0;
                        double eyy = dyu2 - E / 3.0;
                        double ezz = -E / 3.0;
                        double exy = 0.5 * (dyu1 + dxu2);
                        double exz = 0.5 * dxu3;
                        double eyz = 0.5 * dyu3;
                        double sigma_xx = (c23 * (BA - CA * gam) * exx + (1.0 * E - CA * E0)) / (BA - kv(CA));
                        double sigma_yy = (c23 * (BA - CA * gam) * eyy + (1.0 * E - CA * E0)) / (BA - kv(CA));
                        double sigma_zz = (c23 * (BA - CA * gam) * ezz + (1.0 * E - CA * E0)) / (BA - kv(CA));
                        double sigma_xy = (c23 * (BA - CA * gam) * exy) / (BA - kv(CA));
                        double sigma_xz = (c23 * (BA - CA * gam) * exz) / (BA - kv(CA));
                        double sigma_yz = (c23 * (BA - CA * gam) * eyz) / (BA - kv(CA));

                        if (j == M / 2 && x > 0)
                        {
                            foutf << x << " " << u1[i][j] << " " << u2[i][j] << " "//
                                << u3[i][j] << " " << E << " " << Eyz << " " << Exy << " " << Exz << endl;
                        }
                    }

                }

            }
        }


    }

    foutf.close();
        
    ofstream fout;
    fout.open("file.txt");
    ofstream fout1;
    fout1.open("x_.txt");
    ofstream fout2;
    fout2.open("y_.txt");
    ofstream fout3;
    fout3.open("xy_.txt");
    ofstream fout4;
    fout4.open("hole_.txt");
    ofstream fout5;
    fout5.open("file2.txt");

    ofstream fout6;
    fout6.open("file_1.txt");

    fout << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"u1\", \"u2\", \"u3\", ZONE T = \"HP\" " << endl;
    fout1 << "TITLE = \"HP\"  VARIABLES = \"X\", \"u1\", \"u2\", \"u3\", \"sigma_xx\", \"sigma_yy\", \"sigma_zz\",\"sigma_xy\", \"sigma_xz\", \"sigma_yz\",\"E\", \"Eyz\", \"Exx\", \"Eyy\", \"Exy\",\"Exz\", ZONE T = \"HP\" " << endl;
    fout2 << "TITLE = \"HP\"  VARIABLES = \"Y\", \"u1\", \"u2\", \"u3\", \"sigma_xx\", \"sigma_yy\", \"sigma_zz\",\"sigma_xy\", \"sigma_xz\", \"sigma_yz\",\"E\", \"Eyz\", \"Exx\", \"Eyy\", \"Exy\",\"Exz\", ZONE T = \"HP\" " << endl;
    fout3 << "TITLE = \"HP\"  VARIABLES = \"X\", \"r\", \"u1\", \"u2\", \"u3\", \"sigma_xx\", \"sigma_yy\", \"sigma_zz\",\"sigma_xy\", \"sigma_xz\", \"sigma_yz\",\"E\", \"Eyz\", \"Exx\", \"Eyy\", \"Exy\",\"Exz\", ZONE T = \"HP\" " << endl;
    fout4 << "TITLE = \"HP\"  VARIABLES = \"angle\", \"r\", \"u1\", \"u2\", \"u3\", \"sigma_xx\", \"sigma_yy\", \"sigma_zz\",\"sigma_xy\", \"sigma_xz\", \"sigma_yz\",\"E\", \"Eyz\", \"Exx\", \"Eyy\", \"Exy\",\"Exz\", ZONE T = \"HP\" " << endl;
    fout5 << "TITLE = \"HP\"  VARIABLES = \"X\", \"Y\", \"sigma_xx\", \"sigma_yy\", \"sigma_zz\",\"sigma_xy\",\"sigma_xz\",\"sigma_yz\",\"E\",\"Eyz\", ZONE T = \"HP\" " << endl;
    
    for (int i = 0; i < N ; i++)
    {
        for (int j = 0; j < M ; j++)
        {
                fout6 << u1[i][j] << " " << u2[i][j] << " " << u3[i][j] << endl;
        }

    }

    for (int i = 1; i < N - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {
            if (IF[i][j] == false)
            {
                double x = x_min + i * dx;
                double y = y_min + j * dy;
                dxu1 = (u1[i + 1][j] - u1[i - 1][j]) / (2.0 * dx);
                dxu2 = (u2[i + 1][j] - u2[i - 1][j]) / (2.0 * dx);
                dxu3 = (u3[i + 1][j] - u3[i - 1][j]) / (2.0 * dx);
                dyu1 = (u1[i][j + 1] - u1[i][j - 1]) / (2.0 * dy);
                dyu2 = (u2[i][j + 1] - u2[i][j - 1]) / (2.0 * dy);
                dyu3 = (u3[i][j + 1] - u3[i][j - 1]) / (2.0 * dy);
                double E = dxu1 + dyu2;
                double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
                double gam = E / E0;
                double Exx = dxu1;
                double Eyy = dyu2;
                double Exy = 0.5 * (dyu1 + dxu2);
                double Exz = 0.5 * dxu3;
                double Eyz = 0.5 * dyu3;
                double exx = dxu1 - E / 3.0;
                double eyy = dyu2 - E / 3.0;
                double ezz = -E / 3.0;
                double exy = 0.5 * (dyu1 + dxu2);
                double exz = 0.5 * dxu3;
                double eyz = 0.5 * dyu3;
                double sigma_xx = (c23 * (BA - CA * gam) * exx + (1.0 * E - CA * E0)) / (BA - kv(CA));
                double sigma_yy = (c23 * (BA - CA * gam) * eyy + (1.0 * E - CA * E0)) / (BA - kv(CA));
                double sigma_zz = (c23 * (BA - CA * gam) * ezz + (1.0 * E - CA * E0)) / (BA - kv(CA));
                double sigma_xy = (c23 * (BA - CA * gam) * exy) / (BA - kv(CA));
                double sigma_xz = (c23 * (BA - CA * gam) * exz) / (BA - kv(CA));
                double sigma_yz = (c23 * (BA - CA * gam) * eyz) / (BA - kv(CA));

                fout << x << " " << y << " " << u1[i][j] << " " << u2[i][j] << " " << u3[i][j]  << endl;
                fout5 << x << " " << y << " " << sigma_xx << " " << sigma_yy << " " << sigma_zz << " " << sigma_xy << " " << //
                    sigma_xz << " " << sigma_yz << " " << E << " " << Eyz  << endl;
                if (i == j)
                {
                    fout3 << x << " " << sign(x) * sqrt(x*x + y*y) << " " << u1[i][j] << " " << u2[i][j] << " "//
                        << u3[i][j] << " " << sigma_xx << " " << sigma_yy << " " << sigma_zz << " " << sigma_xy << " " << //
                        sigma_xz << " " << sigma_yz << " " << E << " " << Eyz << " " << Exx << " " << Eyy << " " << Exy << " " << Exz << endl;
                }
                if (i == N/2)
                {
                    fout2 << y  << " " << u1[i][j] << " " << u2[i][j] << " "//
                        << u3[i][j] << " " << sigma_xx << " " << sigma_yy << " " << sigma_zz << " " << sigma_xy << " " << //
                        sigma_xz << " " << sigma_yz << " " << E << " " << Eyz << " " << Exx << " " << Eyy << " " << Exy << " " << Exz << endl;
                }
                if (j == M/2)
                {
                    fout1 << x << " "  << u1[i][j] << " " << u2[i][j] << " "//
                        << u3[i][j] << " " << sigma_xx << " " << sigma_yy << " " << sigma_zz << " " << sigma_xy << " " << //
                        sigma_xz << " " << sigma_yz << " " << E << " " << Eyz << " " << Exx << " " << Eyy << " " << Exy << " " << Exz << endl;
                }
            }
            else
            {
                double x = x_min + i * dx;
                double y = y_min + j * dy;
                double r = sqrt(x*x + y*y);
                if (r > R * 1.4 || r < R * 0.3)
                {
                    continue;
                }
                if (IF[i + 1][j] == false || IF[i + 1][j + 1] == false || IF[i][j+1] == false || IF[i - 1][j + 1] == false || IF[i - 1][j] == false//
                    || IF[i - 1][j - 1] == false || IF[i][j - 1] == false || IF[i + 1][j - 1] == false)
                {
                    if (x >= 0 && y >= 0)
                    {
                        dxu1 = (u1[i + 1][j] - u1[i][j]) / (dx);
                        dxu2 = (u2[i + 1][j] - u2[i][j]) / (dx);
                        dxu3 = (u3[i + 1][j] - u3[i][j]) / (dx);
                        dyu1 = (u1[i][j + 1] - u1[i][j]) / (dy);
                        dyu2 = (u2[i][j + 1] - u2[i][j]) / (dy);
                        dyu3 = (u3[i][j + 1] - u3[i][j]) / (dy);
                    }
                    else if (x >= 0 && y < 0)
                    {
                        dxu1 = (u1[i + 1][j] - u1[i][j]) / (dx);
                        dxu2 = (u2[i + 1][j] - u2[i][j]) / (dx);
                        dxu3 = (u3[i + 1][j] - u3[i][j]) / (dx);
                        dyu1 = (u1[i][j] - u1[i][j - 1]) / (dy);
                        dyu2 = (u2[i][j] - u2[i][j - 1]) / (dy);
                        dyu3 = (u3[i][j] - u3[i][j - 1]) / (dy);
                    }
                    else if (x < 0 && y >= 0)
                    {
                        dxu1 = (u1[i][j] - u1[i - 1][j]) / (dx);
                        dxu2 = (u2[i][j] - u2[i - 1][j]) / (dx);
                        dxu3 = (u3[i][j] - u3[i - 1][j]) / (dx);
                        dyu1 = (u1[i][j + 1] - u1[i][j]) / (dy);
                        dyu2 = (u2[i][j + 1] - u2[i][j]) / (dy);
                        dyu3 = (u3[i][j + 1] - u3[i][j]) / (dy);
                    }
                    else if (x < 0 && y < 0)
                    {
                        dxu1 = (u1[i][j] - u1[i - 1][j]) / (dx);
                        dxu2 = (u2[i][j] - u2[i - 1][j]) / (dx);
                        dxu3 = (u3[i][j] - u3[i - 1][j]) / (dx);
                        dyu1 = (u1[i][j] - u1[i][j - 1]) / (dy);
                        dyu2 = (u2[i][j] - u2[i][j - 1]) / (dy);
                        dyu3 = (u3[i][j] - u3[i][j - 1]) / (dy);
                    }
                    double E = dxu1 + dyu2;
                    double E0 = sqrt(c49 * (kv(dxu1) + kv(dyu2) - dxu1 * dyu2) + c13 * (kv(dyu1) + kv(dxu2) + 2.0 * dyu1 * dxu2 + kv(dxu3) + kv(dyu3)));
                    double gam = E / E0;
                    double Exx = dxu1;
                    double Eyy = dyu2;
                    double Exy = 0.5 * (dyu1 + dxu2);
                    double Exz = 0.5 * dxu3;
                    double Eyz = 0.5 * dyu3;
                    double exx = dxu1 - E / 3.0;
                    double eyy = dyu2 - E / 3.0;
                    double ezz = -E / 3.0;
                    double exy = 0.5 * (dyu1 + dxu2);
                    double exz = 0.5 * dxu3;
                    double eyz = 0.5 * dyu3;
                    double sigma_xx = (c23 * (BA - CA * gam) * exx + (1.0 * E - CA * E0)) / (BA - kv(CA));
                    double sigma_yy = (c23 * (BA - CA * gam) * eyy + (1.0 * E - CA * E0)) / (BA - kv(CA));
                    double sigma_zz = (c23 * (BA - CA * gam) * ezz + (1.0 * E - CA * E0)) / (BA - kv(CA));
                    double sigma_xy = (c23 * (BA - CA * gam) * exy) / (BA - kv(CA));
                    double sigma_xz = (c23 * (BA - CA * gam) * exz) / (BA - kv(CA));
                    double sigma_yz = (c23 * (BA - CA * gam) * eyz) / (BA - kv(CA));
                    fout4 << polar_angle(x,y) << " "  << r << " " << u1[i][j] << " " << u2[i][j] << " "//
                        << u3[i][j] << " " << sigma_xx << " " << sigma_yy << " " << sigma_zz << " " << sigma_xy << " " << //
                        sigma_xz << " " << sigma_yz << " " << E << " " << Eyz << " " << Exx << " " << Eyy << " " << Exy << " " << Exz << endl;
                }
            }

        }

    }

    fout.close();
}
