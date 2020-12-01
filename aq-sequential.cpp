#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <omp.h>
#include <stack>
#include <time.h>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <stdio.h>
#include <unistd.h>
#include <bits/stdc++.h>

using namespace std;

double f1(double x){
    double result = (4 * pow(x, 6) - 2 * pow(x, 3) + (7 * x) - 4);
    return result;
}

double f2(double x){
    double result = cos(x) - (3 / pow(x, 5));
    return result;
}

double f3(double x){
    double result = exp(x) + 1 / (pow(x, 2) + 1);
    return result;
}

double integralF1(double a, double b){
    double integral;
    integral = (((4 * pow(b, 7)) / 7) - ((2 * pow(b, 4)) / 4) + ((7 * pow(b, 2)) / 2) - 4 * b)
                    - (((4 * pow(a, 7)) / 7) - ((2 * pow(a, 4)) / 4) + ((7 * pow(a, 2)) / 2) - 4 * a);
    return integral;
}

double integralF2(double a, double b){
    double integral;
    integral = (sin(b) + (3 / (4 * pow(b, 4)))) - (sin(a) + (3 / (4 * pow(a, 4))));
    return integral;
}

double integralF3(double a, double b){
    double integral;
    integral = (exp(b) + atan(b)) - (exp(a) + atan(a));
    return integral;
}

double SimpsonsRule(double a, double b, int func){
    double integral = 0.0;
    double c = (a + b) / 2;
    if (func == 1){
        integral = ((b - a) * (f1(a) + 4 * f1(c) + f1(b))) / 6;
    }
    else if (func == 2){
        integral = ((b - a) * (f2(a) + 4 * f2(c) + f2(b))) / 6;
    }
    else if (func == 3){
        integral = ((b - a) * (f3(a) + 4 * f3(c) + f3(b))) / 6;
    }

    return integral;
}

double getError(double val1, double val2){
    return (abs((val1 - val2) / val1));
}

double AdaptiveQuadrature(double lower, double upper, double error, int func){
    double simpsons = SimpsonsRule(lower, upper, func);
    double actual;
    if (func == 1){
        actual = integralF1(lower, upper);
    }
    else if (func == 2){
        actual = integralF2(lower, upper);
    }
    else if (func == 3){
        actual = integralF3(lower, upper);
    }
    double integral = 0.0;

    if (getError(actual, simpsons) > error){
        double midpoint = (upper + lower) / 2;
        integral += AdaptiveQuadrature(lower, midpoint, error, func);
        integral += AdaptiveQuadrature(midpoint, upper, error, func);
    }
    else{
        integral += simpsons;
    }

    return integral;
}

int main()
{
    double approxf1, approxf2, approxf3, actualf1, actualf2, actualf3, error1, error2, error3, runtime;
    
    //runtime = clock()/(double)CLOCKS_PER_SEC;
    runtime = omp_get_wtime();

    approxf1 = AdaptiveQuadrature(0, 10, 0.02, 1);
    approxf2 = AdaptiveQuadrature(1, 8, 0.02, 2);
    approxf3 = AdaptiveQuadrature(0, 10, 0.02, 3);

    actualf1 = integralF1(0, 10);
    actualf2 = integralF2(1, 8);
    actualf3 = integralF3(0, 10);

    error1 = getError(actualf1, approxf1);
    error2 = getError(actualf2, approxf2);
    error3 = getError(actualf3, approxf3);

    runtime = omp_get_wtime() - runtime;
    //runtime = (clock()/(double)CLOCKS_PER_SEC ) - runtime;

    cout << "Approximate integral of f(x) = 4x^6 - 2x^3 + 7x - 4 from 0 to 10: " << approxf1 << endl;
    cout << "Actual integral of f(x) = 4x^6 - 2x^3 + 7x - 4 from 0 to 10: " << actualf1 << endl;
    cout << "Error: " << error1 << "\n\n";
    cout << "Approximated f(x) = cos(x) - 3x^-5 from 1 to 8: " << approxf2 << endl;
    cout << "Actual f(x) = cos(x) - 3x^-5 from 1 to 8: " << actualf2 << endl;
    cout << "Error: " << error2 << "\n\n";
    cout << "Approximated f(x) = e^x + 1 / (x^2 + 1) from 0 to 10: " << approxf3 << endl;
    cout << "Actual f(x) = e^x + 1 / (x^2 + 1) from 0 to 10: " << actualf3 << endl;
    cout << "Error: " << error3 << "\n\n";
    cout << "Program runs in " << setiosflags(ios::fixed) << setprecision(8) << runtime << " seconds\n";

    return 0;
}