#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <omp.h>
#include <stack>
#include <time.h>
#include <cstdlib>
#include <iomanip>
#include "mpi.h"

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

void AdaptiveQuadrature(double lower, double upper, double error, int func, int rank, int size, MPI_Status status){
    if (rank != 0) {
        double simpsons = SimpsonsRule(lower, upper, func);
        double actual;

        while (rank < rank^2) { //ensure processes are on correct layers
            rank += 1;
        }

        if (func == 1){
            actual = integralF1(lower, upper);
        }
        else if (func == 2){
            actual = integralF2(lower, upper);
        }
        else if (func == 3){
            actual = integralF3(lower, upper);
        }
        double integral_result, left_area, right_area;

        if (getError(actual, simpsons) > error){
            double midpoint = (upper + lower) / 2;
            
            AdaptiveQuadrature(lower, midpoint, error, func, rank, size, status);
            MPI_Send(&left_area, 1, MPI_INT, rank, FROM_WORKER, MPI_COMM_WORLD, &status);
               
            AdaptiveQuadrature(midpoint, upper, error, func, rank, size, status);
            MPI_Send(&right_area, 1, MPI_INT, rank, FROM_WORKER, MPI_COMM_WORLD, &status);
            
            MPI_Recv(&left_area, 1, MPI_INT, rank, FROM_WORKER, MPI_COMM_WORLD, &status);
            MPI_Recv(&right_area, 1, MPI_INT, rank, FROM_WORKER, MPI_COMM_WORLD, &status);
            integral_result = left_area + right_area;
            //Send back to master
            MPI_Send(MPI_Send(&integral_result, 1, MPI_INT, 0, FROM_WORKER, MPI_COMM_WORLD);
        }
        else {
            integral_result = simpsons;
            //Send back to master
            MPI_Send(&integral_result, 1, MPI_INT, 0, FROM_WORKER, MPI_COMM_WORLD);
        }
    }
}

int main()
{
    double approxf1, approxf2, approxf3, actualf1, actualf2, actualf3, error1, error2, error3, runtime;
    int size, rank;

    MPI_Status status;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    if (size < 2 ) {
        printf("Need at least two MPI tasks. Quitting...\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
        exit(1);
    }

    if (rank == 0) {
        runtime = MPI_Wtime();
    }

    source = 0;

    if (rank == 0) {
        approxf1 = AdaptiveQuadrature(0, 10, 0.02, 1, rank, size, status);
        approxf2 = AdaptiveQuadrature(1, 8, 0.02, 2, rank, size, status);
        approxf3 = AdaptiveQuadrature(0, 10, 0.02, 3, rank, size, status);   
    }

    actualf1 = integralF1(0, 10);
    actualf2 = integralF2(1, 8);
    actualf3 = integralF3(0, 10);

    error1 = getError(actualf1, approxf1);
    error2 = getError(actualf2, approxf2);
    error3 = getError(actualf3, approxf3);

    runtime = MPI_Wtime() - runtime;

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