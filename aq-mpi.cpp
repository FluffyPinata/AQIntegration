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
#include <stack>
#include <vector>
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

float master(double lower, double upper, double error, int func, int numtasks, int taskid) {
    double bounds[5] = {lower, upper, error, func, 0}; //Last is where result is stored
    // idk how to send vector into mpi so just convert between array/vector as needed
    vector<double> tempVec(bounds, bounds + 5);
    double integral_result = 0;
    //can't stick an array on a stack so this is the next best option
    stack< vector<double> > tasks;
    //Stores current availability of each task
    bool tasksAvailable[numtasks - 1];
    //Checks if program is still running even if stack is empty
    bool tasksRunning = false;
    tasks.push(tempVec);

    MPI_Status status;

    //either stuffs still in the stack to be assigned or tasks aren't done with what they're working on
    while (!tasks.empty() || tasksRunning) {
        MPI_Recv(bounds, 5, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        tasksAvailable[status.MPI_SOURCE - 1] = true; //Set received task as being done

        if (status.MPI_TAG == 1) { //flag for pushing back to stack
            //Received flag to split into two processes
            vector<double> tempLeftVec(bounds, bounds + 5);
            tasks.push(tempLeftVec);
            MPI_Recv(bounds, 5, MPI_DOUBLE, status.MPI_SOURCE, status.MPI_TAG, MPI_COMM_WORLD, &status);
            vector<double> tempRightVec(bounds, bounds + 5);
            tasks.push(tempRightVec);
        } else if (status.MPI_TAG == 2) { //flag for adding to result
            //Received result, now add it to total
            integral_result += bounds[4];
        }

        for (int i = 0; i < numtasks - 1; i++) {
            //if still work to assign and task is available, loop through tasks to find 
            if (!tasks.empty() && tasksAvailable[i] == true) {
                tempVec = tasks.top();
                tasks.pop();
                bounds[0] = tempVec[0];
                bounds[1] = tempVec[1];
                bounds[2] = tempVec[2];
                bounds[3] = tempVec[3];
                bounds[4] = 0; //This should always be 0 if it gets here but just in case
                MPI_Send(bounds, 5, MPI_DOUBLE, i + 1, 0, MPI_COMM_WORLD);
                tasksAvailable[i] = false;
            }
        }

        //Catching if there are still left over tasks after stack is cleared. basically to stop any premature cancellation
        tasksRunning = false;
        for (int i = 0; i < numtasks - 1; i++) {
            if (tasksAvailable[i] == false) {
                tasksRunning = true;
            }
        }
    }
    //Send message to cancel all slave tasks
    for (int i = 0; i < numtasks - 1; i++) {
        MPI_Send(bounds, 5, MPI_DOUBLE, i + 1, 5, MPI_COMM_WORLD);
    }   

    return integral_result;

}

void slave(int taskid) {
    double bounds[5] = {0, 0, 0, 0, 0};
    double lower, upper, error, func, simpsons, actual;

    MPI_Status status;

    MPI_Send(bounds, 5, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

    while (true) { //will run until it recieves flag to cancel (i.e got a result)
        MPI_Recv(bounds, 5, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        if (status.MPI_TAG == 5) { //Kill task when it recieves cancelation message
            break;
        } else {
            lower = bounds[0];
            upper = bounds[1];
            error = bounds[2];
            func = bounds[3];

            // this is just the code from AdaptiveQuadature function
            simpsons = SimpsonsRule(lower, upper, func);
            if (func == 1){
                actual = integralF1(lower, upper);
            }
            else if (func == 2){
                actual = integralF2(lower, upper);
            }
            else if (func == 3){
                actual = integralF3(lower, upper);
            }

            //Spawn two new tasks and send to be requeued in stack
            if (getError(actual, simpsons) > error) { // 1 is flag for needs to be pushed to stack
                double midpoint = (upper + lower) / 2;
                bounds[0] = lower;
                bounds[1] = midpoint;
                MPI_Send(bounds, 5, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                bounds[0] = midpoint;
                bounds[1] = upper;
                MPI_Send(bounds, 5, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
            }
            //Send result back to master
            else { // 2 is flag for got a result
                bounds[4] = simpsons;
                MPI_Send(bounds, 5, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
            }
        }
    }
}

int main (int argc, char *argv[])
{
    double approxf1, approxf2, approxf3, actualf1, actualf2, actualf3, error1, error2, error3, runtime;
    int taskid, numtasks;

    MPI_Status status;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    if (numtasks < 2 ) {
        printf("Need at least two MPI tasks. Quitting...\n");
        MPI_Abort(MPI_COMM_WORLD, 0);
        exit(1);
    }

    runtime = MPI_Wtime();

    if (taskid == 0) {
        //0 is master task
        approxf1 = master(0, 10, 0.02, 1, numtasks, taskid);
        //approxf2 = master(1, 8, 0.02, 2, numtasks, taskid);
        //approxf3 = master(0, 10, 0.02, 3, numtasks, taskid);
    } else {
        //All others are slaves
        slave(taskid);
    }

    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    if (taskid == 0) {
        //0 is master task
        //approxf1 = master(0, 10, 0.02, 1, numtasks, taskid);
        approxf2 = master(1, 8, 0.02, 2, numtasks, taskid);
        //approxf3 = master(0, 10, 0.02, 3, numtasks, taskid);
    } else {
        //All others are slaves
        slave(taskid);
    }

    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    MPI_Comm_size(MPI_COMM_WORLD,&numtasks);

    if (taskid == 0) {
        //0 is master task
        //approxf1 = master(0, 10, 0.02, 1, numtasks, taskid);
        //approxf2 = master(1, 8, 0.02, 2, numtasks, taskid);
        approxf3 = master(0, 10, 0.02, 3, numtasks, taskid);
    } else {
        //All others are slaves
        slave(taskid);
    }

    if (taskid == 0) {
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
    }

    MPI_Finalize();
    return 0;
}