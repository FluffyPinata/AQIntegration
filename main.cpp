#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <omp.h>
#include <list>
#include <utility>
#include <exception>
#include <cmath>
#include <time.h>
#include <stdlib.h>

using namespace std;

struct Result
{
    double area;
};

double f(const double x);
const Result rectangleMethod(const double, const double, const double, const int);
const Result trapezoidalMethod(const double, const double, const double, const int);

int main()
{
    const short maxThreads = 10;
    short method;
    double x1, x2, dx;

    cout << fixed << setprecision(8) << endl;
    try
    {
        while (true)
        {
            cout << "   X1: "; cin >> x1;
            cout << "   X2: "; cin >> x2;
            cout << "   dx: "; cin >> dx;
            cout << "   Method (1 - rectangle, 2 - trapezoidal): "; cin >> method;

            list<pair<short, Result> > results;
            for (int i = 0; i < maxThreads; i++)
            {
                Result result = (method == 1) ?
                                rectangleMethod(x1, x2, dx, i + 1) :
                                trapezoidalMethod(x1, x2, dx, i + 1);

                pair<short, Result> s_result(i + 1, result);
                results.push_back(s_result);
            }

            cout << endl << "   Results:" << endl;
            for (int i = 0; i < results.size(); i++)
            {
                cout << "   Threads: " << results.front().first;
                cout << ", area: " << results.front().second.area << endl;
                results.pop_front();
            }
            cout << endl;
        }
    }
    catch (exception & e)
    {
        cout << e.what() << endl;
    }
    cin.get();
    return 0;
}

const Result rectangleMethod(const double x1, const double x2, const double dx, const int nThreads)
{
    const int N = static_cast<int>((x2 - x1) / dx);
    double s = 0;

    #pragma omp parallel for num_threads(nThreads) reduction(+: s)
    for (int i = 1; i <= N; i++) s += f(x1 + i * dx);

    s *= dx;

    Result temp;
    temp.area = s;

    return temp;
}

const Result trapezoidalMethod(const double x1, const double x2, const double dx, const int nThreads)
{
    const int N = static_cast<int>((x2 - x1) / dx);
    double s = 0;

    #pragma omp parallel for num_threads(nThreads) reduction(+: s)
    for (int i = 1; i < N; i++) s += f(x1 + i * dx);

    s = (s + (f(x1) + f(x2)) / 2) * dx;

    Result temp;
    temp.area = s;

    return temp;
}

double f(const double x)
{
    return sin(x);
}