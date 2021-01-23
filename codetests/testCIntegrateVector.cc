#include "integratevec.h"
#include <iostream>
using namespace std;

void integrand(int* ndim, double* q, int* numfunc, double* f){
	f[0]=q[0];
	f[1]=q[1]*q[1];
}


int main(void){
	CIntegrateVector junk;
    junk.SetNDim(2);
    junk.SetNumFunc(2);
    junk.SetLimits(0, 0.0, 1.0);
    junk.SetLimits(1, 2.0, 2.2);


    cout << "Check Member access"<< endl;
    cout << "GetNDim :" << 2 << " == " << junk.GetNDim() << " ?" << endl;
    cout << "GetNumFunc :" << 2 << " == " << junk.GetNumFunc() << " ?"<< endl;
    cout << "GetMinPts :" << 10 << " == " << junk.GetMinPts() << " ?"<< endl;
    cout << "GetMaxPts :" << 100 << " == " << junk.GetMaxPts() << " ?"<< endl;
    cout << "GetKey :" << 0 << " == " << junk.GetKey() << " ?"<< endl;
    cout << "GetNW :" << 0 << " == " << junk.GetNW() << " ?"<< endl;
    cout << "GetRestart :" << 0 << " == " << junk.GetRestart() << " ?"<< endl;
    cout << "GetAbsErr :" << 1e-14 << " == " << junk.GetAbsErr() << " ?"<< endl;
    cout << "GetRelErr :" << 1e-6 << " == " << junk.GetRelErr() << " ?"<< endl;
    cout << "GetUpperLimit(0) :" << 1 << " == " << junk.GetUpperLimit(0) << " ?"<< endl;
    cout << "GetLowerLimit(0) :" << 0 << " == " << junk.GetLowerLimit(0) << " ?"<< endl;
    cout << "GetUpperLimit(1) :" << 2.2 << " == " << junk.GetUpperLimit(1) << " ?"<< endl;
    cout << "GetLowerLimit(1) :" << 2 << " == " << junk.GetLowerLimit(1) << " ?"<< endl;
    cout << "GetNumEvals :" << 0 << " == " << junk.GetNumEvals() << " ?"<< endl;
    cout << "GetIFail :" << 0 << " == " << junk.GetIFail() << " ?"<< endl;

    cout << endl << "Integrate" << endl;
    junk.Compute(integrand);
    cout << 1.0 << " == " << junk.GetResults(0) << " +/- " << junk.GetError(0) << " ?"<<endl;
    cout << 0.88266667 << " == " << junk.GetResults(1) << " +/- " << junk.GetError(1) << " ?" << endl;


}
 
