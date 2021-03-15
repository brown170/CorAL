#include "integrator.h"
#include <iostream>
using namespace std;

/*static void integrand(int* ndim, double* q, int* numfunc, double* f){
	cout << q << endl;
	f[0]=q[0];
	f[1]=q[1]*q[1];
}*/

//The integral function to be tested
int f(unsigned ndim, const double *x, void *fdata, unsigned fdim, double *fval) {
    
    double *sigma = ((double *) fdata);//turning fdata into a pointer for an array of double
    fval[0] = (*(sigma+0))*x[0];
    fval[1] =(*(sigma+1))*x[1]*x[1];
    return 0; // success*
}




    int main(void){


	CIntegrateCubature junk;
	//Initialization
    double sigma1[2] = {1,1};//passing parameters for the function
    junk.set_ndim(2);
    junk.set_fdim(2);
    //junk.SetMaxPts(200);
    junk.set_limit(0, 0, 1);
    junk.set_limit(1, 2, 2.2);

    //Member Access
    cout << "Check Member access"<< endl;
    cout << "GetNDim :" << 2 << " == " << junk.get_ndim() << " ?" << endl;
    cout << "GetFDim :" << 2 << " == " << junk.get_fdim() << " ?"<< endl;
    //cout << "GetMinPts :" << 10 << " == " << junk.GetMinPts() << " ?"<< endl;
    //cout << "GetMaxPts :" << 200 << " == " << junk.GetMaxPts() << " ?"<< endl;
    //cout << "GetKey :" << 0 << " == " << junk.GetKey() << " ?"<< endl;
    //cout << "GetNW :" << 0 << " == " << junk.GetNW() << " ?"<< endl;
    //cout << "GetRestart :" << 0 << " == " << junk.GetRestart() << " ?"<< endl;
    cout << "GetAbsErr :" << 1e-14 << " == " << junk.get_abserr() << " ?"<< endl;
    cout << "GetRelErr :" << 1e-6 << " == " << junk.get_relerr() << " ?"<< endl;
    cout << "GetUpperLimit(0) :" << 1 << " == " << junk.get_upper_limit(0) << " ?"<< endl;
    cout << "GetLowerLimit(0) :" << 0 << " == " << junk.get_lower_limit(0) << " ?"<< endl;
    cout << "GetUpperLimit(1) :" << 2.2 << " == " << junk.get_upper_limit(1) << " ?"<< endl;
    cout << "GetLowerLimit(1) :" << 2 << " == " << junk.get_lower_limit(1) << " ?"<< endl;
    cout << "GetNumEvals :" << 0 << " == " << junk.get_neval() << " ?"<< endl;
    //cout << "GetIFail :" << 0 << " == " << junk.GetIFail() << " ?"<< endl;

    cout << endl << "Integrate" << endl;
    
    junk.compute(f,&sigma1);
    cout << 0.1 << " == " << junk.value[0]<< " +/- " << junk.error[0]<< " ?"<<endl;
    cout << 0.88266667 << " == " << junk.value[1]<< " +/- " << junk.error[1] << " ?" << endl;
    
    return 1;
}
 
