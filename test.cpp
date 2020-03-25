#include<iostream>
#include<fstream>
extern "C"{
#include<math.h>
#include<time.h>
#include<unistd.h>
#include<sys/time.h>
}
#include "IMParton.h"
using namespace std;


IMParton proton(1,1);              // Z=1 and A=1, which is proton.
IMParton neutron(0,1);             // Z=0 and A=1, which is neutron.

int main()
{
    //using data set B
    proton.setDataSet(2);             //1 for data set A and 2 data set B
    neutron.setDataSet(2);             //1 for data set A and 2 data set B

    struct timeval tstart,tend;
    gettimeofday(&tstart,NULL);       //get the starting time

    double lnxstep=(log(1)-log(0.1))/100;       //step length of log(x) for output data
    ofstream outdata1("proton_Q2_0.5_xuv.txt");         //output file for xuv of proton at Q^2 = 0.5 GeV^2
    ofstream outdata2("neutron_Q2_0.5_xdv.txt");        //output file for xdv of neutron at Q^2 = 0.5 GeV^2
    ofstream outdata3("proton_Q2_2.5_xsbar.txt");       //output file for xsbar at Q^2 = 2.5 GeV^2
    ofstream outdata4("proton_Q2_10_xgluon.txt");       //output file for xgluon at Q^2 = 10 GeV^2
    ofstream outdata5("proton_Q2_54_dbar_over_ubar.txt");     //output file for dbar/ubar at Q^2 = 54 GeV^2
    //get the parton distribution functions
    for(int i=0;i<400;i++)
    {
        double x = exp(log(1e-4)+i*lnxstep);
        outdata1<<x<<" "<<x*(proton.getPDF(1,x,0.5)-proton.getPDF(-1,x,0.5))<<endl;
        outdata2<<x<<" "<<x*(neutron.getPDF(2,x,0.5)-neutron.getPDF(-2,x,0.5))<<endl;
        outdata3<<x<<" "<<x*proton.getPDF(-3,x,2.5)<<endl;
        outdata4<<x<<" "<<x*proton.getPDF(0,x,10)<<endl;
        outdata5<<x<<" "<<(proton.getPDF(-2,x,54)/proton.getPDF(-1,x,54))<<endl;
    }

    gettimeofday(&tend,NULL);         //get the ending time
    double timeuse=1000000*(tend.tv_sec-tstart.tv_sec)+(tend.tv_usec-tstart.tv_usec);
    //display the total program running time
    cout<<"    Runtime : "<<(timeuse/1000.0)<<" ms."<<endl;
}
