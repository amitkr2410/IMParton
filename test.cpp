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


IMParton proton(1,1);             // Z=1 and A=1, which is proton. And so far proton PDF only

int main()
{
    //using data set B
    proton.setDataSet(2);             //1 for data set A and 2 data set B

    struct timeval tstart,tend;
    gettimeofday(&tstart,NULL);       //get the starting time

    double lnxstep=(log(1)-log(0.1))/100;       //step length of log(x) for output data
    ofstream outdata1("Q2_0.1_xuv.dat");        //output file for xuv at Q^2 = 0.1 GeV^2
    ofstream outdata2("Q2_0.7_xdv.dat");        //output file for xdv at Q^2 = 0.7 GeV^2
    ofstream outdata3("Q2_2.5_xsbar.dat");      //output file for xsbar at Q^2 = 2.5 GeV^2
    ofstream outdata4("Q2_10_xgluon.dat");      //output file for xgluon at Q^2 = 10 GeV^2
    ofstream outdata5("Q2_54_dbar_over_ubar.dat");     //output file for dbar/ubar at Q^2 = 54 GeV^2
    //get the parton distribution functions
    for(int i=0;i<400;i++)
    {
        double x = exp(log(1e-4)+i*lnxstep);
        outdata1<<x<<" "<<x*(proton.getPDF(1,x,0.1)-proton.getPDF(-1,x,0.1))<<endl;
        outdata2<<x<<" "<<x*(proton.getPDF(2,x,0.7)-proton.getPDF(-2,x,0.7))<<endl;
        outdata3<<x<<" "<<x*proton.getPDF(-3,x,2.5)<<endl;
        outdata4<<x<<" "<<x*proton.getPDF(0,x,10)<<endl;
        outdata5<<x<<" "<<(proton.getPDF(-2,x,54)/proton.getPDF(-1,x,54))<<endl;
    }

    gettimeofday(&tend,NULL);         //get the ending time
    double timeuse=1000000*(tend.tv_sec-tstart.tv_sec)+(tend.tv_usec-tstart.tv_usec);
    //display the total program running time
    cout<<"    Runtime : "<<(timeuse/1000.0)<<" ms."<<endl;
}
