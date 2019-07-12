
/************************************************************/
/***  HomeWork3                                           ***/                                                     
/***  ZEINAB KOOHPEIMA, 901315                            ***/
/***   2 neurons                                          ***/
/************************************************************/

#include <cstdlib>
#include <iostream>
#include <cmath>
#include <math.h>
#include<fstream>
#include <time.h>

#define ESP .04
#define I 2.5

using namespace std;

int main(int argc, char *argv[])
{
    double vth,v1,v2,t, m1,m2, n1,n2, h1,h2, s2;
   
   double am=(-2.5)/(1-exp(2.5));
   double bm=4.0;
   double an=0.01*(-10)/(1-exp(1));
   double bn=0.125;
   double ah=0.07;
   double bh=1./(1+exp(3.0));

   m1 = am /(am+bm);
   n1 = an/(an+bn);
   h1 = ah/(ah +bh); 

   m2 = am /(am+bm);
   n2 = an/(an+bn);
   h2 = ah/(ah +bh);

    vth=-50.0;
    s2=0;
    v1=-65.0;
    v2=-65.0;
    
    ofstream outputV1;
   outputV1.open ("Voltge1.txt");
   
    ofstream outputV2;
   outputV2.open ("Voltage2.txt");

    ofstream outputS;
   outputS.open ("S.txt");
    
    for(int i =0 ; i<=1000 ; ++i)
{
     outputV1 << i << '\t' << v1 << '\n'; 
     outputV2 << i << '\t' << v2 << '\n';
     outputS << i << '\t' << s2 << '\n';
    
       double alpha_n1 = (0.01 * v1 + 0.55) / double (1-exp((-v1-55) *.1 ));
       double beta_n1 = 0.125 * exp(.0125 * (-v1-65));
       
       double alpha_m1 =( 0.1 * v1 + 4) / double(1-exp((-v1-40)*.1));
       double beta_m1 = 4 * exp(0.05*(-v1-65));
       
       double alpha_h1 = 0.07 * exp((-v1-60)*0.05);
       double beta_h1 = 1.0 / double(1+exp((-v1-30)*.1));
       
       double alpha_n2 = (0.01 * v2 + 0.55) / double(1-exp((-v2-55)*.1));
       double beta_n2 = 0.125 * exp(.0125*(-v2-65));
       
       double alpha_m2 =( 0.1 * v2 + 4) / double(1-exp((-v2-40)*.1));
       double beta_m2 = 4 * exp(0.05*(-v2-65));
       
       double alpha_h2 = 0.07 * exp((-v2-60)*0.05);
       double beta_h2 = 1.0 / double(1+exp((-v2-30)*.1));
       
       n1 +=  ESP *((alpha_n1)*(1-n1) - (beta_n1) *n1);
       m1 +=  ESP *((alpha_m1)*(1-m1) - (beta_m1)*m1);
       h1 +=  ESP *((alpha_h1)*(1-h1)- (beta_h1)*h1);
       v1 +=  ESP *(I-(120*pow(m1,3)*h1*(v1-55.72)+ 36*pow(n1,4)*(v1+77)+0.3*(v1+54.4)));
       
       n2 +=  ESP *((alpha_n2)*(1-n2) - (beta_n2) *n2);
       m2 +=  ESP *((alpha_m2)*(1-m2) - (beta_m2)*m2);
       h2 +=  ESP *((alpha_h2)*(1-h2)- (beta_h2)*h2);
       s2 +=  ESP *(.5*(1+tanh(5*(v1-vth)))*(1-s2)-.2*s2);
       
       v2 +=  ESP *(I -(120*pow(m2,3)*h2*(v2-55.72)+ 36*pow(n2,4)*(v2+77)+0.3*(v2+54.4))-.05*s2*(v2-20));

}

    outputV1.close(); 
    outputV2.close();
    outputS.close();
    
    system("PAUSE");
    return EXIT_SUCCESS;
}
