//#include <dlib/rand.h>
#include <dlib/svm_threaded.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
using namespace std;
//using namespace dlib;

//typedef matrix<double,2,1> sample_type;


int main(void){
    ofstream datafile;
    datafile.open ("svm_data_test.txt");
    
   // dlib::rand rnd;
    double x;
    double y; 
    double r;
    double theta;
    
    for(int i=0; i < 50; i++){
       r = ((double)rand() / RAND_MAX)*10.0;
       while (r > 10){
         r = ((double)rand() / RAND_MAX)*10.0;
       }
       if (r < 0){r*=-1;}
       theta = (double)rand();
       x = r*cos(theta);
       if (x < 0){x*=-1;}
       y = r*sin(theta);
       if (y<0){y*=-1;}
       datafile << x << " " << y << " " << "+1\n";
    }  
    
    for(int i=0; i < 50; i++){
       r = ((double)rand() / RAND_MAX)*10.0 + 10.0;
       while (r < 10){
         r = ((double)rand() / RAND_MAX)*10.0 + 10.0;
       }
       theta = (double)rand();
       x = r*cos(theta);
       y = r*sin(theta);
       datafile << x << " " << y << " " << "-1\n";
    }  
    datafile.close();

} 
