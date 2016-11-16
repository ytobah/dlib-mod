#include <dlib/rand.h>
#include <dlib/svm_threaded.h>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;
using namespace dlib;

typedef matrix<double,2,1> sample_type;


int main(void){
    ofstream datafile;
    datafile.open ("multi_data_test.txt");
    
    const long num = 50;

    sample_type m;

    dlib::rand rnd;


    // make some samples near the origin
    double radius = 0.5;
    for (long i = 0; i < num+10; ++i)
    {
        double sign = 1;
        if (rnd.get_random_double() < 0.5)
            sign = -1;
        m(0) = 2*radius*rnd.get_random_double()-radius;
        m(1) = sign*sqrt(radius*radius - m(0)*m(0));

        // add this sample to our set of training samples 
//        samples.push_back(m);
 //       labels.push_back(1);
        datafile << m(0) << " " << m(1) << " 1\n";
    }

    // make some samples in a circle around the origin but far away
    radius = 10.0;
    for (long i = 0; i < num+20; ++i)
    {
        double sign = 1;
        if (rnd.get_random_double() < 0.5)
            sign = -1;
        m(0) = 2*radius*rnd.get_random_double()-radius;
        m(1) = sign*sqrt(radius*radius - m(0)*m(0));

        // add this sample to our set of training samples 
    //    samples.push_back(m);
     //   labels.push_back(2);
        datafile << m(0) << " " << m(1) << " 2\n";
    }

    // make some samples in a circle around the point (25,25) 
    radius = 4.0;
    for (long i = 0; i < num+30; ++i)
    {
        double sign = 1;
        if (rnd.get_random_double() < 0.5)
            sign = -1;
        m(0) = 2*radius*rnd.get_random_double()-radius;
        m(1) = sign*sqrt(radius*radius - m(0)*m(0));

        // translate this point away from the origin
        m(0) += 25;
        m(1) += 25;


        datafile << m(0) << " " << m(1) << " 3\n";
    }
    datafile.close();
}
