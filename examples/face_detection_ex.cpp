// The contents of this file are in the public domain. See LICENSE_FOR_EXAMPLE_PROGRAMS.txt
/*

    This example program shows how to find frontal human faces in an image.  In
    particular, this program shows how you can take a list of images from the
    command line and display each on the screen with red boxes overlaid on each
    human face.

    The examples/faces folder contains some jpg images of people.  You can run
    this program on them and see the detections by executing the following command:
        ./face_detection_ex faces/*.jpg

    
    This face detector is made using the now classic Histogram of Oriented
    Gradients (HOG) feature combined with a linear classifier, an image pyramid,
    and sliding window detection scheme.  This type of object detector is fairly
    general and capable of detecting many types of semi-rigid objects in
    addition to human faces.  Therefore, if you are interested in making your
    own object detectors then read the fhog_object_detector_ex.cpp example
    program.  It shows how to use the machine learning tools which were used to
    create dlib's face detector. 


    Finally, note that the face detector is fastest when compiled with at least
    SSE2 instructions enabled.  So if you are using a PC with an Intel or AMD
    chip then you should enable at least SSE2 instructions.  If you are using
    cmake to compile this program you can enable them by using one of the
    following commands when you create the build project:
        cmake path_to_dlib_root/examples -DUSE_SSE2_INSTRUCTIONS=ON
        cmake path_to_dlib_root/examples -DUSE_SSE4_INSTRUCTIONS=ON
        cmake path_to_dlib_root/examples -DUSE_AVX_INSTRUCTIONS=ON
    This will set the appropriate compiler options for GCC, clang, Visual
    Studio, or the Intel compiler.  If you are using another compiler then you
    need to consult your compiler's manual to determine how to enable these
    instructions.  Note that AVX is the fastest but requires a CPU from at least
    2011.  SSE4 is the next fastest and is supported by most current machines.  
*/


#include <dlib/image_processing/frontal_face_detector.h>
#include <dlib/gui_widgets.h>
#include <dlib/image_io.h>
#include <iostream>

using namespace dlib;
using namespace std;


//copy paste the following lines: start here
#include <fstream>
#include "assert.h"
#include "Operators.h"
#include "operatorFile_parser.h"
#include "setSubType.h"
#include "operandFile_parser.h"
#include "globals.h"
//#include "foo.h"
using namespace std;
extern hw_ac **myOp;   
//end here


// ----------------------------------------------------------------------------------------

//int main(int argc, char** argv)
//{  
int main(int argc, char* argv[]){
//    cout<< "test" << argc <<endl;
    string resultFolderName; 
    string resultFileName; 
    string operatorFileName;
    if (argc < 4) {
        cout<< "provide the name of the file that you want the result to be written to"<<endl;
        cout<< "Example: resultFolderName.txt resultFile.txt operatorFile.txt"<<endl; 
        return 0; 
    }else{
        cout << "argv[2] " << argv[2] << endl;
        resultFolderName= argv[2]; 
        resultFileName = argv[3]; 
        operatorFileName = argv[4]; 
    }
    assign_global_variables(resultFolderName, operatorFileName);
    string resultFileNameCompleteAddress = resultFileName;
    ofstream resultFile;
    resultFile.open(resultFileNameCompleteAddress.c_str(), ios_base::app);
    resultFile<<"*****************start******"<<endl; 
    //end here

    try
    {
        if (argc < 5)
        {
            cout << "Give some image files as arguments to this program." << endl;
            return 0;
        }

        frontal_face_detector detector = get_frontal_face_detector();
        image_window win;

        // Loop over all the images provided on the command line.
        for (int i = 1; i < (argc - 3); ++i)
        {
           // printf("HERE1\n");
            cout << "processing image " << argv[i] << endl;
           // printf("HERE2\n");
            array2d<unsigned char> img;
            load_image(img, argv[i]);
            // Make the image bigger by a factor of two.  This is useful since
            // the face detector looks for faces that are about 80 by 80 pixels
            // or larger.  Therefore, if you want to find faces that are smaller
            // than that then you need to upsample the image as we do here by
            // calling pyramid_up().  So this will allow it to detect faces that
            // are at least 40 by 40 pixels in size.  We could call pyramid_up()
            // again to find even smaller faces, but note that every time we
            // upsample the image we make the detector run slower since it must
            // process a larger image.
            pyramid_up(img);

            // Now tell the face detector to give us a list of bounding boxes
            // around all the faces it can find in the image.
            std::vector<rectangle> dets = detector(img);

            cout << "Number of faces detected: " << dets.size() << endl;
            // Now we show the image on the screen and the face detections as
            // red overlay boxes.
            win.clear_overlay();
            win.set_image(img);
            win.add_overlay(dets, rgb_pixel(255,0,0));

            cout << "Hit enter to process the next image..." << endl;
            cin.get();
        }
    }
    catch (exception& e)
    {
        cout << "\nexception thrown!" << endl;
        cout << e.what() << endl;
    }
}

// ----------------------------------------------------------------------------------------

