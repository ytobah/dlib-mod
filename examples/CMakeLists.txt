#
# This is a CMake makefile.  You can find the cmake utility and
# information about it at http://www.cmake.org
#


cmake_minimum_required(VERSION 2.8.4)

PROJECT(examples)

include(../dlib/cmake)
include_directories(/home/ytobah/simple_apx_app/inc/BookKeeping /home/ytobah/simple_apx_app/inc/OperatorHeaders)  #adding the include folder



# Tell CMake to compile a program.  We do this with the ADD_EXECUTABLE()
# statement which takes the name of the output executable and then a list of
# .cpp files to compile.  Here each example consists of only one .cpp file but
# in general you will make programs that const of many .cpp files.
ADD_EXECUTABLE(assignment_learning_ex assignment_learning_ex.cpp)
# Then we tell it to link with dlib.
TARGET_LINK_LIBRARIES(assignment_learning_ex dlib)




# Since there are a lot of examples I'm going to use a macro to simply this
# CMakeLists.txt file.  However, usually you will create only one executable in
# your cmake projects and use the syntax shown above.
MACRO(add_example name)
    ADD_EXECUTABLE(${name} ${name}.cpp)
   TARGET_LINK_LIBRARIES(${name} dlib )
   
   # add_definitions(-std=c++11 -g) #to add the -std=c++11 flag 
   # project(app_project)
   # file(GLOB all_src1 "../../simple_apx_app/src/CSrc/*.cpp")
   # file(GLOB all_src2 "../../simple_apx_app/src/Operators/*.cpp")
   # file(GLOB all_src3  "../../simple_apx_app/src/BookKeeping/*.cpp") #includes all the files in the folder src/
   # add_executable(${name} ${name}.cpp ${all_src1} ${all_src2} ${all_src3}) #which files to build
#install(TARGETS tool_exe DESTINATION bin) #where to place the binary
# include_directories(../../simple_apx_app/inc/BookKeeping ../../simple_apx_app/inc/OperatorHeaders)  #adding the include folder
ENDMACRO()

# The deep learning toolkit requires a C++11 capable compiler.  
if (COMPILER_CAN_DO_CPP_11)
   add_example(dnn_introduction_ex)
   add_example(dnn_introduction2_ex)
   add_example(dnn_inception_ex)
   add_example(dnn_imagenet_ex)
   add_example(dnn_imagenet_train_ex)
endif()


SET(GCC_COVERAGE_COMPILE_FLAGS "-std=c++11")
add_definitions(${GCC_COVERAGE_COMPILE_FLAGS})

#
#here we apply our macros 
#add_example(3d_point_cloud_ex)
#add_example(bayes_net_ex)
#add_example(bayes_net_from_disk_ex)
#add_example(bayes_net_gui_ex)
#add_example(bridge_ex)
#add_example(bsp_ex)
#add_example(compress_stream_ex)
#add_example(config_reader_ex)
#add_example(custom_trainer_ex)
#add_example(dir_nav_ex)
#add_example(empirical_kernel_map_ex)
#add_example(face_detection_ex)
#add_example(face_landmark_detection_ex)
#add_example(fhog_ex)
#add_example(fhog_object_detector_ex)
#add_example(file_to_code_ex)
#add_example(graph_labeling_ex)
#add_example(gui_api_ex)
#add_example(hough_transform_ex)
#add_example(image_ex)
#add_example(integrate_function_adapt_simp_ex)
#add_example(iosockstream_ex)
#add_example(kcentroid_ex)
#add_example(kkmeans_ex)
#add_example(krls_ex)
#add_example(krls_filter_ex)
#add_example(krr_classification_ex)
#add_example(krr_regression_ex)
#add_example(learning_to_track_ex)
#add_example(least_squares_ex)
#add_example(linear_manifold_regularizer_ex)
#add_example(logger_custom_output_ex)
#add_example(logger_ex)
#add_example(logger_ex_2)
#add_example(matrix_ex)
#add_example(matrix_expressions_ex)
#add_example(max_cost_assignment_ex)
#add_example(member_function_pointer_ex)
#add_example(mlp_ex)
#add_example(model_selection_ex)
#add_example(mpc_ex)
#add_example(multiclass_classification_ex)
#add_example(multithreaded_object_ex)
#add_example(object_detector_advanced_ex)
#add_example(object_detector_ex)
#add_example(one_class_classifiers_ex)
#add_example(optimization_ex)
#add_example(parallel_for_ex)
#add_example(pipe_ex)
#add_example(pipe_ex_2)
#add_example(quantum_computing_ex)
#add_example(queue_ex)
#add_example(rank_features_ex)
#add_example(running_stats_ex)
#add_example(rvm_ex)
#add_example(rvm_regression_ex)
#add_example(sequence_labeler_ex)
#add_example(sequence_segmenter_ex)
#add_example(server_http_ex)
#add_example(server_iostream_ex)
#add_example(sockets_ex)
#add_example(sockstreambuf_ex)
#add_example(std_allocator_ex)
#add_example(surf_ex)
##add_example(svm_c_ex)
#add_example(svm_ex)
#add_example(svm_pegasos_ex)
#add_example(svm_rank_ex)
#add_example(svm_sparse_ex)
#add_example(svm_struct_ex)
#add_example(svr_ex)
#add_example(thread_function_ex)
#add_example(thread_pool_ex)
#add_example(threaded_object_ex)
#add_example(threads_ex)
#add_example(timer_ex)
#add_example(train_object_detector)
#add_example(train_shape_predictor_ex)
#add_example(using_custom_kernels_ex)
#add_example(video_tracking_ex)
#add_example(xml_parser_ex)
#

find_package(OpenCV QUIET)
if (OpenCV_FOUND)
   include_directories(${OpenCV_INCLUDE_DIRS})

   ADD_EXECUTABLE(webcam_face_pose_ex webcam_face_pose_ex.cpp)
   TARGET_LINK_LIBRARIES(webcam_face_pose_ex dlib ${OpenCV_LIBS} )
else()
   message("OpenCV not found, so we won't build the webcam_face_pose_ex example.")
endif()


if (DLIB_LINK_WITH_SQLITE3)
   add_example(sqlite_ex)
endif()



set_property(GLOBAL PROPERTY RULE_MESSAGES OFF)
cmake_minimum_required(VERSION 2.8)

add_definitions(-std=c++11 -g) #to add the -std=c++11 flag 
project(app_project)
file(GLOB all_src1 "../../simple_apx_app/src/CSrc/*.cpp")
file(GLOB all_src2 "../../simple_apx_app/src/Operators/*.cpp")
file(GLOB all_src3  "../../simple_apx_app/src/BookKeeping/*.cpp") #includes all the files in the folder src/
add_executable(svm_c_ex_2 ${all_src1} ${all_src2} ${all_src3} svm_c_ex.cpp) #which files to build
#install(TARGETS tool_exe DESTINATION bin) #where to place the binary
 TARGET_LINK_LIBRARIES(svm_c_ex_2 dlib )
add_executable(multiclass_classification_ex_2 ${all_src1} ${all_src2} ${all_src3} multiclass_classification_ex.cpp)
TARGET_LINK_LIBRARIES(multiclass_classification_ex_2 dlib)
include_directories(../../simple_apx_app/inc/BookKeeping ../../simple_apx_app/inc/OperatorHeaders)  #adding the include folder

