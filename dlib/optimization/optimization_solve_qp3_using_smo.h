// Copyright (C) 2007  Davis E. King (davis@dlib.net)
// License: Boost Software License   See LICENSE.txt for the full license.
#ifndef DLIB_SOLVE_QP3_USING_SMo_Hh_
#define DLIB_SOLVE_QP3_USING_SMo_Hh_

#include "optimization_solve_qp3_using_smo_abstract.h"
#include <cmath>
#include <limits>
#include <sstream>
#include "../matrix.h"
#include "../algs.h"

//copy paste the following lines: start here
#include <fstream>
#include "assert.h"
#include "Operators.h"
#include "operatorFile_parser.h"
#include "setSubType.h"
#include "operandFile_parser.h"
#include "globals.h"
using namespace std;
//#include "foo.h"
extern hw_ac **myOp;   
// end here
#include <typeinfo>
//#include<iostream>
namespace dlib 
{

// ----------------------------------------------------------------------------------------

    class invalid_qp3_error : public dlib::error 
    { 

    public: 
        invalid_qp3_error(
            const std::string& msg, 
            double B_,
            double Cp_,
            double Cn_
        ) : 
            dlib::error(msg), 
            B(B_),
            Cp(Cp_),
            Cn(Cn_)
        {};

        const double B;
        const double Cp;
        const double Cn;
    };

// ----------------------------------------------------------------------------------------

    template <
        typename matrix_type
        >
    class solve_qp3_using_smo
    {
    public:
        typedef typename matrix_type::mem_manager_type mem_manager_type;
        typedef typename matrix_type::type scalar_type;
        typedef typename matrix_type::layout_type layout_type;
        typedef matrix<scalar_type,0,0,mem_manager_type,layout_type> general_matrix;
        typedef matrix<scalar_type,0,1,mem_manager_type,layout_type> column_matrix;


        template <
            typename EXP1,
            typename EXP2,
            typename EXP3,
            long NR
            >
        unsigned long operator() ( 
            const matrix_exp<EXP1>& Q,
            const matrix_exp<EXP2>& p,
            const matrix_exp<EXP3>& y,
            const scalar_type B,
            const scalar_type Cp,
            const scalar_type Cn,
            matrix<scalar_type,NR,1,mem_manager_type, layout_type>& alpha,
            scalar_type eps
        ) 
        {
            DLIB_ASSERT(Q.nr() == Q.nc() && y.size() == Q.nr() && p.size() == y.size() && 
                        y.size() > 0 && is_col_vector(y) && is_col_vector(p) &&
                        sum((y == +1) + (y == -1)) == y.size() &&
                        Cp > 0 && Cn > 0 &&
                        eps > 0,
                "\t void solve_qp3_using_smo::operator()"
                << "\n\t invalid arguments were given to this function"
                << "\n\t Q.nr():                     " << Q.nr() 
                << "\n\t Q.nc():                     " << Q.nc() 
                << "\n\t is_col_vector(p):           " << is_col_vector(p) 
                << "\n\t p.size():                   " << p.size() 
                << "\n\t is_col_vector(y):           " << is_col_vector(y) 
                << "\n\t y.size():                   " << y.size() 
                << "\n\t sum((y == +1) + (y == -1)): " << sum((y == +1) + (y == -1)) 
                << "\n\t Cp:                         " << Cp
                << "\n\t Cn:                         " << Cn
                << "\n\t eps:                        " << eps 
                );



            set_initial_alpha(y, B, Cp, Cn, alpha);


            const scalar_type tau = 1e-12;

            typedef typename colm_exp<EXP1>::type col_type;

            // initialize df.  Compute df = Q*alpha + p
            df = p;
            for (long r = 0; r < df.nr(); ++r)
            {
                if (alpha(r) != 0)
                {
                   //out << typeid(matrix_cast<scalar_type>(colm(Q,r))).name() << endl;
                   df += alpha(r)*matrix_cast<scalar_type>(colm(Q,r)); //ORNIG


                 //  int a = alpha(r);
                  // int b = matrix_cast<scalar_type>(colm(Q,r));                    
                   //int dtmp = mySlOp[0]->calc(a,b);
                   //int c = p;
                   //int dtmp = mySlOp[0]->calc(alpha(r),matrix_cast<scalar_type>(colm(Q,r)));
                   //df = mySlOp[1]->calc(c,dtmp);
                }
            }

            unsigned long count = 0;
            // now perform the actual optimization of alpha
            long i=0, j=0;
            while (find_working_group(y,alpha,Q,df,Cp,Cn,tau,eps,i,j))
            {
                ++count;
                const scalar_type old_alpha_i = alpha(i);
                const scalar_type old_alpha_j = alpha(j);

                optimize_working_pair(alpha,Q,y,df,tau,i,j, Cp, Cn );

                // update the df vector now that we have modified alpha(i) and alpha(j)
           //       scalar_type delta_alpha_i = alpha(i) - old_alpha_i;   //ORIG (clean)
               
                  //cout << typeid(old_alpha_i).name() << endl;
              
                double a = alpha(i);
                double b = (-1)*old_alpha_i;
             //   printf("Op0");
                  scalar_type delta_alpha_i = myOp[0]->calc(a,b);   //Addition (subtraction)
                


  //              scalar_type delta_alpha_j = alpha(j) - old_alpha_j; //ORIG (clean)
                a = alpha(j);
                b = (-1)*old_alpha_j;
//                printf("Op1\n");
                 scalar_type delta_alpha_j = myOp[1]->calc(a,b);    //Addition(subtraction)



                col_type Q_i = colm(Q,i);
                col_type Q_j = colm(Q,j);
                for(long k = 0; k < df.nr(); ++k){
//                       df(k) += Q_i(k)*delta_alpha_i + Q_j(k)*delta_alpha_j;//ORIG
//                    printf("Q_i(k) is: ");                  
//                    cout << typeid(Q_i(k)).name() << endl;
//                    printf("delta_alpha_i is: ");
//                    cout << typeid(delta_alpha_i).name() << endl;
//                    printf("Q_j(k) is: ");
//                    cout << typeid(Q_j(k)).name() << endl;
//                    printf("delta_alpha_j is: ");
//                    cout << typeid(delta_alpha_j).name() << endl;
//                    printf("df(k) is: ");
//                    cout << typeid(df(k)).name() << endl;
//
 //                   a = delta_alpha_i*Q_i(k);
  //                  b = delta_alpha_j*Q_j(k);
                    double c = Q_i(k);
                    a = delta_alpha_i;
//                    printf("Op2\n");
                    a = myOp[2]->calc(c,a);                           //Multiplicaton
                    c = Q_j(k);
                    b = delta_alpha_j;
//                    printf("Op3\n");
                    b = myOp[3]->calc(c,b);                          //Multiplication
 //                   printf("Op4\n");
                    a = myOp[4]->calc(a,b);                          //Addition
                    b = df(k);
  //                  printf("Op5\n");
                    df(k)= myOp[5]->calc(a,b);                        //Addition
//
//                 

//                    df(k) += a + b;

                }


            }

            return count;
        }

        const column_matrix& get_gradient (
        ) const { return df; }

    private:

    // -------------------------------------------------------------------------------------

        template <
            typename scalar_type,
            typename scalar_vector_type,
            typename scalar_vector_type2
            >
        inline void set_initial_alpha (
            const scalar_vector_type& y,
            const scalar_type B,
            const scalar_type Cp,
            const scalar_type Cn,
            scalar_vector_type2& alpha
        ) const
        {
            alpha.set_size(y.size());

            set_all_elements(alpha,0);

            // It's easy in the B == 0 case
            if (B == 0)
                return;

            const scalar_type C = (B > 0)?  Cp : Cn;

            scalar_type temp = std::abs(B)/C;
            long num = (long)std::floor(temp);
            long num_total = (long)std::ceil(temp);

            const scalar_type B_sign = (B > 0)? 1 : -1;

            long count = 0;
            for (long i = 0; i < alpha.nr(); ++i)
            {
                if (y(i) == B_sign)
                {
                    if (count < num)
                    {
                        ++count;
                        alpha(i) = C;
                    }
                    else 
                    {
                        if (count < num_total)
                        {
                            ++count;
                            alpha(i) = C*(temp - std::floor(temp)); //ORNIG
                            
                            printf("alpha(i) is: ");
                            cout << typeid(alpha(i)).name() << endl;
                            printf("C is: ");
                            cout << typeid(C).name() << endl;
                            printf("temp is: ");
                            cout << typeid(temp).name() << endl;
                            printf("std::floor (temp) is: ");
                            cout << typeid(std::floor(temp)).name() << endl;

                          

                        }
                        break;
                    }
                }
            }

            if (count != num_total)
            {
                std::ostringstream sout;
                sout << "Invalid QP3 constraint parameters of B: " << B << ", Cp: " << Cp << ", Cn: "<< Cn;
                throw invalid_qp3_error(sout.str(),B,Cp,Cn);
            }
        }

    // ------------------------------------------------------------------------------------

        template <
            typename scalar_vector_type,
            typename scalar_type,
            typename EXP,
            typename U, typename V
            >
        inline bool find_working_group (
            const V& y,
            const U& alpha,
            const matrix_exp<EXP>& Q,
            const scalar_vector_type& df,
            const scalar_type Cp,
            const scalar_type Cn,
            const scalar_type tau,
            const scalar_type eps,
            long& i_out,
            long& j_out
        ) const
        {
            using namespace std;

            long ip = 0;
            long jp = 0;


            typedef typename colm_exp<EXP>::type col_type;
            typedef typename diag_exp<EXP>::type diag_type;

            scalar_type ip_val = -numeric_limits<scalar_type>::infinity();
            scalar_type jp_val = numeric_limits<scalar_type>::infinity();

            // loop over the alphas and find the maximum ip and in indices.
            for (long i = 0; i < alpha.nr(); ++i)
            {
                if (y(i) == 1)
                {
                    if (alpha(i) < Cp)
                    {
                        if (-df(i) > ip_val)
                        {
                            ip_val = -df(i);
                            ip = i;
                        }
                    }
                }
                else
                {
                    if (alpha(i) > 0.0)
                    {
                        if (df(i) > ip_val)
                        {
                            ip_val = df(i);
                            ip = i;
                        }
                    }
                }
            }

            scalar_type Mp = -numeric_limits<scalar_type>::infinity();

            // Pick out the column and diagonal of Q we need below.  Doing
            // it this way is faster if Q is actually a symmetric_matrix_cache
            // object.
            col_type Q_ip = colm(Q,ip);
            diag_type Q_diag = diag(Q);



            // now we need to find the minimum jp indices
            for (long j = 0; j < alpha.nr(); ++j)
            {
                if (y(j) == 1)
                {
                    if (alpha(j) > 0.0)
                    {
                        scalar_type b = ip_val + df(j);
                        if (df(j) > Mp)
                            Mp = df(j);

                        if (b > 0)
                        {
                            scalar_type a = Q_ip(ip) + Q_diag(j) - 2*y(ip)*Q_ip(j); 
                            if (a <= 0)
                                a = tau;
                            scalar_type temp = -b*b/a;
                            if (temp < jp_val)
                            {
                                jp_val = temp;
                                jp = j;
                            }
                        }
                    }
                }
                else
                {
                    if (alpha(j) < Cn)
                    {
                        scalar_type b = ip_val - df(j);
                        if (-df(j) > Mp)
                            Mp = -df(j);

                        if (b > 0)
                        {
                            scalar_type a = Q_ip(ip) + Q_diag(j) + 2*y(ip)*Q_ip(j); 
                            if (a <= 0)
                                a = tau;
                            scalar_type temp = -b*b/a;
                            if (temp < jp_val)
                            {
                                jp_val = temp;
                                jp = j;
                            }
                        }
                    }
                }
            }

            // if we are at the optimal point then return false so the caller knows
            // to stop optimizing
            if (Mp + ip_val < eps)
                return false;


            i_out = ip;
            j_out = jp;

            return true;
        }

    // ------------------------------------------------------------------------------------

        template <
            typename EXP,
            typename EXP2,
            typename T, typename U
            >
        inline void optimize_working_pair (
            T& alpha,
            const matrix_exp<EXP>& Q,
            const matrix_exp<EXP2>& y,
            const U& df,
            const scalar_type& tau,
            const long i,
            const long j,
            const scalar_type& Cp,
            const scalar_type& Cn
        ) const
        {
            const scalar_type Ci = (y(i) > 0 )? Cp : Cn;
            const scalar_type Cj = (y(j) > 0 )? Cp : Cn;

            if (y(i) != y(j))
            {
                scalar_type quad_coef = Q(i,i)+Q(j,j)+2*Q(j,i);
                if (quad_coef <= 0)
                    quad_coef = tau;
                scalar_type delta = (-df(i)-df(j))/quad_coef;


//                  scalar_type diff = alpha(i) - alpha(j);  //ORIG (clean) 
//                printf("alpha(i) is: ");
//                cout << typeid(alpha(i)).name() << endl;
//                printf("alpha(j) is: ");
//                cout << typeid(alpha(j)).name() << endl; 
//                printf("diff is: ");
//                cout << typeid(diff).name() << endl;

                  double a = alpha(i);
                  double b = alpha(j)*-1;
              //  printf("Op6\n");
                  scalar_type diff = myOp[6]->calc(a,b); //Addition (subtraction)

//                  alpha(i) += delta; //ORIG (clean)
               // printf("Op7\n");
                  alpha(i) = myOp[7]->calc(alpha(i),delta); //Addition
 //                 alpha(j) += delta; //ORIG (clean)
               // printf("Op8\n");
                  alpha(j) = myOp[8]->calc(alpha(j),delta); //Addition

              //  printf("delta: ");
               // cout << typeid(delta).name() << endl;


                if (diff > 0)
                {
                    if (alpha(j) < 0)
                    {
                        alpha(j) = 0;
                        alpha(i) = diff;

                    }
                }
                else
                {
                    if (alpha(i) < 0)
                    {
                        alpha(i) = 0;
                        
                    
                        
                       alpha(j) = -diff; //ORIG
                       a = alpha(j);
                       b = diff*-1;
                       
//                      alpha(j) = myOp[9]->calc(a,b);  //Addition (subtraction) INCORRECT
                    }
                }

                if (diff > Ci - Cj)
                {
                    if (alpha(i) > Ci)
                    {
                        alpha(i) = Ci;
//                        alpha(j) = Ci - diff; //ORIG

                        b = diff*-1;
                        alpha(j) = myOp[10]->calc(Ci,b); //Addition (subtraction)

                    }
                }
                else
                {
                    if (alpha(j) > Cj)
                    {
                        alpha(j) = Cj;
//                        alpha(i) = Cj + diff; //ORIG
                        
                        alpha(i) = myOp[11]->calc(Cj, diff); //Addition

                    }
                }
            }
            else
            {
       //        scalar_type quad_coef = Q(i,i)+Q(j,j)-2*Q(j,i); //ORIG (clean)
             
             //       printf("Q(i,i) is: ");
             //       cout << typeid(Q(i,i)).name() << endl;                    
             //       printf("Q(j,j) is: ");
             //       cout << typeid(Q(j,j)).name() << endl;
             //       printf("Q(j,i) is: ");
             //       cout << typeid(Q(j,i)).name() << endl;
             //       printf("quad_coef is: ");
             //       cout << typeid(quad_coef).name() << endl;
                    
                float tmpflt = myOp[12]->calc(Q(i,i),Q(j,j)); //Addition
                float tmpflt2 = -2*(Q(j,i));               
                scalar_type quad_coef = myOp[13]->calc(tmpflt, tmpflt2); //Addition (subtraction)


                if (quad_coef <= 0)
                    quad_coef = tau;
                
               

               // scalar_type delta = (df(i)-df(j))/quad_coef; //ORIG (clean)
                
                double a = df(i);
                double b = df(j)*-1;
                a = myOp[14]->calc(a,b); //Addition (subtraction)
                
                b = quad_coef;
                b = 1/b;
             //   scalar_type delta = myOp[15]->calc(a,b); // Multiplication (division);
                  scalar_type delta = (a)/quad_coef;  
                


//                scalar_type sum = alpha(i) + alpha(j);  //ORIG (clean)
                scalar_type sum = myOp[16]->calc(alpha(i),alpha(j)); //Addition


               // alpha(i) -= delta; //ORIG (clean)
                b = delta * -1;
                alpha(i) = myOp[17]->calc(alpha(i),b); //Addition (subtraction)

              //  alpha(j) += delta; //ORIG (clean)
                alpha(j) = myOp[18]->calc(alpha(j),delta); //Addition
                

                if(sum > Ci)
                {
                    if(alpha(i) > Ci)
                    {
                        alpha(i) = Ci;
                      //  alpha(j) = sum - Ci; //ORIG
                        b = Ci * -1;
                        alpha(j) = myOp[19]->calc(sum, b); //Addition (subtraction)
                        
                        

                    }
                }
                else
                {
                    if(alpha(j) < 0)
                    {
                        alpha(j) = 0;
                        alpha(i) = sum;
                    }
                }

                if(sum > Cj)
                {
                    if(alpha(j) > Cj)
                    {
                        alpha(j) = Cj;
                   //     alpha(i) = sum - Cj;//ORIG 
                        b = Cj * -1;
                        alpha(i) = myOp[20]->calc(sum, b); //Addition (subtraction)
                    }
                }
                else
                {
                    if(alpha(i) < 0)
                    {
                        alpha(i) = 0;
                        alpha(j) = sum;
                    }
                }

            }
        }

    // ------------------------------------------------------------------------------------

        column_matrix df; // gradient of f(alpha)
    };

// ----------------------------------------------------------------------------------------

}

#endif // DLIB_SOLVE_QP3_USING_SMo_Hh_


