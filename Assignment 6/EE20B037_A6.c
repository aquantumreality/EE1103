// Date : 14 March 2021
// Name : G Abhiram (EE20B037)
// Purpose: To find the solution of Landau-Lifshitz-Gilbert equation, an ordinary differential equation.
// Input : Step size (h) , r (radius parameter for spherical coordinates)
/* Output(s):  
             theta_error_vs_iteration.txt - stores values of error in theta, no. of iteration (for the first part of the problem)
             rms_vs_step_size.txt - stores values of RMS error for varying step size (for the second part of the problem)
             switchtime_vs_alpha.txt - stores values of switch time for varying alpha (for the third part of the problem)
             theta_phi.txt - stores values of theta and phi using the gold standard
             ode.txt - stores values of the x, y and z-components of reduced magnetisation vector, i.e. mx, my and mz  (for the fourth part of the problem) */

#include <stdio.h>  // Header file for Standard Input Output
#include <stdlib.h> // Header file for Standard Library
#include <math.h>   // Header file for mathematical functions like pow etc.

#define gamma -1.76e11 // Gyromagnetic ratio in the Landau-Lifshitz-Gilbert Equation
#define pi 3.14159


double Hz = pow(10,6)/4*pi;
double alpha, currentvalue, value; // initializing variables alpha, currentvalue and value

// defining functions for the first derivative of theta and phi
double theta_dot(double theta, double alpha); 
double phi_dot(double theta, double alpha);

double nextvalue_RK45(double currentvalue, double alpha, float h, int flag); /* for solving ODE using the Runge-Kutta-Fehlberg Method, but
in the current problem we're using RK4 method as our standard method */ 
double nextvalue_RK4(double currentvalue, double alpha, float h, int flag); // for solving ODE using the fourth order Runge Kutta RK4 Method
double func_dot(double value, double alpha, int flag); 
float h, r; // initializing h (step size) and r (radius parameter for spherical coordinate system)
double mx, my, mz; // initializing variables for the x, y and z-components of reduced magnetisation vector
void run_for_a_specific_step(float h, float r); 
void run_for_different_steps(); 
void run_for_different_alphas(float h, float r);
void run_for_alpha_and_step(float h, float r); 
void iterate_withstep(double alpha, float h, int flag, float r);

int main(int argc, char **argv){  // the main function of the program calls the functions as needed.

     const float h = atof(argv[1]);
     const float r = atof(argv[2]);
     run_for_a_specific_step(h,r);
     run_for_different_steps();
     run_for_different_alphas(h,r);
     run_for_alpha_and_step(h,r);
     return 0;
}

// Except for part 3 where we are supposed to plot the variation of switching time vs alpha, generally alpha = 0.1 in all other cases.

/* Purpose : To calculate the error in theta for a specific stepsize and store it in theta_error_vs_iteration.txt (Part 1)
   Input : Step size, r */

void run_for_a_specific_step(float h,float r) {
    int i;
    remove("theta_error_vs_iteration.txt");
    FILE *file1; // declaring the file pointer
    file1 = fopen("theta_error_vs_iteration.txt", "a"); /* “a” – Searches file. If the file is opened successfully fopen( ) loads it into memory and 
    sets up a pointer that points to the last character in it. If the file doesn’t exist, a new file is created. */
    fprintf(file1, "iter, error\n"); // Writing to a file
    fclose(file1);
    iterate_withstep(0.1, h, 3, r); // alpha = 0.1 in this case
}

/* Purpose : To calculate RMS error for each step size by varying the stepsize and store it in rms_vs_step_size.txt (Part 2) (No Inputs)*/

void run_for_different_steps() {
    int i;
    remove("rms_vs_step_size.txt");
    FILE *file0; // declaring the file pointer
    file0 = fopen("rms_vs_step_size.txt", "a"); /* “a” – Searches file. If the file is opened successfully fopen( ) loads it into memory and 
    sets up a pointer that points to the last character in it. If the file doesn’t exist, a new file is created. */
    fprintf(file0, "step_size, root_mean_squared_error\n"); // Writing to a file
    fclose(file0); // Closing the file 
    for(i = 5; i < 501 ; i=i+5) {
        iterate_withstep(0.1, i * 0.0000000000000000010, 0,r); // in this case alpha = 0.1, step size is varied as multiples of 10^-18.
        // break;
    }
}
/* Purpose : It runs the iterate_with_step_size function and stores values of switch time
            for different values of alpha in switch_time_vs_alpha.txt file (Part 3)
   Input :  Step size, r */
void run_for_different_alphas(float h, float r) {
    int i;
    remove("switchtime_vs_alpha.txt");
    FILE *file0; // declaring the file pointer
    file0 = fopen("switchtime_vs_alpha.txt", "a"); 
    fprintf(file0, "switchtime, alpha\n"); // Writing to a file
    fclose(file0); // Closing the file
    for(i = 1; i < 100 ; i++) {
       iterate_withstep(0.005 * i, h, 1, r); // to observe variation of switching time with alpha, alpha is varied as multiples of 5x10^-3. 
        // break;
    }
}
/* Purpose: It runs the iterate_with_step_size function, stores values of theta and phi 
            in theta_phi.txt file and values of mx, my and mz in ode.txt file.
   Inputs : Step size, r */
void run_for_alpha_and_step(float h, float r) {
    remove("theta_phi.txt");
    FILE *file0,*file1; // declaring the file pointer(s)
    file0 = fopen("theta_phi.txt", "a");
    file1 = fopen("ode.txt","a");
    fprintf(file0, "iter, theta, phi\n"); // Writing to a file
    fprintf(file1,"mx, my, mz\n");
    fclose(file0); // Closing the file
    fclose(file1);
    iterate_withstep(0.1, h, 2, r);
}
/* Purpose: It runs get_next_value_RK4 (or get_next_value_RKF45 functions) N times and finds the error in theta(for different no. of iterations),
            RMS error(for varying step size), switching time(for different alpha) and the data is stored in respective files.
   Inputs : alpha, step value and flag  */
void iterate_withstep(double alpha, float h, int flag, float r) 
{
    double theta, theta_est, phi, phi_est, theta_error, phi_error, theta_final;
    theta = theta_est = 179.0/180 * 3.14159;
    theta_final = 1.0/180 * 3.14159;
    const int N = 5000000;
    // const long long int N = (theta_org - theta_final)/step;
     //printf("N = %lld\n", N);
    phi = phi_est = 1.0/180 * 3.14159;
    double rms_theta_error = 0, rms_phi_error = 0;
    int i, count = 1;
    double switch_time = 0;
    FILE *file3, *file1, *file2;
    file3 = fopen("theta_phi_values.txt", "a");
    file1 = fopen("theta_error_vs_iteration.txt", "a");
    file2 = fopen("ode.txt","a");
    for(i = 1; i < N; i++) 
    {
        double x = theta ; double y = phi;          // here we are finding error using RK4 method and Euler's method
        theta = nextvalue_RK4(x, alpha, h, 1);      // RK4 Method which is taken as the gold standard and true value
        theta_est = x + h * theta_dot(x,alpha);     // Euler's method for finding the estimated value of theta
        //theta_est =  get_next_value_RKF45(x, alpha, step, 1);  // here we have commented out the part where error is found using
        phi = nextvalue_RK4(y, alpha, h, 0);         // RK4 Method which is taken as the gold standard and true value
        phi_est= y + h * phi_dot(y,alpha);          // Euler's method for finding the estimated value of phi
        
        // calculating error in theta, error in phi when using Euler's method as compared to RK4 method.
        // Then we use the same errors to calculate the RMS error in theta, phi.  
        theta_error = theta - theta_est;
        rms_theta_error += pow(theta_error,2);
        phi_error = phi - phi_est;
        rms_phi_error += pow(phi_error,2);
        // here we are defining (mx,my,mz) the three components of magnetisation vector using 3D spherical coordinates,
        // inputs for computing mx,my,mz : r , theta and phi
        double mx = r * sin (theta) * cos (phi);
        double my = r * sin (theta) * sin (phi);
        double mz = r * cos (theta);
       
        if (flag == 3) {
            if(theta_error<0)
                theta_error=-1.0*theta_error;
            fprintf(file1, "%d, ", i);
            fprintf(file1, "%0.20f\n", (theta_error));
        }
        if(flag == 2) {
            fprintf(file3, "%d, %f, %f\n", i, theta, phi);
            fprintf(file2,"%f,%f,%f\n",mx,my,mz); // writing mx, my, mz to the ode.txt file
        }
        if(cos(theta) >= 0 && count)  // finding the switch time i.e. time when z-component of magnetization changes sign or theta becomes 90 degrees.
         {
            switch_time = i;          // i = number of iterations
            count --;
         }    
    if(theta < 0 || theta > 6.28318)  // breaking the loop when theta crosses 0 or two*pi, since theta(in radians) can take values lying in between [0,6.28318]
        break;
    }
    fclose(file1);
    fclose(file3);
    switch_time = switch_time * h;   // finding the switching time
    if(flag)
    printf("alpha = %0.3f, i = %d,last theta = %lf, %0.32f\n", alpha, i, theta, switch_time);
    rms_theta_error = pow(rms_theta_error/(i), 0.5);   // calculating rms error
    rms_phi_error = pow(rms_phi_error/(i), 0.5);
    
    // exit;
    if(flag == 0) {
        FILE *file0;
        file0 = fopen("rms_vs_step_size.txt", "a");   // storing step and rms error 
        printf("step = %0.20f, rms = %0.9f\n", h, rms_theta_error);
        fprintf(file0, "%0.20f, ", h);
        fprintf(file0, "%0.9f\n", rms_theta_error);
        fclose(file0);
    }

    if (flag == 1) {
        FILE *file2;
        file2 = fopen("switchtime_vs_alpha.txt", "a");  // storing switching time and alpha 
        fprintf(file2, "%0.20f, ", switch_time);
        fprintf(file2, "%0.9f\n", alpha);
        fclose(file2);
    }
}
// function to implement the RK45 method, but we used RK4 method as our gold standard.
double nextvalue_RK45(double currentvalue, double alpha, float h, int flag) {
    double k1, k2, k3, k4, k5, k6;
    k1 = func_dot(currentvalue, alpha, flag);
    k2 = func_dot(currentvalue + (0.2 * k1), alpha, flag);
    k3 = func_dot(currentvalue + ((0.075*k1) + (0.225*k2)), alpha, flag);
    k4 = func_dot(currentvalue + ((0.3*k1) - (0.9*k2) + (1.2*k3)), alpha, flag);
    k5 = func_dot(currentvalue - ((0.2037)*k1) + ((2.5)*k2) - ((2.5926)*k3) + ((1.2963)*k4), alpha, flag);
    k6 = func_dot(currentvalue + ((0.0295)*k1) + ((0.3418)*k2) + ((0.0416)*k3) + ((0.4003)*k4) + ((0.0618)*k5), alpha, flag);
    return currentvalue + h*((0.0979)*k1 + (0.4025)*k3 + (0.2104)*k4 + (0.2891)*k6);
}
// function to implement the RK4 method.
double nextvalue_RK4(double currentvalue, double alpha, float h, int flag) {
    double k1, k2, k3, k4, k5, k6;
    k1 = func_dot(currentvalue, alpha, flag);
    k2 = func_dot(currentvalue + 1.0/2 * k1, alpha, flag);
    k3 = func_dot(currentvalue + 1.0/2 * k2, alpha, flag);
    k4 = func_dot(currentvalue + k3, alpha, flag);
    return currentvalue + h/6.0 * (k1 + 2.0 * k2 + 2.0 * k3 + k4); 
}

double func_dot(double value, double alpha, int flag) {
    if(flag) {
        return theta_dot(value, alpha);
    }
    return phi_dot(value, alpha);
}

double theta_dot(double theta, double alpha) {
    return (gamma * Hz * alpha * sin(theta))/(alpha * alpha + 1); // returns the value of theta dot
}

double phi_dot(double theta, double alpha) {
    return (gamma * Hz)/(alpha * alpha + 1); // returns the value of phi dot
}

