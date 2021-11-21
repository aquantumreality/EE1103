/**********************************************************************************
* Name: Sai Shashank GP & Abhiram Rao G                                           *
* Date: 03-02-2021                                                                *
* Group: 2B                                                                       *
* Purpose: To calculate approxiate value of pi using random numbers concept       *
* Input: A positive integer(preferably large)                                     *
* Output: Approximate value of pi                                                 *
**********************************************************************************/

#include <stdio.h>     // Standard input output librabry
#include <stdlib.h>    // Standard library
#include <time.h>      // To use time() function
#include <math.h>      // To use acos() function
#include <ctype.h>     // To use RAND_MAX value

double ran_expo(double lambda){
    double u;

    u = rand() / (RAND_MAX + 1.0);

    return -log(1- u) / lambda;
}

int main(int argc, char **argv)  // Two arguments inside the main function
{
    srand(time(0));   // Produces different random number every instant
    int N, i, z;
    float x, y, p;
    char type;
    z = 0;   // Initialising variable z
    N = atoi(argv[1]);   // Taking command line argument input from user and assigning it to N using argv[] and atoi()
    type = argv[2][0];
    for(i = 1; i < N + 1; i ++)    // For loop to check all the N coordinates whether they are inside the circle or not
    {
        if (type == 'e') 
        {
            x = (2.0*ran_expo(2)) - 1;    // Assigning random value to x between 0 and 1
            y = (2.0*ran_expo(2)) - 1;    // Assigning random value to y between 0 and 1
        }
        else if (type == 'n')
        {
            x = (2.0*rand()/RAND_MAX) - 1;
            y = (2.0*rand()/RAND_MAX) - 1;
        }
        
        if(x*x + y*y <= 1)      // Checking whether inside the circle or not
        {
            z = z + 1;      // Keeping count of coordinates inside the circle
        }
        else
        {
            continue;   // Telling computer to do nothing if coordinate is outside the circle
        }
    }
     p = (float) z/N;     // calculating probability of the coordinate falling inside circle which is pi/4
    printf("Approximate value of pi is %f\n", 4.0 * p);   // printing out value of pi obtained by the program
    float PI = acos(-1);            // Initialising value pf pi using acos()
    float e = fabs((4.0 * p) - PI);   // calculating error in our calculated pi
    printf("Error is %f", e);     // Printing error in pi 
    return 0;
}