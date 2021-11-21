/*****************************************************
 
 Name: Sai Shashank GP and Abhiram Rao G
 Date: 10 - 02 - 2021
 Group: 2B
 Purpose: To isolate and examine a gaussian noise in sine wave
 Inputs: Amplitude(float), Frequency(float), No.of cycles(Integer)
 Outputs: Mean, Standard deviation of gaussian noise and percentages of observations that lie within 2sigma, 4sigma, 6sigma
 
 *****************************************************/

#include <stdio.h>      // Standard input output library
#include <stdlib.h>     // Standard library
#include <time.h>       // To use time() function
#include <math.h>       // To use mathematical functions
#include <ctype.h>

int main(int argc, char **argv)      // Two arguments inside the main() function
{
    float A = atof(argv[1]);         // Taking input of amplitude
    float freq =  atof(argv[2]);     // Taking input of frequency
    int cycles = atoi(argv[3]);      // Taking input of no.of cycles
    float noise_mean;
    float noise_stdev;
    int p = 6;                       // Initialising no.of points per cycle to be calculated
    int N = cycles*p;
    float x[N], y[N], z[N], ERR[N];
    float c = 0.0; 
    float d = 0.0;
    float e = 0.0;
    
    float mean, stdev, err, variance;
    srand(time(0));                    // Produces a random number every instant
    for (int i=0; i < N; i++)
    {
        float U1 = (float) rand()/RAND_MAX;    // U1 is a random number in (0, 1)
        float U2 = (float) rand()/RAND_MAX;    // U2 is a random number in (0, 1)
        x[i] = A * sin(i*2*M_PI*freq/p);       // Defining x[i] as a sine wave array
        z[i] = 0.2 * A * sqrt(-2.0 * log(U1)) * cos(2.0 * M_PI * U2);    // Defining z[i] as a gaussian noise using box muller transformation
        y[i] = x[i] + z[i];         // Defining y[i] as a noisy sine wave
    }
    for (int j=0; j < N; j++)
    {
        ERR[j] = y[j] - x[j];       // Calculating error function
        err += ERR[j];              // summation of all errors
        variance += ERR[j]*ERR[j];  // summation of squares of all errors
    }
    mean = err/N;      // calcuating mean
    variance = (variance/N)-(mean*mean);   // calculating variance
    stdev = sqrt(variance);       // calculating standard deviation
    printf("mean of error for mu=0 = %f\n", mean);    // Printing mean error for default mu
    printf("standard deviation for stdev = 0.2 = %f\n", stdev);  // Printing standard deviation of error for default stdev
    for (int k=0; k < N; k++)
    {
        if((fabs(ERR[k]-mean)<= 1.0*stdev))      // Checking whether given error is within 2 sigma
        c++;
        if((fabs(ERR[k]-mean)<= 2.0*stdev))      // Checking whether given error is within 4 sigma
        d++;
        if((fabs(ERR[k]-mean)<= 3.0*stdev))      // Checking whether given error is within 6 sigma
        e++;
    }
    printf("Points within 1 sigma:%f\n", c*100.0/N);     // Printing percentages of errors that lie within 2 sigma
    printf("Points within 2 sigma:%f\n", d*100.0/N);     // Printing percentages of errors that lie within 4 sigma
    printf("Points within 3 sigma:%f\n", e*100.0/N);     // Printing percentages of errors that lie within 6 sigma

    c = 0;
    d = 0;
    e = 0;

    printf("Enter the value of stdev : ");
    scanf("%f",&noise_stdev);                  // Asking user for fixing noise standard deviation
    printf("Enter the value of mean : ");
    scanf("%f",&noise_mean);                   // Asking user for fixing noise mean

    for (int i=0; i < N; i++)
    {
        float U1 = (float) rand()/RAND_MAX;     // U1 is a random number in (0, 1)
        float U2 = (float) rand()/RAND_MAX;     // U2 is a random number in (0, 1)
        z[i] = 2 * noise_stdev * A * sqrt(-2.0 * log(U1)) * cos(2.0 * M_PI * U2) + 2 * noise_mean * A;     // Defining z[i] as a gaussian noise using box muller tram=nsformation
        y[i] = x[i] + z[i];      // Defining y[i] as a noisy sine wave
    }
    for (int j=0; j < N; j++)
    {
        ERR[j] = y[j] - x[j];     // Calculating error function
        err += ERR[j];            // summation of all errors
        variance += ERR[j]*ERR[j];// summation of squares of all errors
    }
    mean = err/N;        // calcuating mean
    variance = (variance/N)-(mean*mean);    // calculating variance
    stdev = sqrt(variance);      // calculating standard deviation
    printf("mean of error for mu=%f = %f\n",noise_mean,mean);     // Printing mean error for Desired mu
    printf("standard deviation for stdev = %f = %f\n",noise_stdev,stdev);   // Printing standard deviation of error for Desired stdev
    for (int k=0; k < N; k++)
    {
        if((fabs(ERR[k]-mean)<= 1.0*stdev))         // Checking whether given error is within 2 sigma
        c++;
        if((fabs(ERR[k]-mean)<= 2.0*stdev))         // Checking whether given error is within 4 sigma
        d++;
        if((fabs(ERR[k]-mean)<= 3.0*stdev))         // Checking whether given error is within 6 sigma
        e++;
    }
    printf("Points within 1 sigma:%f\n", c*100.0/N);         // Printing percentages of errors that lie within 2 sigma
    printf("Points within 2 sigma:%f\n", d*100.0/N);         // Printing percentages of errors that lie within 4 sigma
    printf("Points within 3 sigma:%f\n", e*100.0/N);         // Printing percentages of errors that lie within 6 sigma
    
    return 0;
}
