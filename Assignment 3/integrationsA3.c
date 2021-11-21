/**
 * Names	:	Sai Shashank GP & Abhiram Rao G
 * Date		:	19-02-2021
 * Purpose	:	This program generates a sine wave, adds random gaussian noise to it,
 * 				and tabulates the error in the 3 methods of integration-
 * 				midpoint rule, trapezoid rule, simpson's rule.
 * 				Also compares with the romberg method.
 * Inputs	:	M	: The number of cycles of the sine wave
 * 				N	: The number of points per each cycle
 * Outputs	:	table.dat	: The data of the errors
 */

#include <stdio.h>							// standard io library
#include <math.h>							// standard math library for sin, log, sqrt
#include <stdlib.h>							// standard library for rand
#include <time.h>							// standard time library for srand(time(0))

double noise();								// makes gaussian noise using AMP
void  tabulate(double *x, int N, int T);	// tabulates the results
double midp(double *x, double h, int T);	// returns error of integration through mid-point
double trap(double *x, double h, int T);	// returns error of integration through trapezoid
double simp(double *x, double h, int T);	// returns error of integration through simpson's 
double berg(double *x, int n, int T);		// the mystical romberg method

/* Purpose	:	This is the main method of the program.
 * 				It parses the inputs and generates sine wave accordingly,
 * 				Then it calls the other functions as needed
 * Inputs	:	M, N
 * Outputs  :   none
 */
int main(int argc, char **argv) {
	const int M	= atoi( argv[1] );					// number of cycles
	const int N	= atoi( argv[2] );					// points per cycle
	const int T = M * N;							// total points including midpts
	double scale = 2 * 3.141592653589793238 / T;	// scale of x-axis: 2pi / T

	srand( time(0) );
	double *x = malloc( sizeof(double) * T );		// declare x: noisy sine array

	for (int i = 0; i < T; i++)
		x[i] = sin(i * scale * M) + noise();		// generate noisy sine wave

	tabulate(x, N, T);

	free(x);										// conserve resources, free memory
	return 0;
}

/* Purpose	:	This function generates gaussian noise of mean = 0, stdev = temp
 * 				through the Marsaglia Polar implementation of Box-Muller transform
 * Inputs	:	temp: Hard-coded stdev = 0.05
 * Outputs	:	returns gaussian noise as a double
 */
double noise() {
	static double storage;						// stores extra noise, if any
	static int flag = 0;						// keeps track of if noise has been used
	double u, v, r;								// temporary variables
	double temp = 0.05;							// stdev is 0.05

	flag = 1 - flag;							// flip flag (either 0 or 1)

	if (flag == 0) return storage;				// return already generated noise, if it exists
	
	do {										// otherwise generate pair of noises
		u = ( 2.0 * rand() / RAND_MAX ) - 1.0;	// u, v are uniformly distributed 
		v = ( 2.0 * rand() / RAND_MAX ) - 1.0;	// in [-1, 1]
		r = u * u  +  v * v;					// distance from origin
	} while (r >= 1 || r == 0);					// repeat if outside unit circle / at origin
		
	temp *= sqrt( -2.0 * log(r) / r );			// mathemagic
	
	storage = u * temp;							// store one result
	return ( v * temp );						// return other result
}

/* Purpose	:	This function tabulates the data by
 * 				printing to the console as well as to a file
 * 				table.dat for graphing
 * Inputs	:	x: noisy sine wave
 * 				N: points per cycle
 * 				T: total points
 * Outputs	:	(n, M, T, S) to console and table.dat
 */
void tabulate(double *x, int N, int T) {
	FILE *fp;										// we will write to table.dat
	fp = fopen("table.dat", "w");
	for (int n = 4; n <= N; n++) {
		/* romberg is slow for large n that is needed for the other methods, */
		/* so keep berg commented out unless it is to be tabulated currently */
		double b = berg(x, n, T);					// the mysterious romberg's method

		double h = (double) N / n;					// interval width
		double m = midp(x, h, T);
		double t = trap(x, h, T);
		double s = simp(x, h, T);
		
		printf("%d  \t%f\t%f\t%f\t%f\n", n, m, t, s);
		fprintf(fp, "%d  \t%f\t%f\t%f\t%f\n", n, m, t, s);

		if (m < 0.0001 && t < 0.0001 && s < 0.0001)
			break;									// break loop if all converge
	}
	fclose(fp);										// conserve resources, close streams
}

/* Purpose	:	This function calculates the integral of the noisy sine wave
 * 				through the midpoint method and returns the error in it
 * Inputs	:	x: noisy sine wave
 * 				h: width of intervals to take
 * 				T: size of array x
 * Outputs	:	returns error in midpoint method
 */
double midp(double *x, double h, int T) {
	double scale = 2 * 3.141592653589793238 / T;	// scale of x-axis: 2pi / T
	double error = 0;

	for (double i = 0; (int)(i + h) < T; i += h) {
		error += x[ (int)(i + h/2) ];
	}
	error *= scale * h;

	return (error > 0)? error : -error;
}

/* Purpose	:	This function calculates the integral of the noisy sine wave
 * 				through the trapezoid method and returns the error in it
 * Inputs	:	x: noisy sine wave
 * 				h: width of intervals to take
 * 				T: size of array x
 * Outputs	:	returns error in trapezoid method
 */
double trap(double *x, double h, int T) {
	double scale = 2 * 3.141592653589793238 / T;	// scale of x-axis: 2pi / T
	double error = 0;

	for (double i = 0; (int)(i + h) < T; i += h) {
		error += x[ (int)(i) ];
		error += x[ (int)(i + h) ];
	}
	error *= scale * h/2;

	return (error > 0)? error : -error;
}

/* Purpose	:	This function calculates the integral of the noisy sine wave
 * 				through the simpson's method and returns the error in it
 * Inputs	:	x: noisy sine wave
 * 				h: width of intervals to take
 * 				T: size of array x
 * Outputs	:	returns error in simpson's method
 */
double simp(double *x, double h, int T) {
	double scale = 2 * 3.141592653589793238 / T;	// scale of x-axis: 2pi / T
	double error = 0;

	for (double i = 0; (int)(i+h) < T; i += h) {
		error += x[ (int)(i) ];
		error += x[ (int)(i + h/2) ] * 4;
		error += x[ (int)(i + h) ];
	}
	error *= scale * h/6;

	return (error > 0)? error : -error;
}

/* Purpose	:	This is the romberg method of integration, which is slow for
 * 				large number of steps, but gives accurate answer pretty fast
 * 				Adapted from: https://en.wikipedia.org/wiki/Romberg%27s_method
 * Inputs	:	x: noisy sine array
 * 				n: maximum number of steps to be taken
 * 				T: size of array x
 * Outputs	:	returns error in romberg's method
 */
double berg(double *x, int n, int T) {
	double scale = 2 * 3.141592653589793238 / T;		// scale of x-axis: 2pi / T
	double b = (T - 1) * scale;							// upper limit
	double R1[n], R2[n]; 								// buffers
	double *Rp = &R1[0], *Rc = &R2[0]; 					// Rp is previous row, Rc is current row
	Rp[0] = x[T - 1] * b * 0.5;						// first trapezoidal step

	for (size_t i = 1; i < n; ++i) {
		b /= 2.0;
		double c = 0;
		size_t ep = 1 << (i - 1); 						// 2^(n-1)
		for (size_t j = 1; j <= ep; ++j) {
			c += x[ (int)((2*j - 1) * b * scale) ];
		}
		Rc[0] = b * c + 0.5 * Rp[0]; 					// R(i,0)

		for (size_t j = 1; j <= i; ++j) {
			double t = pow(4, j);
			Rc[j] = (t * Rc[j-1] - Rp[j-1]) / (t - 1);	// compute R(i,j)
		}

		// swap Rn and Rc as we only need the last row
		double* rt = Rp;
		Rp = Rc;
		Rc = rt;
	}

	double temp = Rp[n - 1];
	return (temp > 0)? temp : -temp;					// return our positive error
}
