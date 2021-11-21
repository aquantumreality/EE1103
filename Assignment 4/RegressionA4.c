//************************************************************************************************
//Name: Sai Shashank GP & Abhiram Rao G
//Date: 24-02-2021
//Group: 2B
//Purpose: To find a suitable linear fit to the low and high voltage parts of the graph.(C Program to find m and c for a straight line given)
//Input: Data from the datafile(reading data from the csv file and reading it into an array).
//Output: The approximate linear fits to the I-V graph of a MOSFET along with the value of the threshold voltage (Vth). 
//************************************************************************************************
#include <stdio.h> // Standard input output library

struct Line { double m, c; }; // creating a data .structure and declaring structure variables using {}

void readData(double m[12][401]); // declaring a 2D array of sixe 12x401
struct Line bestApproximate(double x[], double y[], int a, int b); // creates a new variable of type struct Line

//Driver main function 
int main() {
	double m[12][401], Vth[12]; // declaring the Vth array 
	struct Line l1[12], l2[12]; 
	readData(m);
	FILE *fp; // declaring the file pointer
	fp = fopen("output.dat", "w"); // opening an existing file, "w" searches file.
	// If the file exists, its contents are overwritten. If the file doesn’t exist, a new file is created. 
 
	for (int i = 1; i < 12; i++) {
		l1[i] = bestApproximate(m[0], m[i], 0, 7); 
		l2[i] = bestApproximate(m[0], m[i], 133, 400); 
		Vth[i] = (l2[i].c - l1[i].c) / (l1[i].m - l2[i].m); // array for calculating Vth for each value of given Gate Voltage
		printf("y = %ex + %e\n", l1[i].m, l1[i].c); // printing out the linear fit for low voltages
		printf("y = %ex + %e\n", l2[i].m, l2[i].c); // printing out the linear fit for high coltages
		printf("Vth = %e\n", Vth[i]); // printing out the Threshold Voltage
		fprintf(fp, "%d %f\n", i, Vth[i]); // writing to a file
	}

	return 0; 
}

void readData(double m[12][401]) {
	FILE *fp; // declaring the file pointer
	fp = fopen("sgfet.csv", "r"); // opening an exisitng csv file, "r" Searches file.
	// If the file is opened successfully fopen( ) loads it into memory and sets up a pointer which points to the first character in it.
	char waste[100];
	fscanf(fp, "%s", waste);
	for (int i = 0; i < 401; i++) {
		fscanf(fp, "%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf", &m[0][i]
		, &m[1][i], &m[2][i], &m[3][i], &m[4][i], &m[5][i], &m[6][i]
		, &m[7][i], &m[8][i], &m[9][i], &m[10][i], &m[11][i]); // fscanf reads from the csv file
	}
	/*for (int i = 0; i < 401; i++) {
		for (int j = 0; j < 12; j++)
			printf("%e, ", m[i][j]);
		printf("\n");
	}*/
}

// function to calculate m and c that best fit points 
// represented by x[] and y[] 
struct Line bestApproximate(double x[], double y[], int a, int b) 
{ 
	int n = b - a;
	double sum_x = 0, sum_y = 0, sum_xy = 0, sum_x2 = 0; // initialising the variables for computing the required summations
	for (int i = a; i < b; i++) { 
		sum_x += x[i]; 
		sum_y += y[i]; 
		sum_xy += x[i] * y[i]; 
		sum_x2 += (x[i] * x[i]); 
	} // for loop to calculate the summation

	struct Line bestLine;
	bestLine.m = (n * sum_xy - sum_x * sum_y) / (n * sum_x2 - (sum_x * sum_x)); // calculating the required slope
	bestLine.c = (sum_y - bestLine.m * sum_x) / n; // calculating the required y-intercept

	return bestLine; 
} 