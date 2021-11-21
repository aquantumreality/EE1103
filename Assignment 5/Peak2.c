#include <conio.h> 
#include <stdio.h> 
#include <string.h> 
#include <stdlib.h>
#include <math.h>

// Driver Code 
int main() 
{ 
	// Substitute the full file path 
	// for the string file_path 
	FILE * fp;
    fp=fopen("data.txt","r");
    float y[9800],x[9800];
    for(int i=0;i<9800;i++)
    {
        fscanf(fp,"%f %f",(y+i),(x+i));
    }
    fclose(fp);
    
	int counter=0,width1=0,Timeperiod1=0,width2=0,Timeperiod2=0,holder,lmax,Peaks[10000],Peaks1[10000];
    double Avg[10000],Var[10000],ULimit[10000],LLimit[10000];
	// Window size = 25 & p = 0.75 
    Avg[12]=0.00;
    Var[12]=0.00;
    for(int i=0;i<25;i++)
    {
        Avg[12]=Avg[12]+(x[i]/25);		// Average of the data in a window
        Var[12]=Var[12]+(pow(x[i],2)/25);		// Variance of the data in a window
    }
    Var[12]=Var[12]-(pow(Avg[12],2));
    for(int k=13;k<=10000;k++)
    {
        Avg[k]=0.00;
        Var[k]=0.00;
        Avg[k]=Avg[k-1]+((x[k+12]-x[k-13])/25);
        Var[k]=Var[k-1]+((pow(x[k+12],2)-pow(x[k-13],2))/25)+(pow(Avg[k-1],2)-pow(Avg[k],2));
        ULimit[k]=Avg[k]+(0.75*(pow(Var[k],0.50)));		// Lower error band
        LLimit[k]=Avg[k]-(0.75*(pow(Var[k],0.50)));		// Higher error band
    }
    for(int l=0;l<13;l++)
    {
        ULimit[l]=ULimit[13];
        LLimit[l]=LLimit[13];
		Avg[l]=Avg[13];
    }
	for(int m=14;m<=10000;m++)		// For finding peaks
	{
		if(x[m-1]>LLimit[m-1] && x[m]<LLimit[m] && x[m+1]<LLimit[m+1] && abs(m-Peaks[counter-1])>10)
		{
			Peaks[counter]=m;
			counter++;
		}
	}
	for(int n=0;n<counter;n++)
	{
		holder=Peaks[n];
		for(int p=Peaks[n]+1;p<Peaks[n]+10;p++)
		{
			if(x[p]<x[holder])
			holder=p;
		}
		Peaks[n]=holder;
	}
	for(int q=0;q<counter;q++)
	{
		lmax=Peaks[q];
		for(int r=Peaks[q]-25;r<Peaks[q]+25;r++)
		{
			if(y[r]>y[lmax])
			lmax=r;
		}
		Peaks1[q]=lmax;
	}
	for(int s=2;s<counter;s++)
	{
		Timeperiod1=Timeperiod1+((Peaks[s]-Peaks[s-2])/(counter-2));		// For column 2
		Timeperiod2=Timeperiod2+((Peaks1[s]-Peaks1[s-2])/(counter-2));		// For column 1
	}
	for(int t=2;t<counter;t+=2)
	{
		width1=width1+(2*(Peaks[t]-Peaks[t-1])/(counter-1));		// For column 2
		width2=width2+(2*(Peaks1[t]-Peaks1[t-1])/(counter-1));		// For column 1
	}
	printf("The time interval between two droplets for column 2 and column 1 data are %d and %d units respectively\n",Timeperiod1,Timeperiod2);
	printf("The width of the droplet for column 2 and column 1 data are %d and %d units respectively",width1,width2);

	return 0;
} 

