//---------------------------------------------------------------------
// findmax.cpp
// Used to find the maximum reaching area of a participant and to write
// maximum values to test.csv
//---------------------------------------------------------------------
#include "HapticAPI2.h"
#include "HapticUtility.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <windows.h>  
#include <cstdlib>
#include <algorithm>
#include <iostream> // cin
#include <fstream>
using namespace std;

// Set path to where you want to save arm weight and workspace data
const char* setuppath = "test.csv";

#define IPADDRESS "10.30.203.36" //"10.30.203.26"

#define PosX 0
#define PosY 1
#define PosZ 2

#define maxForce 110.66
#define support_level -0.5

long dev = 0; 
char response[100];
double CurrentPosition[3];


//---------------------------------------------------------------------
//                            M A I N
//---------------------------------------------------------------------
int main(int argc, char** argv)
{
	int returnValue;
	double CurrentForce[3];
	double springPos[3] = {0.0,0.05,-0.17};
	double springStiffness = 1000.0, springDamping = 10.0, springMaxForce = 100.0;
	double xmin, xmax, ymin, ymax;

	// Open the HapticMaster device
	dev = haOpenDevice(IPADDRESS);

	// Check if device connected, then initialize
	if (dev == HARET_ERROR) {
		printf("--- ERROR: Unable to connect to device: %s\n", IPADDRESS);
		return HARET_ERROR;
	}
	else {
		InitializeDevice(dev);
	}

	// printf("Press ENTER to move to position.\n");
	// std::cin.get();
	printf("Moving...\n");
	returnValue = haSendString(dev, "create spring mySpring", response);
	haSendCommand(dev, "set mySpring pos", springPos[PosX], springPos[PosY], springPos[PosZ], response);
	haSendCommand(dev, "set mySpring stiffness", springStiffness, response);
	haSendCommand(dev, "set mySpring dampfactor", springDamping, response);
	haSendCommand(dev, "set mySpring maxforce", springMaxForce, response);
	returnValue = haSendString(dev, "set mySpring enable", response);

	printf("Press ENTER to weigh.\n");
	std::cin.get();
	printf("Weighing...\n");
	haSendCommand(dev, "get measforce", response);
	ParseFloatVec(response, CurrentForce[PosX], CurrentForce[PosY], CurrentForce[PosZ]);
	double armWeight = -CurrentForce[PosZ];
	printf("Arm weight is %f Newtons.\n", armWeight);
	
	printf("Adjust person in seat and press ENTER to allow movement.\n ");
	cin.get();
	returnValue = haSendString(dev, "set mySpring disable", response);

	printf("Move closest to the BODY, then press ENTER.\n");
	std::cin.get();
	printf("Measuring BODY...\n");
	haSendCommand(dev, "get measpos", response);
	ParseFloatVec(response, CurrentPosition[PosX], CurrentPosition[PosY], CurrentPosition[PosZ]);
	xmax = CurrentPosition[PosX];
	printf("xmax is %f.\n", xmax);

	// Set bias force
	printf("Press ENTER to set bias force.\n");
	std::cin.get();
	printf("Setting bias force to %f of the max abduction force.\n", support_level);
	float bias = maxForce * support_level;
	returnValue = haSendString(dev, "create biasforce myBiasForce", response);
	printf("bias: %f \n", bias);
	returnValue = haSendCommand(dev, "set myBiasForce force", 0.0, 0.0, bias, response);
	haSendCommand(dev, "set myBiasForce enable", response);
	printf("Bias force: %s\n", response);

	printf("Move to RIGHT CORNER, then press ENTER.\n");
	std::cin.get();
	printf("Measuring RIGHT CORNER...\n");
	haSendCommand(dev, "get measpos", response);
	ParseFloatVec(response, CurrentPosition[PosX], CurrentPosition[PosY], CurrentPosition[PosZ]);
	ymax = CurrentPosition[PosY];
	float xmin1 = CurrentPosition[PosX];
	printf("ymax is %f.\n", ymax);

	printf("Move to LEFT CORNER, then press ENTER.\n");
	std::cin.get();
	printf("Measuring LEFT CORNER...\n");
	haSendCommand(dev, "get measpos", response);
	ParseFloatVec(response, CurrentPosition[PosX], CurrentPosition[PosY], CurrentPosition[PosZ]);
	ymin = CurrentPosition[PosY];
	float xmin2 = CurrentPosition[PosX];
	printf("ymin is %f.\n", ymin);

	if (xmin1 > xmin2)
	{
		xmin = xmin1;
	}
	else {
		xmin = xmin2;
	}

	printf("xmin is %f.\n", xmin);

	haSendCommand(dev, "remove all", response);
	printf("remove all ==> %s\n", response);
	haSendCommand(dev, "set state stop", response);
	printf("set state stop ==> %s\n", response);
	if (haCloseDevice(dev)) {
		printf("--- ERROR: closing device\n");
	}
	
	ofstream myfile;
	myfile.open(setuppath);
	//myfile << "weight, ymin, ymax, xmin, xmax\n";
	myfile << armWeight << ',' << ymin << ',' << ymax << ',' << xmin << ',' << xmax << '\n';
	myfile.close();

	return 0;
}