//---------------------------------------------------------------------
// stretch.cpp
// Used to allow free motion on the table for stretching between trials
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

#define IPADDRESS "10.30.203.36" //"10.30.203.26"

#define PosX 0
#define PosY 1
#define PosZ 2

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
	double springPos[3] = { -0.04,0.05,-0.17 }; // x was 0.0
	double springStiffness = 4000.0, springDamping = 10.0, springMaxForce = 100.0;

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

	// Fix arm in place using springs 
	printf("Moving...\n");
	returnValue = haSendString(dev, "create spring mySpring", response);
	haSendCommand(dev, "set mySpring pos", springPos[PosX], springPos[PosY], springPos[PosZ] + 0.011, response);
	haSendCommand(dev, "set mySpring stiffness", springStiffness, response);
	haSendCommand(dev, "set mySpring dampfactor", springDamping, response);
	haSendCommand(dev, "set mySpring maxforce", springMaxForce, response);
	returnValue = haSendString(dev, "set mySpring enable", response);

	printf("Press ENTER once arm in position.\n");
	cin.get();
	printf("Table enabled...\n");
	//Set up Haptic Table top
	haSendCommand(dev, "create block myFloor", response);
	haSendCommand(dev, "set myFloor pos", springPos[PosX], springPos[PosY], springPos[PosZ], response);
	haSendCommand(dev, "set myFloor size", 0.7, 0.7, 0.01, response); //set to be a little larger than workspace
	haSendCommand(dev, "set myFloor stiffness", 20000.0, response);
	haSendCommand(dev, "set myFloor enable", response); //comment out if you want to stretch without the table

	printf("Press ENTER to start stretch.\n");
	cin.get();
	printf("Stretch in progress...\n");
	haSendString(dev, "set mySpring disable", response);

	printf("Press ENTER to end stretch.\n");
	cin.get();
	printf("Moving to stop position...\n");
	returnValue = haSendString(dev, "create spring mySpring", response);
	haSendCommand(dev, "set mySpring pos", springPos[PosX], springPos[PosY], springPos[PosZ] + 0.011, response);
	haSendCommand(dev, "set mySpring stiffness", springStiffness, response);
	haSendCommand(dev, "set mySpring dampfactor", springDamping, response);
	haSendCommand(dev, "set mySpring maxforce", springMaxForce, response);
	returnValue = haSendString(dev, "set mySpring enable", response);

	printf("Press ENTER to lock. \n");
	cin.get();
	haSendCommand(dev, "remove all", response);
	printf("remove all ==> %s\n", response);
	haSendCommand(dev, "set state stop", response);
	printf("set state stop ==> %s\n", response);
	if (haCloseDevice(dev)) {
		printf("--- ERROR: closing device\n");
	}

	return 0;
}