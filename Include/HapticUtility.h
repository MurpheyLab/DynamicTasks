//---------------------------------------------------------------------
// HapticUtility.h
// Contains utility functions for initializing the HapticMaster and
// formatting returned data strings
//---------------------------------------------------------------------
#ifndef _HAPTIC_UTILITY_H_
#define _HAPTIC_UTILITY_H_

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <windows.h>
//#include <glut.h>
#include "HapticAPI.h"

const int iNrSegments = 10;
const double AxisLength = 0.25; //meter

double LeftSideArc = -0.5; //radians
double RightSideArc = 0.5; //radians

double MinRadius = 0.28; //meter
double MaxRadius = 0.64; //meter

double CenterOffset = -0.45; //meter

double TotalHeight = 0.4; //meter
double PosHeight = 0.2; //meter
double NegHeight = -0.2; //meter

#ifndef Pi
#define Pi 3.14159265358979
#endif

#ifndef Pi2
#define Pi2 6.28318530717958
#endif

#ifndef DegPerRad
#define DegPerRad 180/Pi
#endif

#ifndef RadPerDeg
#define RadPerDeg Pi/180
#endif





//---------------------------------------------------------------------
//               I N I T I A L I Z E   D E V I C E
//
// This function initializes the HapticMASTER device.
// It first searches for the end positions. When all ends are found,
// the HapticMASTER is set to the force_sensitive state.
//---------------------------------------------------------------------
void InitializeDevice(long dev) {

	char outputString[200] = "";
	// set the inertia to the default value
	if (haDeviceSendString(dev, "set inertia 4.0", outputString)) {
  		printf("--- ERROR: set inertia: %s\n",outputString);
	}

	char isCalibratedStr[10] = "";
	if (haDeviceSendString(dev, "get position_calibrated", outputString)) {
		printf("--- ERROR: get position_calibrated: %s\n", outputString);
	} else {
		strcpy_s(isCalibratedStr, outputString);
	}
	printf("%s\n", isCalibratedStr);

	//strcmp(isCalibratedStr, "false;")
	if (isCalibratedStr=="false;") {
   		if (haDeviceSendString(dev, "set state init", outputString)) {
    		printf("--- ERROR: set state init: %s\n", outputString);
   		} else {
    		printf( "Initializing the HapticMASTER. Please wait...\n" );
		}

		if (haDeviceSendString(dev, "get state", outputString)) {
    		printf("--- ERROR: get state: %s\n", outputString);
		}
		printf("%s\n", outputString);

		while (strcmp(outputString, "stop")!=1) {
			if (haDeviceSendString (dev, "get state", outputString)) {
      			printf("--- ERROR: get state: %s\n", outputString);
			}
		}
	}

	printf("Setting to state Force\n");
	if (haDeviceSendString(dev, "set state force", outputString)) {
  		printf("--- ERROR: set state force: %s\n", outputString);
	}
	printf("%s\n", outputString);
}

//---------------------------------------------------------------------
//                 P A R S E   F L O A T   V E C
//
// Break a string vector into its components.
// Returns 0 if was succesful.
//---------------------------------------------------------------------
int ParseFloatVec(const char* inputString, double& xValue, double& yValue, double& zValue ) {
  char xValueStr[100] = "";
	char yValueStr[100] = "";
	char zValueStr[100] = "";

	char c;
	unsigned int srcIndex = 1;

	unsigned int dstIndex = 0;
	c = inputString[srcIndex];
	while (c != ',' && srcIndex < strlen(inputString)) {
		xValueStr[dstIndex] = inputString[srcIndex];
		srcIndex++;
		dstIndex++;
		c = inputString[srcIndex];
	}
	xValueStr[dstIndex] = '\0';
	xValue = atof(xValueStr);

	srcIndex++;
	dstIndex = 0;

	c = inputString[srcIndex];
	while (c != ',' && srcIndex < strlen(inputString)) {
		yValueStr[dstIndex] = inputString[srcIndex];
		srcIndex++;
		dstIndex++;
		c = inputString[srcIndex];
	}
	yValueStr[dstIndex] = '\0';
	yValue = atof(yValueStr);

	srcIndex++;
	dstIndex = 0;

	c = inputString[srcIndex];
	while (c != ']' && srcIndex < strlen(inputString)) {
		zValueStr[dstIndex] = inputString[srcIndex];
		srcIndex++;
		dstIndex++;
		c = inputString[srcIndex];
	}
	zValueStr[dstIndex] = '\0';
	zValue = atof(zValueStr);

	return 0;
}

//---------------------------------------------------------------------
//                   B R E A K   R E S P O N S E
//---------------------------------------------------------------------
int BreakResponse(char* outputString, const char* inputString, unsigned int wantedStringNumber) {
   unsigned int beginIdx = 0;
   unsigned int endIdx = 0;
   unsigned int currIdx = 0;
   char c;
   bool seperatorFound = false;

   for (unsigned int i = 0; i < wantedStringNumber; i++) {
      seperatorFound = false;
      beginIdx = currIdx;
      while (!seperatorFound) {
         c = inputString[currIdx];

         if ( c == ';' || c == '\0' ) {
            seperatorFound = true;
         } else {
            currIdx++;
         }
      } // while
      endIdx = currIdx;
      currIdx++;
   } // for

   for (unsigned int i = beginIdx; i < endIdx; i++) {
      outputString[i-beginIdx] = inputString[i];
   }
   outputString[endIdx-beginIdx] = '\0';

   return 0;
}

//---------------------------------------------------------------------
//          H A   S E N D   C O M M A N D   -   n o   p a r a m s
//
// Overloaded function to send a command to the HapticMASTER.
// This version excpect no extra parameters for the command.
//---------------------------------------------------------------------
int  haSendCommand(long inDev,
	const char* inCommand,
	char* outCommand) {
	return (haDeviceSendString(inDev, inCommand, outCommand));
}

//---------------------------------------------------------------------
//          H A   S E N D   C O M M A N D   -   1   d o u b l e
//
// Overloaded function to send a command to the HapticMASTER.
// This version excpect one double as a parameter for the command.
//---------------------------------------------------------------------
int  haSendCommand(long inDev,
	const char* inCommand,
	double inDouble1,
	char* outCommand) {
	char tempString[100];
	sprintf_s(tempString, "%s %g", inCommand, inDouble1);
	return (haDeviceSendString(inDev, tempString, outCommand));
}

//---------------------------------------------------------------------
//         H A   S E N D   C O M M A N D   -   3   d o u b l e s
//
// Overloaded function to send a command to the HapticMASTER.
// This version excpect 3 doubles as parameters for the command.
//---------------------------------------------------------------------
int  haSendCommand(long inDev,
	const char* inCommand,
	double inDouble1,
	double inDouble2,
	double inDouble3,
	char* outCommand) {
	char tempString[100];
	sprintf_s(tempString, "%s [%g,%g,%g]", inCommand, inDouble1, inDouble2, inDouble3);
	return (haDeviceSendString(inDev, tempString, outCommand));
}

//---------------------------------------------------------------------
//         H A   S E N D   C O M M A N D   -   4   d o u b l e s
//
// Overloaded function to send a command to the HapticMASTER.
// This version excpect 3 doubles as parameters for the command.
//---------------------------------------------------------------------
int  haSendCommand(long inDev,
	const char* inCommand,
	double inDouble1,
	double inDouble2,
	double inDouble3,
	double inDouble4,
	char* outCommand) {
	char tempString[100];
	sprintf_s(tempString, "%s [%g,%g,%g,%g]", inCommand, inDouble1, inDouble2, inDouble3, inDouble4);
	return (haDeviceSendString(inDev, tempString, outCommand));
}


//---------------------------------------------------------------------
//                  S E T  B I A S  F O R C E
//---------------------------------------------------------------------
/*int BiasForce(float support) {
	char* msg;
	char* force;
	char* end_msg;

	msg = "set myBiasForce force [0.0,0.0,";
	force = gcvt(support, 4, force);
	end_msg = "]";

	char* command = (char*) malloc(1 + strlen(msg) + strlen(force) + strlen(end_msg));
    strcpy(command, msg);
    strcat(command, force);
	strcat(command, end_msg);
    printf("%s", command);

	return 0;
}*/

#endif
