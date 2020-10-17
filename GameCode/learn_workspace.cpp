//---------------------------------------------------------------------
// ballinbowl_all.cpp
// Interactive bowl visualization moved using either the HapticMASTER 
// or a PlayStation 4 Controller.
//---------------------------------------------------------------------
// Game header files
#include "ballinbowl.hpp"
#include "parameters.hpp"
// HapticMaster header files
#include "HapticAPI.h"
#include "HapticUtility.h"
// OpenGL header files
#include "glut.h"
#include "glfw3.h"
// Other header files
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define _USE_MATH_DEFINES // this has to come before including the math header file
#include <math.h> // sine, cosine
#include "cmath"
#include <cstdlib>
#include <iostream> // cin
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
#include "paths.hpp"
#include <windows.h> // Sleep()

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//using namespace std;

//Screen dimension constants
const int SCREEN_WIDTH = 6000;
const int SCREEN_HEIGHT = 3000;

// Define possible support levels
//float support_level[3] = { -0.0, -0.1, -0.3 }; // fraction of max shoulder abduction loading (0-1)
const int max_support = 0; // support level index that represents the most load

// Define file for logging flag locations
string logflags = dirpath + "Flags_S" + to_string(subject_num) + "_SL" + to_string(support_num) + ".csv"; 
ofstream logflagfile;

// Define file for logging arm weight and workspace min/max
string setuppath2 = dirpath + "setup_S" + to_string(subject_num) + ".csv";

// Define file for storing path of participant
string logpath = dirpath + "Workspace_Path_S" + to_string(subject_num) + "_SL" + to_string(support_num) + ".csv";
ofstream logpathfile;

#define IPADDRESS "10.30.203.26" //"10.30.203.26"

// Game variables
#define PosX 0
#define PosY 1
#define PosZ 2
#define table_z -0.14 // -0.17
double z_tolerance = table_z + 0.01; // -0.16 // if lower, person can eat flags while slacking
#define tol_delta 0.008 // should be the same as the ball-in-bowl code
#define size_of_flag 0.008 // this value is only used when not using textures

float xmin = 0.0; float xmax = 0.0; float ymin = 0.0; float ymax = 0.0;
float armWeight = 10.0; float maxForce = 10.0;

// Parameters for getting responses from the robot or joystick
long dev = 0;
char response[100];
double CurrentPosition[3];
double CurrentJoystick[2];
double CurrentForce[3];


// Global game parameters to be updated during play
//float time_current = 0.0;
float trial_flag = 0.0;
float scoring_enabled = 0.0;

// Variables to update time during play
double deltaT = 0.01;
double deltaTdefault = 0.05; // 20Hz
double timer_val = 0.0;
double timer_beg;
double prevtime = 0.0;
float sys_time = 0.0;
float beg_sys_time = 0.0;

// flags
float flag_pos[690][3];
int flags_total = 0;
int goal_status[690];

// store path
float path_pos[8000][3];
int path_total = 0;
int path_count = 0;

//---------------------------------------------------------------------
//                  O P E N G L   M A T E R I A L S
//---------------------------------------------------------------------
// General OpenGL Material Parameters
GLfloat Specular[] = { 1.00, 1.00, 1.00, 1.00 };
GLfloat Emissive[] = { 0.00, 0.00, 0.00, 1.00 };
GLfloat Shininess = { 128.00 };
GLfloat SpecularOff[] = { 0.00, 0.00, 0.00, 0.00 };
GLfloat EmissiveOff[] = { 0.50, 0.50, 0.50, 0.00 };
GLfloat ShininessOff = { 0.00 };

//---------------------------------------------------------------------
//                     D R A W   T A B L E
//
// This Function Is Called To Draw A Plane
//---------------------------------------------------------------------
void DrawTable(void) 
{
	int i, j;
	GLfloat v[4][3];
	GLfloat width = 1.5, length = 0.3;

	for (i = 0; i < 4; ++i) { // level with x-y plane
		v[i][2] = 0;
	}

	v[0][0] = -width;
	v[0][1] = length;
	v[1][0] = width;
	v[1][1] = length;
	v[2][0] = width;
	v[2][1] = -length;
	v[3][0] = -width;
	v[3][1] = -length;

	glColor3f(0.0f, 0.0f, 0.0f);
	glBegin(GL_QUADS);
	glVertex3fv(v[0]);
	glVertex3fv(v[1]);
	glVertex3fv(v[2]);
	glVertex3fv(v[3]);
	glEnd();
}


//---------------------------------------------------------------------
//                       D R A W   P A T H
//
// This Function draws the path of the person as a series of red flags
//---------------------------------------------------------------------
void DrawPath(void)
{
	for (int i = 0; i < path_total; i++)
	{
		glColor3f(1.0f, 0.0f, 0.0f); // red
		glPushMatrix();
		glTranslatef(path_pos[i][PosX], path_pos[i][PosY], 0.0 + size_of_flag / 2.0);
		glutSolidCube(size_of_flag); // change this to make the cubes bigger or smaller
		glPopMatrix();
	}
}

//---------------------------------------------------------------------
//                       D R A W   F L A G S
//
// This Function Is Called To Draw Goal Flags (Cubes)
//---------------------------------------------------------------------
void DrawFlags(void) 
{
	for (int i = 0; i < flags_total+1; i++)
	{
		if (goal_status[i] == 0)
		{
			glColor3f(0.0f, 0.99f, 0.0f); // greenS
		}
		else
		{
			glColor3f(0.0f, 0.0f, 0.0f); // black
		}
		glPushMatrix();
		//glTranslatef(goals[i][PosX], goals[i][PosY], goals[i][PosZ] + size_of_flag / 2.0);
		glTranslatef(flag_pos[i][PosX], flag_pos[i][PosY], flag_pos[i][PosZ] + size_of_flag / 2.0);
		glutSolidCube(size_of_flag); // change this to make the cubes bigger or smaller
		glPopMatrix();
	}

}


//---------------------------------------------------------------------
//                    D E F I N E   T A S K
//
// This Function Is Called To Define the Location of Goal Flags
//---------------------------------------------------------------------
void DefineTask(void)
{
	int i = 0;
	float resolution = 0.02;
	for (float xx = -0.25; xx <=0.2; xx=xx+resolution) {
		for (float yy = -0.29; yy <= 0.29; yy=yy+resolution) {
			flag_pos[i][PosX] = xx;
			flag_pos[i][PosY] = yy;
			flag_pos[i][PosZ] = 0.0;
			i++;
		}
	}
	flags_total = i-1; // max index for flag_pos
	printf("Total flags is %d.\n", flags_total);
}

//---------------------------------------------------------------------
// This Function finds the point of intercection of two lines
//---------------------------------------------------------------------
void lineintersection(float Ax, float Ay, float Bx, float By, float Cx, float Cy, float Dx, float Dy, double intersectionP[2]) {
	

	// find equation for line AB 
	double m1 = ((double)By - (double)Ay) / ((double)Bx - (double)Ax);
	double b1 = -m1 * (double)Ax + (double)Ay;

	// find equation for line CD
	double m2 = ((double)Dy - (double)Cy) / ((double)Dx - (double)Cx);
	double b2 = -m2 * (double)Cx + (double)Cy;

	//printf("lines  %f %f %f %f \n", m1,b1,m2,b2);
	if (m1 == m2) {
		// The lines are parallel. 
		// set equal to 1000
		intersectionP[0] = 1000;
		intersectionP[1] = 1000;
	}
	else {
		intersectionP[0] = (b2 - b1) / (m1 - m2);
		intersectionP[1] = m1 * intersectionP[0] + b1;
	}
}

//---------------------------------------------------------------------
//                     C H E C K   F L A G S
//
// This Function Is Called at exit to check if each flag is incircled by
// the person's path. Is does this by drawing a line between the flag and 
// a point outside the workspace and checking to see if it intersects the 
// path.
//---------------------------------------------------------------------
void CheckFlags(void)
{
	//// Chose flags in a circle around the person's current position
	float direction_x[10];
	float direction_y[10];
	int numdir = 10;

	for (int i = 0; i < numdir; i++) {
		direction_x[i] = 2 * cos(2 * 3.1428 * ((float)i / numdir));
		direction_y[i] = 2 * sin(2 * 3.1428 * ((float)i / numdir));
	}

	// define variables for the point of intersection
	double intersect[2];

	// Loop through possible flag locations
	for (int i = 0; i < flags_total+1; i++) 
	{
		// If the flag has already been found to be inside the player's path,
		// don't check the flag again
		if (goal_status[i] == 0) 
		{
			// initialize save_point variable, which increases by 1 for every
			// direction that crosses the player's path. 
			int save_point = 0; 

			// Loop through the directions defined by the point outside the 
			// workspace
			for (int j = 0; j < numdir; j++) 
			{
				
				// Loop through all the path line-segments
				for (int k = 1; k < path_total; k++) 
				{

					lineintersection(flag_pos[i][PosX], flag_pos[i][PosY],
						flag_pos[i][PosX]+direction_x[j], flag_pos[i][PosY]+direction_y[j],
						path_pos[k][PosX], path_pos[k][PosY],
						path_pos[k - 1][PosX], path_pos[k - 1][PosY],
						intersect);

					// check if the intersection point is within each line segment:
					// (path_pos[k],path_pos[k-1]) and (flag_pos[i],direction)
					// true for each if the length of the line segment is greater than distance
					// of the intersection point to each line segment point
					double seg_path_len = sqrt(pow(path_pos[k][PosX] - path_pos[k - 1][PosX], 2) +
						pow(path_pos[k][PosY] - path_pos[k - 1][PosY], 2));
					double seg_path_p1 = sqrt(pow(path_pos[k][PosX] - intersect[0], 2) +
						pow(path_pos[k][PosY] - intersect[1], 2));
					double seg_path_p2 = sqrt(pow(intersect[0] - path_pos[k - 1][PosX], 2) +
						pow(intersect[1] - path_pos[k - 1][PosY], 2));
					double seg_flag_len = sqrt(pow(flag_pos[i][PosX] - direction_x[j], 2) +
						pow(flag_pos[i][PosY] - direction_y[j], 2));
					double seg_flag_p1 = sqrt(pow(flag_pos[i][PosX] - intersect[0], 2) +
						pow(flag_pos[i][PosY] - intersect[1], 2));
					double seg_flag_p2 = sqrt(pow(intersect[0] - direction_x[j], 2) +
						pow(intersect[1] - direction_y[j], 2));


					float mult = 1.05; // to account for numerical error
					float mult2 = 1.001; // to account for numerical error
					if ((seg_path_p1 + seg_path_p2 < mult * seg_path_len) &&
						(seg_flag_p1 + seg_flag_p2 < mult2 * seg_flag_len))
					{
						save_point++;
						break; // don't find more than one crossing in every direction
					}
				}
			}


			// if the line crosses the path almost every direction, save flag to file.			
			if (save_point > numdir - 1) {
				goal_status[i] = 1; // set achieved goal to disappear (turn transparent)
				logflagfile << flag_pos[i][PosX] << "," << flag_pos[i][PosY] << "," << flag_pos[i][PosZ] << "\n";

				// Update xmin, ymin, ymax if needed
				if (flag_pos[i][PosX] < xmin) {
					xmin = flag_pos[i][PosX];
				}
				if (flag_pos[i][PosY] < ymin) {
					ymin = flag_pos[i][PosY];
				}
				if (flag_pos[i][PosY] > ymax) {
					ymax = flag_pos[i][PosY];
				}
			}

		}
	}
}

//---------------------------------------------------------------------
//                    I N I T   O P E N G L
//
// This Function Initializes the OpenGl Graphics Engine
//---------------------------------------------------------------------
void InitOpenGl(void)
{
	glShadeModel(GL_SMOOTH);
	glLoadIdentity();
	GLfloat GrayLight[] = { 0.75, 0.75, 0.75, 1.0 };
	GLfloat LightPosition[] = { 1.0, 2.0, 1.0, 0.0 };
	GLfloat LightDirection[] = { 0.0, 0.0, -1.0, 0.0 };
	glLightfv(GL_LIGHT0, GL_POSITION, LightPosition);
	glLightfv(GL_LIGHT0, GL_AMBIENT, GrayLight);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, GrayLight);
	glLightfv(GL_LIGHT0, GL_SPECULAR, GrayLight);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
	glEnable(GL_COLOR_MATERIAL);
	glClearColor(1.0, 1.0, 1.0, 1.0); // white background
}

//---------------------------------------------------------------------
//                       D R A W   B O W L   B O T T O M
//
// This Function Is Called To Draw Labels for the time and score
//---------------------------------------------------------------------
void DrawBowlBottom(void)
{
	// here we check for: 1. trial is running, 2. person is above the haptic table
	if ((trial_flag == 1.0) && (CurrentPosition[PosZ] > z_tolerance))
	{
		scoring_enabled = 1.0;
		glColor3f(0.01f, 0.5f, 1.0f);
	}
	else
	{
		scoring_enabled = 0.0;
		glColor3f(0.0f, 0.0f, 0.0f);
	}
	glPushMatrix();
	glTranslatef(CurrentPosition[PosX], CurrentPosition[PosY], -0.005);
	glutSolidCube(size_of_flag+1.5*tol_delta); // change this to make the cubes bigger or smaller
	glPopMatrix();
}

//---------------------------------------------------------------------
//                         D I S P L A Y
//
// This Function Is Called By OpenGL To Redraw The Scene
// Here's Where You Put The End Effector And Block Drawing FuntionCalls
//---------------------------------------------------------------------
void Display(void)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glPushMatrix();
	//gluLookAt(0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0); // for top down view 
	//gluLookAt(0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0); // for participant view point
	//gluLookAt(0.5, 0.0, 0.5, -0.2, 0.0, 0.0, 0.0, 0.0, 1.0); // for participant view point
	//gluLookAt(0.645, 0.0, 0.5, -0.1, 0.0, 0.0, 0.0, 0.0, 1.0); // for participant view point
	
	glutPostRedisplay();
	gluLookAt(0.68, 0.0, 0.4, -0.015, 0.0, 0.0, 0.0, 0.0, 1.0);
	DrawTable();
	DrawBowlBottom();
	DrawFlags();
	DrawPath();
	glPopMatrix();
	glutSwapBuffers();
}

//---------------------------------------------------------------------
//                          R E S H A P E
//
// The Function Is Called By OpenGL Whenever The Window Is Resized
//---------------------------------------------------------------------
void Reshape(int iWidth, int iHeight)
{
	glViewport(0, 0, (GLsizei)iWidth, (GLsizei)iHeight);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	float fAspect = (float)iWidth / iHeight;
	gluPerspective(30.0, fAspect, 0.05, 20.0);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
}

//---------------------------------------------------------------------
//                         K E Y B O A R D
//
// This Function Is Called By OpenGL WhenEver A Key Was Hit
//---------------------------------------------------------------------
void Keyboard(unsigned char ucKey, int iX, int iY)
{
	//float max_force = 0.0;
	switch (ucKey)
	{
	case 's': // penalty_timer trial
		// Allow movement 
		printf("S pressed: Trial starting...\n");
		//time_current = 0.0;
		trial_flag = 1.0;
		if (mode == 0) // if using HapticMASTER
		{
			haDeviceSendString(dev, "set myBiasForce enable", response);
			haDeviceSendString(dev, "set mySpring disable", response);
			
		}
		beg_sys_time = clock() / (float)CLOCKS_PER_SEC;
		break;
	case 'c': // remove flags from inside the workspace
		printf("c pressed - checking flags \n");
		CheckFlags();
		break;
	case 'b': // enable bias force
		printf("B pressed \n");
		haDeviceSendString(dev, "set myBiasForce enable", response);
		printf("%s\n", response);
		break;
	case 'd': // disable bias force
		printf("D pressed \n");
		haDeviceSendString(dev, "set myBiasForce disable", response);
		printf("%s\n", response);
		break;
	case 27: // ESC key pressed - stop motion and simulation
		printf("ESC pressed: exit routine starting \n");
		
		// Check to see if any flags are in the enclosed space again, then close file
		CheckFlags();
		logflagfile.close();

		// If it is the heaviest loading level, save arm weight and min/max values
		if (support_num == max_support)
		{
			ofstream myfile;
			myfile.open(setuppath2);
			myfile << armWeight << ',' << ymin << ',' << ymax << ',' << xmin << ',' << xmax << ',' << maxForce << '\n';
			myfile.close();
		}
		
		if (mode == 0) // if using HapticMASTER
		{
			haSendCommand(dev, "remove all", response);
			printf("remove all ==> %s\n", response);
			haSendCommand(dev, "set state stop", response);
			printf("set state stop ==> %s\n", response);
			if (haDeviceClose(dev)) 
			{
				printf("--- ERROR: closing device\n");
			}
		}
		exit(0);
		break;
	}
}

//---------------------------------------------------------------------
//                        T I M E R   C B
//
// This Is A Timer Function Which Calls The HapticMASTER To Get The
// New EndEffector Position (or obtains joystick input) and Updates 
// The Bowl and Ball Position
//---------------------------------------------------------------------
void TimerCB(int iTimer)
{	
	// Get current position
	if (mode == 0) // if using HapticMASTER
	{
		haSendCommand(dev, "get measpos", response);
		ParseFloatVec(response, CurrentPosition[PosX], CurrentPosition[PosY], CurrentPosition[PosZ]);
	}
	else if (mode ==1) // if using a joystick
	{
		int axesCount;
		const float* axes = glfwGetJoystickAxes(GLFW_JOYSTICK_1, &axesCount);
		CurrentJoystick[PosX] = axes[1];
		CurrentJoystick[PosY] = axes[0];
		CurrentPosition[PosX] = CurrentPosition[PosX] + CurrentJoystick[PosX] / 300;
		CurrentPosition[PosY] = CurrentPosition[PosY] + CurrentJoystick[PosY] / 300;
		CurrentPosition[PosZ] = 1; // always above table
	}
	else if (mode == 2) //simulated trajectory 
	{
		CurrentPosition[PosX] = 0.11 * cos(sys_time) + .08 + .04 * sin(10 * sys_time);
		CurrentPosition[PosY] = 0.17 * sin(sys_time) + .05 + .04 * cos(10 * sys_time);
		CurrentPosition[PosZ] = 1;
	}

	// Update simulation once trial has started
	//if (trial_flag == 1.0) 
	//{
	//	time_current = time_current + deltaT;
	//}

	if (trial_flag == 1.0) {
		prevtime = sys_time;
		sys_time = clock() / (float)CLOCKS_PER_SEC - beg_sys_time;
	}
	deltaT = sys_time - prevtime;
	if (deltaT == 0.0) {
		deltaT = deltaTdefault;
	}

	// Store every 10th flag into path_pos- only draws if scoring is enabled
	if (scoring_enabled == 1) {
		if (path_count == 10) {
			path_pos[path_total][PosX] = CurrentPosition[PosX];
			path_pos[path_total][PosY] = CurrentPosition[PosY];
			path_pos[path_total][PosZ] = CurrentPosition[PosZ];
			logpathfile << sys_time << "," << CurrentPosition[PosX] << "," << CurrentPosition[PosY] << "," << CurrentPosition[PosZ] << "\n";
			path_count = 0;
			path_total++;
		}
		path_count++;
	}

    // Poll for and process events
    glfwPollEvents( );

	// Set The Timer For This Function Again
	glutTimerFunc(10, TimerCB, 1);
}

//---------------------------------------------------------------------
//                            M A I N
//---------------------------------------------------------------------
int main(int argc, char** argv)
{
	int returnValue; // for accepting returns from the robot
	
	// Set up the starting position for the spring
	double springPos[3] = { -0.0,0.05, table_z }; //-.17
	double springStiffness = 4000.0, springDamping = 10.0, springMaxForce = 100.0; // this spring is for the home position
	double superSpringStiffness = 4000.0, superSpringDamping = 10.0, superSpringMaxForce = 100.0; // this spring is for measuring maxForce

	// Open file to log flag data, path data, and define an array of flags
	DefineTask();
	logflagfile.open(logflags);
	logpathfile.open(logpath);
	logpathfile << "system_time, x, y, z \n";


	if (mode == 0) // if using HapticMASTER
	{

		// Open the HapticMaster device
		dev = haDeviceOpen(IPADDRESS);

		// Check if device connected, then initialize
		if (dev == HARET_ERROR) {
			printf("--- ERROR: Unable to connect to device: %s\n", IPADDRESS);
			return HARET_ERROR; // comment this out to test visualization without the HapticMaster
		}
		else {
			InitializeDevice(dev);
		}

		haSendCommand(dev, "set inertia 4.0", response);

		// Fix arm in place (home position) using springs 
		//printf("Press ENTER to move to home position.\n");
		//std::cin.get();
		printf("Moving...\n");
		returnValue = haDeviceSendString(dev, "create spring mySpring", response);
		haSendCommand(dev, "set mySpring pos", springPos[PosX], springPos[PosY], springPos[PosZ] + 0.011, response);
		haSendCommand(dev, "set mySpring stiffness", springStiffness, response);
		haSendCommand(dev, "set mySpring dampfactor", springDamping, response);
		haSendCommand(dev, "set mySpring maxforce", springMaxForce, response);
		returnValue = haDeviceSendString(dev, "set mySpring enable", response);

		if (support_num == max_support)
		{

			// Switch to super spring for measuring max force
			printf("Press ENTER to lower position for measuring max force.\n");
			std::cin.get();
			returnValue = haDeviceSendString(dev, "create spring mySuperSpring", response);
			haSendCommand(dev, "set mySuperSpring pos", springPos[PosX], springPos[PosY], -0.20778, response);
			haSendCommand(dev, "set mySuperSpring stiffness", superSpringStiffness, response);
			haSendCommand(dev, "set mySuperSpring dampfactor", superSpringDamping, response);
			haSendCommand(dev, "set mySuperSpring maxforce", superSpringMaxForce, response);
			printf("Moving down...\n");
			returnValue = haDeviceSendString(dev, "set mySpring disable", response);
			returnValue = haDeviceSendString(dev, "set mySuperSpring enable", response);

			// Turn on large bias force so that the participant cannot move the robot up
			printf("Press ENTER to turn on large bias force.\n");
			std::cin.get();
			returnValue = haDeviceSendString(dev, "create biasforce myBiasForce0", response);
			returnValue = haSendCommand(dev, "set myBiasForce0 force", 0.0, 0.0, -200.0, response); // tested using -100.0
			haDeviceSendString(dev, "set myBiasForce0 enable", response);

			// Take 3 measurements for max force
			printf("Press 'f' then ENTER to take max force measurement for 5s.\n");
			printf("Press 's' then ENTER to skip max force measurement.\n");
			char key;
			float max_force;
			for (int j = 0; j < 4; j++) {
				std::cin >> key;
				std::cin.get();
				if (key == 'f') {
					max_force = 0;
					for (int i = 0; i < 500; i++) // take measurements every 10ms for 5 seconds
					{
						haSendCommand(dev, "get measforce", response);
						ParseFloatVec(response, CurrentForce[PosX], CurrentForce[PosY], CurrentForce[PosZ]);
						if (CurrentForce[PosZ] > max_force) {
							max_force = CurrentForce[PosZ];
						}
						Sleep(10);
					}
					printf("Max force is %f Newtons.\n", max_force);
				}
				else {
					printf("Skipping force measurement\n");
					break;
				}
			}
			printf("Please calculate average isometric force measurement.\n");
			printf("Enter the average isometric force measurement, then press ENTER.\n");
			std::cin >> maxForce;
			std::cin.get();
			printf("Setting the max force to %f \n", maxForce);

			// Return to home position
			printf("Press ENTER to return to home position.\n");
			std::cin.get();
			haDeviceSendString(dev, "set myBiasForce0 disable", response);
			printf("Moving up...\n");
			returnValue = haDeviceSendString(dev, "set mySuperSpring disable", response);
			returnValue = haDeviceSendString(dev, "set mySpring enable", response);


			// Weigh arm
			printf("Press 'f' then ENTER to weigh arm.\n");
			printf("Press 's' then ENTER to skip max force measurement.\n");
			//char key;
			float meas_weight;
			for (int j = 0; j < 4; j++) {
				std::cin >> key;
				std::cin.get();
				if (key == 'f') {
					haSendCommand(dev, "get measforce", response);
					ParseFloatVec(response, CurrentForce[PosX], CurrentForce[PosY], CurrentForce[PosZ]);
					meas_weight = -CurrentForce[PosZ];
					printf("Arm Weight is %f Newtons.\n", meas_weight);
				}
				else {
					printf("Skipping arm weight measurement.\n");
					break;
				}
			}
			printf("Please calculate average arm weight.\n");
			printf("Enter the average arm weight measurement, then press ENTER.\n");
			std::cin >> armWeight;
			std::cin.get();
			printf("Setting the arm weight to %f \n", armWeight);

			// Release arm
			printf("Adjust person in seat and press ENTER to allow movement.\n ");
			cin.get();
			returnValue = haDeviceSendString(dev, "set mySpring disable", response);

			// Measure xmax (body sheild position)
			printf("Move closest to the BODY, then press ENTER.\n");
			std::cin.get();
			printf("Measuring BODY...\n");
			haSendCommand(dev, "get measpos", response);
			ParseFloatVec(response, CurrentPosition[PosX], CurrentPosition[PosY], CurrentPosition[PosZ]);
			xmax = CurrentPosition[PosX];
			printf("xmax is %f.\n", xmax);

			// Move back to home
			printf("Press ENTER to move arm back to starting place\n");
			std::cin.get();
			returnValue = haDeviceSendString(dev, "set mySpring enable", response);
		}
		else {
			// read file to get the participants xmax for shield
			#pragma warning (disable : 4996) // because strcpy, fopen, and fscan are depreciated

			// convert string file name to char array
			char setuppath_char[30];
			strcpy(setuppath_char, setuppath2.c_str());

			// read values from first row of file
			FILE* setupfile = fopen(setuppath_char, "r");
			int val2 = fscanf(setupfile, "%f,%f,%f,%f,%f,%f\n", &armWeight, &ymin, &ymax, &xmin, &xmax, &maxForce);
			printf("Read csv for shield position (xmax = %f) and max force (%f) \n", xmax, maxForce);
		}
		
		printf("Press ENTER once arm is in position.\n");
		std::cin.get();

		// Set bias force (not enable yet)
		printf("Setting bias force to %f of the max abduction force.\n", support_level[support_num]);
		float bias = maxForce * support_level[support_num] + armWeight; // changed to subtract arm weight
		returnValue = haDeviceSendString(dev, "create biasforce myBiasForce", response);
		printf("Bias: %f \n", bias);
		printf("Weight: %f \n", armWeight);
		printf("Max: %f \n", maxForce);
		returnValue = haSendCommand(dev, "set myBiasForce force", 0.0, 0.0, bias, response);
		printf("Bias force: %s\n", response);

		//Set up haptic table top 
		haSendCommand(dev, "create block myFloor", response);
		haSendCommand(dev, "set myFloor pos", springPos[PosX], springPos[PosY], springPos[PosZ], response);
		haSendCommand(dev, "set myFloor size", 0.7, 0.7, 0.01, response); //set to be a little larger than workspace
		haSendCommand(dev, "set myFloor stiffness", 80000.0, response);
		haSendCommand(dev, "set myFloor enable", response);
		printf("Table enabled...\n");

		//Set up haptic block to protect person
		haSendCommand(dev, "create block myShield", response);
		haSendCommand(dev, "set myShield pos", xmax, springPos[PosY], springPos[PosZ], response);
		haSendCommand(dev, "set myShield size", 0.01, 0.7, 0.7, response); //set to be a little larger than workspace
		haSendCommand(dev, "set myShield stiffness", 20000.0, response);
		haSendCommand(dev, "set myShield enable", response);
		printf("Personal shield enabled...\n");

		// Initialize start position 
		haSendCommand(dev, "get measpos", response);
		ParseFloatVec(response, CurrentPosition[PosX], CurrentPosition[PosY], CurrentPosition[PosZ]);
	
	}
	else
	{
		// Set initial conditions
		CurrentPosition[PosX] = 0; CurrentPosition[PosY] = 0; CurrentPosition[PosZ] = 0;
	}

	if (mode == 1)
	{
		// Check for joystick
		if (!glfwJoystickPresent(GLFW_JOYSTICK_1))
		{
			std::cout << "ERROR: No Joystick Connected " << std::endl;
			exit(-1);
		}
	}

	printf("Visualization starting...\n");

	// OpenGL Initialization Calls
	glutInit(&argc, argv);
	glfwInit();
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(SCREEN_WIDTH, SCREEN_HEIGHT);

	// Create The OpenGlWindow and Initialize Things
	glutCreateWindow("Workspace Setup");
	InitOpenGl();
	glutReshapeFunc(Reshape);
	glutDisplayFunc(Display);
	glutKeyboardFunc(Keyboard);
	glutTimerFunc(10, TimerCB, 1);
	glutMainLoop();
	glfwTerminate();

	return 0;
}