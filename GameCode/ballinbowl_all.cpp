//---------------------------------------------------------------------
// ballinbowl_all.cpp
// Interactive bowl visualization moved using either the HapticMASTER 
// or a PlayStation 4 Controller.
//---------------------------------------------------------------------
// Game header files
#include "ballinbowl.hpp"
// #include "parameters.hpp"
// HapticMaster header files
#include "HapticAPI.h"
#include "HapticUtility.h"
// OpenGL header files
#include <glut.h>
#include <glfw3.h>
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

// Visualization mode
#define TEXTURES_ON // uncomment this line to use textures

// Option visualization upgrade (textures)
#ifdef TEXTURES_ON
#include "texture.h"
static Renderer* renderer; //this is a static pointer to a Renderer used in the glut callback functions
#else
#include "paths.hpp"
#endif

//---------------------------------------------------------------------
//---------------------------------------------------------------------
//using namespace std;

//Screen dimension constants
const int SCREEN_WIDTH = 3500;
const int SCREEN_HEIGHT = 2000;

// Define possible support levels
//float support_level[3] = { -0.0, -0.1, -0.3 }; // fraction of max shoulder abduction loading (0-1)

// Define possible tasks & setup variables
const char* tasklist[6] = { taskpath0, taskpath1, taskpath2, taskpath3, taskpath4, taskpath5};

// Define file for logging data - this is assuming game version 'f'
const char* taskpath = tasklist[task_num];
//string logpath = dirpath + "OutputData_S" + to_string(subject_num) + "_Trial" + to_string(trial_num) + "_Freq" + to_string(freq_num) + "_SL0_F" + to_string(feedback_forces) + "_B" + to_string(ball_moving) + ".csv";
string logpath = dirpath + "OutputData_S" + to_string(subject_num) + "_Trial" + to_string(trial_num) + "_Freq" + to_string(freq_num) + "_SL" + to_string(support_num) + "_A" + to_string(arm) + ".csv";

// Define file for logging data - this is assuming game version 'i'
//string logpath = "OutputData_S" + to_string(subject_num) + "_Freq" + to_string(freq) + "_SL" + to_string(support_num) + "_Trial" + to_string(trial_num) + ".csv";

ofstream logfile;

#define IPADDRESS "10.30.203.26" //"10.30.203.26"

#define PosX 0
#define PosY 1
#define PosZ 2

double deltaT = 0.01;
double deltaTdefault = 0.05; // 0.05; // 20Hz
double timer_val = 0.0;
double timer_beg;
double prevtime = 0.0;
#define R_vis 0.07 // radius of bowl for visualization purposes
#define factor 20.0 // 60.0 without feedback forces, 20.0 with feedback forces
#define table_z -0.14 // -0.17
double z_tolerance = table_z + 0.01; // -0.16 // if lower, person can eat flags while slacking

#define tol_delta 0.008
#define size_of_flag 0.008 

int pb_vec[num_freqs_tested][2];

static GLUquadricObj * q; // water object

GLfloat xmin = 10.0; GLfloat xmax = 10.0; GLfloat ymin = 10.0; GLfloat ymax = 10.0;
GLfloat armWeight = 10.0;  GLfloat maxForce = 10.0;

long dev = 0;
char response[100];
double CurrentPosition[3];
double CurrentJoystick[2];
double CurrentForce[3];
double BallPosition[3];
double ballXprev[3] = { 0.0,0.0,0.0 }; double ballYprev[3] = { 0.0,0.0,0.0 };
char ball_direction, ball_direction_temp; // used for water droplets
int ball_intensity, ball_intensity_temp; // used for ball color and water droplets

double springPos[3];

int score = 0;
float sys_time = 0.0;
float beg_sys_time = 0.0;
float trial_flag = 0.0;
float scoring_enabled = 0.0;

// pysical system parameters
double gravity = 9.81;
double mass = 1.0;
BallBowl sys(mass, -damping, -gravity, R, deltaTdefault); // m, B, g, radius, dt in seconds
double xprev[3] = {0.0,0.0,0.0}; double yprev[3] = {0.0,0.0,0.0};
double thetaXprev[3] = {0.0,0.0,0.0}; double thetaYprev[3] = {0.0,0.0,0.0};
double deltaTvec[3] = { deltaTdefault, deltaTdefault, deltaTdefault };
float ballXaccPREV = 0.0; float ballYaccPREV = 0.0;
double ballXaccVEC[3] = { 0.0,0.0,0.0 }; double ballYaccVEC[3] = { 0.0,0.0,0.0 };
int top_flag = 0;

// define ball
#define radius 0.01
#define RR_vis (R_vis - radius) // radius of bowl-ball for visualization purposes
//#define RR R //(R - radius) 

// ball energy
double ball_energy = 0.0;
double max_energy = mass * gravity * percent_height * R;
double mid_energy = mass * gravity * percent_height * R / 2;
//
//float max_energy = mass * gravity * percent_height * RR_vis;
//float mid_energy = mass * gravity * percent_height * RR_vis / 2.0;

// water
clock_t drop_timer;
bool reset_drop_timer = true;
double drop_duration;

// flags
float flag_pos[500][3]; // list of all reachable flags 
int num_reachable = 0;

GLfloat goals[num_csv_flags][3], goals_rand[num_rand_flags][3]; 
int goal_status[num_csv_flags]; 

double get_random() { return (double)rand() / (double)RAND_MAX; } // function for generating a random number [0,1]

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

#ifndef TEXTURES_ON
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
//                       D R A W   F L A G S
//
// This Function Is Called To Draw Goal Flags (Cubes)
//---------------------------------------------------------------------
void DrawFlags(void) 
{
	if (game_version == 'f')
	{
		for (int i = 0; i < num_csv_flags; i++)
		{
			if (goal_status[i] == 0)
			{
				glColor3f(0.0f, 0.99f, 0.0f);
			}
			else
			{
				glColor3f(0.0f, 0.0f, 0.0f);
			}
			glPushMatrix();
			glTranslatef(goals[i][PosX], goals[i][PosY], goals[i][PosZ] + size_of_flag / 2.0);
			glutSolidCube(size_of_flag); // change this to make the cubes bigger or smaller
			glPopMatrix();
		}
	}
	else
	{
		for (int i = 0; i < num_rand_flags; i++)
		{
			glColor3f(0.0f, 0.99f, 0.0f);
			glPushMatrix();
			glTranslatef(goals_rand[i][PosX], goals_rand[i][PosY], goals_rand[i][PosZ] + size_of_flag / 2.0);
			glutSolidCube(size_of_flag); // change this to make the cubes bigger or smaller
			glPopMatrix();
		}
	}
}
#endif

//---------------------------------------------------------------------
//                     D R A W   B A L L
//
// This Function Is Called To Draw The Ball
//---------------------------------------------------------------------
void DrawBall(void)
{
	
	BallPosition[PosX] = CurrentPosition[PosX] + RR_vis * sin(sys.Xcurr[0]) * cos(sys.Xcurr[2]);
	BallPosition[PosY] = CurrentPosition[PosY] + RR_vis * sin(sys.Xcurr[2]) * cos(sys.Xcurr[0]);
	BallPosition[PosZ] = RR_vis - RR_vis * cos(sys.Xcurr[0]) * cos(sys.Xcurr[2]);
	//BallPosition[PosZ] = RR_vis - cos(asin(sqrt(pow(sin(sys.Xcurr[0]),2) + pow(sin(sys.Xcurr[2]),2))))*RR_vis;
	
	//ball_energy = (mass * gravity * BallPosition[PosZ]) + (0.5 * mass * (pow(RR_vis * sys.Xcurr[1], 2) + pow(RR_vis * sys.Xcurr[3], 2)));
	ball_energy = (mass * gravity * BallPosition[PosZ] * R/RR_vis) + (0.5 * mass * (pow(R*sys.Xcurr[1], 2) + pow(R * sys.Xcurr[3],2)));


	if (ball_energy > max_energy * 4 && ball_moving==1) {
		glColor3f(1.0f, 0.0f, 0.0f); // red
		ball_intensity = 2; 
		//printf("ball energy %f \n", ball_energy/ max_energy);
	}
	else if (ball_energy > max_energy * 3.25 && ball_moving == 1) {

		glColor3f(1.0f, 0.25f, 0.0f); // red - orange
		ball_intensity = 1;
	}
	else if (ball_energy > max_energy* 2.5 && ball_moving == 1) {
		glColor3f(1.0f, 0.5f, 0.0f); // orange
		ball_intensity = 1;
	}
	else if (ball_energy > max_energy * 1.75 && ball_moving == 1) {
		glColor3f(1.0f, 0.75f, 0.0f); // yellow - orange
		ball_intensity = 1;
	}
	else if (ball_energy > max_energy || ball_moving == 0) {
		glColor3f(1.0f, 1.0f, 0.0f); // yellow
		ball_intensity = 0;
	}
	else {
		glColor3f(0.0f, 1.0f, 0.0f && ball_moving == 1); // green
		ball_intensity = 0;
	}


	//// New Protocol - based on ball energy
	//if (ball_energy > mid_energy && ball_energy <= max_energy)
	//{
	//	glColor3f(1.0f, 0.5f, 0.0f); // orange
	//	//printf("ball height: %f \n", BallPosition[PosZ]);
	//	ball_intensity = 1; // DON'T CHANGE   ball_intensity is on a range from 0 - 2
	//}
	//else if (ball_energy > max_energy)
	//{
	//	glColor3f(1.0f, 0.0f, 0.0f); // red
	//	ball_intensity = 2; // DON'T CHANGE
	//}
	//else 
	//{
	//	glColor3f(1.0f, 1.0f, 0.0f); // yellow
	//	printf("ball height: %f \n", BallPosition[PosZ]);
	//	ball_intensity = 0; // DON'T CHANGE
	//} 

	glPushMatrix();
	if (ball_moving == 1) {
		glTranslatef(BallPosition[PosX], BallPosition[PosY], BallPosition[PosZ] + radius);
	} else {
		glTranslatef(CurrentPosition[PosX], CurrentPosition[PosY], radius);
	}
	glutSolidSphere(radius, 20, 20);
	glPopMatrix();
}

//---------------------------------------------------------------------
//                     D R A W   B O W L
//
// This Function Is Called To Draw The Bowl
//---------------------------------------------------------------------
void DrawBowl(void) 
{
	int i, j;
	GLfloat r = R_vis;
	// change these to render more or less circles
	const int scalex = 60; 
	const int scaley = 60;
	GLfloat v[scalex * scaley][3];

	for (i = 0; i < scalex; ++i) 
	{
		for (j = 0; j < scaley; ++j) 
		{
			v[i * scaley + j][0] = CurrentPosition[PosX] + r * cos(j * 2 * M_PI / scaley) * cos(i * M_PI / (2 * scalex));
			v[i * scaley + j][1] = CurrentPosition[PosY] + r * sin(j * 2 * M_PI / scaley) * cos(i * M_PI / (2 * scalex));
			v[i * scaley + j][2] = r - r * sin(i * M_PI / (2 * scalex));
		}
	}

	if (scoring_enabled==1.0) 
	{
		glColor3f(0.3f, 0.0f, 0.51f);
	}
	else 
	{
		glColor3f(0.0f, 0.0f, 0.0f);
	}

	glBegin(GL_LINES);
	for (i = 26; i < scalex - 1; ++i) // change the penalty_timer condition for i to make the bowl shorter or taller (was i=1 for full bowl)
	{ 
		for (j = 0; j < scaley; ++j) 
		{
			glVertex3fv(v[i * scaley + j]);
			glVertex3fv(v[i * scaley + (j + 1) % scaley]);
			glVertex3fv(v[(i + 1) * scaley + (j + 1) % scaley]);
			glVertex3fv(v[(i + 1) * scaley + j]);
		}
	}
	glEnd();

}

//---------------------------------------------------------------------
//                     D R A W   W A T E R
//
// This Function Is Called To Draw The Water In The Bowl
//---------------------------------------------------------------------
void DrawWater(void)
{
	const GLint numberOfSides = 360; 
	const GLint numberOfVertices = numberOfSides + 2;	

	GLfloat circleVerticesX[numberOfVertices];
	GLfloat circleVerticesY[numberOfVertices];
	GLfloat circleVerticesZ[numberOfVertices];

	// DON'T CHANGE must all be 0
	circleVerticesX[0] = 0; 
	circleVerticesY[0] = 0; 
	circleVerticesZ[0] = 0;

	//double r = sqrt(pow(R_vis,2) + pow((1.0-percent_height)*R_vis,2));
	double r = 0.05;
	for (int i = 1; i < numberOfVertices; i++)
	{
		circleVerticesX[i] = r * cos(i * 2 * M_PI / numberOfSides);
		circleVerticesY[i] = r * sin(i * 2 * M_PI / numberOfSides);
		circleVerticesZ[i] = circleVerticesZ[0];
	}

	// combine X, Y, Z coordinates to single array
	GLfloat allCircleVertices[numberOfVertices * 3];
	
	for (int i = 0; i < numberOfVertices; i++)
	{
		allCircleVertices[i*3] = circleVerticesX[i];
		allCircleVertices[(i*3) + 1] = circleVerticesY[i];
		allCircleVertices[(i*3) + 2] = circleVerticesZ[i];
	}

	glPushMatrix();
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnable(GL_BLEND);
	glColor4f(0.31f, 0.54f, 0.95f, 0.35f);

	if(!q)
	{
		q = gluNewQuadric();
		gluQuadricDrawStyle(q, GLU_FILL);
	}

	double ball_angle_x = fmod(sys.Xcurr[2]*180/M_PI, 360);
	double ball_angle_y = -1*fmod(sys.Xcurr[0]*180/M_PI, 360);
	double water_angle_limit = 60.0; // ADJUSTABLE max angle that the water can reach w.r.t. the x-y plane

	// set limits on water angle about the x-axis
	if ((ball_angle_x > water_angle_limit && ball_angle_x < 180.0) || (ball_angle_x <= -180.0 && ball_angle_x > -360 + water_angle_limit))
	{
		ball_angle_x = water_angle_limit;
	}
	else if ((ball_angle_x >= 180.0 && ball_angle_x < 360.0 - water_angle_limit) || (ball_angle_x < -water_angle_limit && ball_angle_x > -180.0))
	{
		ball_angle_x = -water_angle_limit;
	}
	else
	{
		ball_angle_x = fmod(sys.Xcurr[2]*180/M_PI, 360);
	}	

	// set limits on water angle about the y-axis
	if ((ball_angle_y > water_angle_limit && ball_angle_y < 180.0) || (ball_angle_y <= -180.0 && ball_angle_y > -360 + water_angle_limit))
	{
		ball_angle_y = water_angle_limit;
	}
	else if ((ball_angle_y >= 180.0 && ball_angle_y < 360.0 - water_angle_limit) || (ball_angle_y < -water_angle_limit && ball_angle_y > -180.0))
	{
		ball_angle_y = -water_angle_limit;
	}
	else
	{
		ball_angle_y = -1*fmod(sys.Xcurr[0]*180/M_PI, 360);
	}

	// determine direction of the ball in the y direction (Left or Right)
	if ((ball_angle_x > -water_angle_limit/3 && ball_angle_x < water_angle_limit/3 && ball_angle_y > -water_angle_limit/3 && ball_angle_y < water_angle_limit/3) || ball_angle_x > 360.0 - water_angle_limit/3 || ball_angle_x < -360.0 + water_angle_limit/3)
	{
		ball_direction = 'n';
	}
	else if ((ball_angle_x > 0.0 && ball_angle_x < 180.0)|| (ball_angle_x < -180.0 && ball_angle_x > -360.0))
	{
		ball_direction = 'r';
	}
	else
	{
		ball_direction = 'l';
	}

	// draw water hemisphere   
	glTranslatef(CurrentPosition[PosX], CurrentPosition[PosY], R_vis);
	glRotatef(ball_angle_y, 0.0, 1.0, 0.0);
	glRotatef(ball_angle_x, 1.0, 0.0, 0.0);
	glEnable(GL_CLIP_PLANE0);
	double clipEq[4] = { 0.0, 0.0, -1.0, -(1 - percent_height) * R_vis }; // ADJUSTABLE plane that removes the top half of the sphere of water
   	glClipPlane(GL_CLIP_PLANE0, clipEq);
	gluSphere(q, R_vis, 16, 16);
	glDisable(GL_CLIP_PLANE0);

	// draw water surface
	glTranslatef(0.0, 0.0, -(1-percent_height)*R_vis ); // z offset from bowl sphere center
	glVertexPointer(3, GL_FLOAT, 0, allCircleVertices);
	glDrawArrays(GL_TRIANGLE_FAN, 0, numberOfVertices);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisable(GL_BLEND);
	glPopMatrix();
}

//---------------------------------------------------------------------
//                    D E F I N E   T A S K
//
// This Function Is Called To Define the Location of Goal Flags
//---------------------------------------------------------------------
void DefineTask(void)
{
	// game version determines if flag positions are read from a csv file or randomly chosen
	if (game_version == 'f')
	{
		int i = 0;
		float x_temp, y_temp, z_temp;
		GLfloat norm[num_csv_flags][3];
		string line, field;


		// read in values from .csv
		FILE* taskfile;
		#pragma warning (disable : 4996)
		taskfile = fopen(taskpath, "r");

		#pragma warning (disable : 4996)
		while (fscanf(taskfile, "%f,%f,%f\n", &x_temp, &y_temp, &z_temp) == 3)
		{
			norm[i][PosX] = x_temp;
			norm[i][PosY] = y_temp;
			norm[i][PosZ] = z_temp;
			i++;
		}
		fclose(taskfile);

		// scale based on reachable workspace
		for (i = 0; i < num_csv_flags; i++)
		{
			goals[i][PosX] = xmax * norm[i][PosX] + xmin * (1 - norm[i][PosX]);
			goals[i][PosY] = ymax * norm[i][PosY] + ymin * (1 - norm[i][PosY]);
			goals[i][PosZ] = norm[i][PosZ];
		}
	}
	else
	{
		int i = 0;
		float x_temp, y_temp, z_temp;

		// read in values from .csv
		FILE* flagsfile;
		#pragma warning (disable : 4996)
		flagsfile = fopen(flagspath, "r");

		#pragma warning (disable : 4996)
		while (fscanf(flagsfile, "%f,%f,%f\n", &x_temp, &y_temp, &z_temp) == 3)
		{
			flag_pos[i][PosX] = x_temp;
			flag_pos[i][PosY] = y_temp;
			flag_pos[i][PosZ] = z_temp;
			i++;
		}
		fclose(flagsfile);
		num_reachable = i-1; // max index for flag_pos

		// select first set of flags
		int rand_int;
		for (i = 0; i < num_rand_flags; i++)
		{
			double temp = get_random()*num_reachable;
			rand_int = (int) (temp + 0.5);
			goals_rand[i][PosX] = flag_pos[rand_int][PosX];
			goals_rand[i][PosY] = flag_pos[rand_int][PosY];
			goals_rand[i][PosZ] = flag_pos[rand_int][PosZ];
		}
    }
}

//---------------------------------------------------------------------
//                     C H E C K   F L A G S
//
// This Function Is Called By TimerCB to check if a flag has been
// reached
//---------------------------------------------------------------------
void CheckFlags(void)
{
    GLfloat rand_x, rand_y;

	// here we check for: 1. ball is in bowl, 2. trial is running, 3. person is above the haptic table
	//if (((ball_energy < max_energy) || (ball_moving == 0)) && (trial_flag == 1.0) && ((CurrentPosition[PosZ] > z_tolerance) || (support_level[support_num]>0.5)))
	if (((ball_energy < max_energy) || (ball_moving == 0)) && (trial_flag == 1.0) && (CurrentPosition[PosZ] > z_tolerance))

	{
		scoring_enabled = 1.0;
	}
	else 
	{
		scoring_enabled = 0.0;
	}

	// 
	// game version determines how flags are scored
	//
	if (game_version == 'f')
	{
		for (int i = 0; i < num_csv_flags; i++) 
		{
			// ADJUSTABLE first two if-statements determine the size and shape of the flag hit box  
			if ((abs(CurrentPosition[PosX] - goals[i][PosX]) < tol_delta) && (abs(CurrentPosition[PosY] - goals[i][PosY]) < tol_delta) && (goal_status[i]==0) && (scoring_enabled==1.0))
			// eats flags when only x and y overlap
			{
				goal_status[i] = 1; // set achieved goal to disappear (turn transparent)
				score++;
			}
		}
	}
	else if (game_version == 'i')
	{
		for (int i = 0; i < num_rand_flags; i++) 
		{
			if ((abs(CurrentPosition[PosX] - goals_rand[i][PosX]) < tol_delta) && (abs(CurrentPosition[PosY] - goals_rand[i][PosY]) < tol_delta) && (scoring_enabled==1.0)) //&& (goal_rand_status[i]==0)
			{
				// replace eaten flag with new, randomly placed flag
				double temp = get_random() * num_reachable;
				int rand_int = (int)(temp + 0.5);

				goals_rand[i][PosX] = flag_pos[rand_int][PosX];
				goals_rand[i][PosY] = flag_pos[rand_int][PosY];
				goals_rand[i][PosZ] = flag_pos[rand_int][PosZ];
				score++;

				float rand_num = get_random() - 0.5;
				float sign = rand_num / abs(rand_num);
				float x_ang = sign * get_random() * angle_disturbance; // change to be both positive and negative
				rand_num = get_random() - 0.5;
				sign = rand_num / abs(rand_num);
				float y_ang = sign * asin(sqrt(pow(sin(angle_disturbance), 2) - pow(sin(x_ang), 2)));
				//float y_ang = sign*sqrt(pow(angle_disturbance, 2) - pow(x_ang, 2));
				//printf("sign: %f x_angle: %f y_angle: %f", sign, x_ang, y_ang);

				sys.Xcurr = { x_ang,0.0,y_ang,0.0,CurrentPosition[PosX],0.0,CurrentPosition[PosY],0.0 }; // M_PI / 3.0

				//sys.Xcurr = { M_PI / 3.0,0.0,M_PI / 3.0,0.0,CurrentPosition[PosX],0.0,CurrentPosition[PosY],0.0 };
			}
		}
	}
	else
	{
		cout << "Not a valid game version. Change game_version to 'f' for reading flag coords from csv file, 'i' for infinite flags, and 'p' for penalty for dropping the ball" << endl;
		exit(0);
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
//                       D R A W   L A B E L S
//
// This Function Is Called To Draw Labels for the time and score
//---------------------------------------------------------------------

void DrawText(GLfloat x, GLfloat y, const char* text)
{
    glPushMatrix();
    glTranslatef(x, y, -0.95);
    glScalef(1/3050.0, 1/3050.0, 1/3050.0);
    for( const char* p = text; *p; p++)
    {
        glutStrokeCharacter(GLUT_STROKE_ROMAN, *p);
    }
    glPopMatrix();
}

void DrawLabels(void) {

	string score_string = to_string(score);
	char* cscore = const_cast<char*>(score_string.c_str());

	int time_current_rounded = round(TF - sys_time);
	string time_current_string = to_string(time_current_rounded);
	char* ctime = const_cast<char*>(time_current_string.c_str());

	//glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// thicken the text
	for (float i = 0.0; i <= 0.002; i = i + 0.0001) {
		for (float j = 0.0; j <= 0.002; j = j + 0.0001) {
			glColor3f(0.0f, 0.7f, 1.0f);
			DrawText(-0.44 + i, 0.19 - j, "SCORE");

			if (score < 10)
			{
				DrawText(-0.385 + i, 0.14 - j, cscore);
			}
			else
			{
				DrawText(-0.395 + i, 0.14 - j, cscore);
			}

			glColor3f(1.0f, 0.0f, 0.0f);
			DrawText(0.325 + i, 0.19 - j, "TIME");

			if (time_current_rounded < 10)
			{
				DrawText(0.355 + i, 0.14 - j, ctime);
			}
			else
			{
				DrawText(0.345 + i, 0.14 - j, ctime);
			}
		}
	}
}

//----------------------------------------------------------------------------------
//               D R A W   P E R S O N A L   B E S T
//
// This Function Is Called To Draw A Score Board With Participant's Personal Bests
//----------------------------------------------------------------------------------
void DrawPersonalBestBoard(void) {

	int i, j;
	GLfloat v[4][3];
	GLfloat bottom = 0.07, width = 0.9, height = 0.04;

	for (i = 0; i < 4; ++i)
	{ // parallel with x-y plane
		v[i][0] = 0.25;
	}

	v[0][2] = -0.1;
	v[0][1] = -width / 2;
	v[1][2] = height;
	v[1][1] = -width / 2;
	v[2][2] = -0.1;
	v[2][1] = width / 2;
	v[3][2] = height;
	v[3][1] = width / 2;

	glColor3f(1.0f, 1.0f, 1.0f); // white
	//glColor3f(0.0f, 0.0f, 0.0f); // black
	glBegin(GL_QUADS);
	glVertex3fv(v[0]);
	glVertex3fv(v[1]);
	glVertex3fv(v[2]);
	glVertex3fv(v[3]);
	glEnd();
}

void DrawText3D(GLfloat x, GLfloat y, GLfloat z, const char* text)
{
	glPushMatrix();
	glTranslatef(x, y, z);
	glScalef(1 / 8050.0, 1 / 8050.0, 1 / 8050.0);
	for (const char* p = text; *p; p++)
	{
		glutStrokeCharacter(GLUT_STROKE_ROMAN, *p);
	}
	glPopMatrix();
}

void DrawPersonalBestLabel(void) {
	string score_string = to_string(pb_vec[freq_num][condition]);
	char* cscore = const_cast<char*>(score_string.c_str());

    glLoadIdentity();

	// thicken the text
	for( float i = 0.0; i <= 0.0012; i = i + 0.0001) {
		for (float j = 0.0; j <= 0.0012; j = j + 0.0001) {
			glColor3f(0.0f, 0.0f, 0.0f);
			DrawText3D(-0.085 + i, -0.08 - j, -0.36, "PERSONAL BEST:");

			DrawText3D(0.07 + i, -0.08 - j, -0.36, cscore);

		}
    }
}


//---------------------------------------------------------------------
//                     D R A W   T I M E R   B A R
//
// This Function Is Called To Draw A Timer Bar
//---------------------------------------------------------------------
void DrawTimerBar(void) 
{
	int i, j;
	GLfloat v[4][3];
	GLfloat bottom = 0.07, width = 0.053, x_offset = 0.3, length = 0.89 * (1 - sys_time / TF);

	for (i = 0; i < 4; ++i) 
	{ // parallel with x-y plane
		v[i][2] = 0.0;
	}

	v[0][0] = bottom;
	v[0][1] = x_offset + width;
	v[1][0] = bottom - length;
	v[1][1] = x_offset + width;
	v[2][0] = bottom - length;
	v[2][1] = x_offset;
	v[3][0] = bottom;
	v[3][1] = x_offset;

	glColor3f(1.0f, 0.0f, 0.0f);
	glBegin(GL_QUADS);
	glVertex3fv(v[0]);
	glVertex3fv(v[1]);
	glVertex3fv(v[2]);
	glVertex3fv(v[3]);
	glEnd();
}

//---------------------------------------------------------------------
//                     D R A W   S C O R E   B A R
//
// This Function Is Called To Draw A Score Bar
//---------------------------------------------------------------------
void DrawScoreBar(void) 
{
	int i, j;
	GLfloat v[4][3];
	GLfloat bottom = 0.07, width = 0.053, x_offset = 0.3, length = 0.89 * score / num_rand_flags;

	for (i = 0; i < 4; ++i) { // parallel with x-y plane
		v[i][2] = 0.0;
	}

	v[0][0] = bottom;
	v[0][1] = -x_offset - width;
	v[1][0] = bottom - length;
	v[1][1] = -x_offset - width;
	v[2][0] = bottom - length;
	v[2][1] = -x_offset;
	v[3][0] = bottom;
	v[3][1] = -x_offset;

	glColor3f(0.0f, 0.7f, 1.0f);
	glBegin(GL_QUADS);
	glVertex3fv(v[0]);
	glVertex3fv(v[1]);
	glVertex3fv(v[2]);
	glVertex3fv(v[3]);
	glEnd();
}

//---------------------------------------------------------------------
//                       D R A W   B O W L   B O T T O M
//
// This Function Is Called To Draw Labels for the time and score
//---------------------------------------------------------------------
void DrawBowlBottom(void)
{
	if (scoring_enabled==1.0) 
	{
		glColor3f(0.01f, 0.5f, 1.0f);
	}
	else 
	{
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
	
	glutPostRedisplay();

	#ifdef TEXTURES_ON
	// the order of these elements MATTERS
	renderer->DrawTable();
	DrawBall();
	DrawBowl();
	DrawTimerBar();
	if (trial_flag == 0.0 && personalbests_flag==1) { DrawPersonalBestBoard(); }
	DrawBowlBottom();
	if (game_version == 'f')
	{
		renderer->DrawFlags(goals, goal_status, size_of_flag);
	}
	else
	{
		renderer->DrawRandomFlags(goals_rand, size_of_flag);
	}
	if (ball_moving == 1) {	DrawWater();}
	DrawLabels();
	if (trial_flag == 0.0 && personalbests_flag == 1) { DrawPersonalBestLabel(); }
	// water drop flag is true if timer is reset
	if (reset_drop_timer)
	{
		drop_timer = clock();
		reset_drop_timer = false;
		ball_direction_temp = ball_direction;
		ball_intensity_temp = ball_intensity;
	}
	else
	{
		// reset timer after a set time has passed
		drop_duration = (clock() - drop_timer) / (double)CLOCKS_PER_SEC;
		if (drop_duration > 0.3)
		{
			reset_drop_timer = true;
		}
	}
	if (ball_moving == 1) {renderer->DrawWaterDrops(CurrentPosition, ball_direction_temp, ball_intensity_temp, drop_duration, RR_vis); }
	#else
	//gluLookAt(0.0, 0.0, 3.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0); // for top down view 
	//gluLookAt(0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0); // for participant view point
	//gluLookAt(0.5, 0.0, 0.5, -0.2, 0.0, 0.0, 0.0, 0.0, 1.0); // for participant view point
	//gluLookAt(0.645, 0.0, 0.5, -0.1, 0.0, 0.0, 0.0, 0.0, 1.0); // for participant view point
	gluLookAt(0.68, 0.0, 0.4, -0.015, 0.0, 0.0, 0.0, 0.0, 1.0);
	DrawTable();
	DrawBall();
	DrawBowl();
	DrawBowlBottom();
	DrawTimerBar();
	if (trial_flag == 0.0 && personalbests_flag == 1) { DrawPersonalBestBoard(); }
	DrawFlags();
	if (ball_moving == 1) { DrawWater(); }
	DrawLabels();
	if (trial_flag == 0.0 && personalbests_flag == 1) { DrawPersonalBestLabel(); }
	#endif 


	// no score bar for infinite game version
	if (game_version == 'f' || game_version == 'p')
	{
		DrawScoreBar();
	}

	glPopMatrix();
	glutSwapBuffers();

	// draw image when participant not lifting
	if (trial_flag == 1.0 && CurrentPosition[PosZ] < z_tolerance) {
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		glPushMatrix();
		glutPostRedisplay();
		renderer->DrawLift();
		//DrawPersonalBestBoard;
		//DrawPersonalBestLabel;
		glPopMatrix();
		glutSwapBuffers();
	}
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
	switch (ucKey)
	{
	case 's': // penalty_timer trial
		// Allow movement and penalty_timer timer function
		printf("S pressed: Trial and timer starting...\n");
		
		if (mode == 0) // if using HapticMASTER
		{
			haDeviceSendString(dev, "set myBiasForce enable", response);
			haDeviceSendString(dev, "set mySpring disable", response);
		}

		while (CurrentPosition[PosZ] < z_tolerance) {
			haSendCommand(dev, "get measpos", response);
			ParseFloatVec(response, CurrentPosition[PosX], CurrentPosition[PosY], CurrentPosition[PosZ]);
		}
		trial_flag = 1.0;
		beg_sys_time = clock() / (float)CLOCKS_PER_SEC;
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
		logfile.close();

		// Start normal shutdown routine
		if (mode == 0) // if using HapticMASTER
		{
			// Move the person back to their home position by enabling springs
			printf("Moving to home position\n");
			haDeviceSendString(dev, "set mySpring enable", response);

			// Disable springs when ready 
			double pos_tol = 0.01;
			while ((abs(CurrentPosition[PosX] - springPos[PosX]) > pos_tol) || (abs(CurrentPosition[PosY] - springPos[PosY]) > pos_tol) || (abs(CurrentPosition[PosZ] - springPos[PosZ] - 0.011) > pos_tol)) {
				haSendCommand(dev, "get measpos", response);
				ParseFloatVec(response, CurrentPosition[PosX], CurrentPosition[PosY], CurrentPosition[PosZ]);
			}
			haDeviceSendString(dev, "set mySpring disable", response);

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
	#ifdef TEXTURES_ON
	float offset = 1.0f;
	renderer->t += offset;
	#endif

	// Get current position
	if (mode == 0) // if using HapticMASTER
	{
		haSendCommand(dev, "get measpos", response);
		ParseFloatVec(response, CurrentPosition[PosX], CurrentPosition[PosY], CurrentPosition[PosZ]);
	}
	else // if using a joystick
	{
		int axesCount;
		const float* axes = glfwGetJoystickAxes(GLFW_JOYSTICK_1, &axesCount);
		CurrentJoystick[PosX] = axes[1];
		CurrentJoystick[PosY] = axes[0];
		CurrentPosition[PosX] = CurrentPosition[PosX] + CurrentJoystick[PosX] / 300;
		CurrentPosition[PosY] = CurrentPosition[PosY] + CurrentJoystick[PosY] / 300;
		CurrentPosition[PosZ] = 1; // always above table
	}

	//printf("delta t: %f.\n", deltaT);
	if (trial_flag == 1.0) {
		prevtime = sys_time;
		sys_time = clock() / (float)CLOCKS_PER_SEC - beg_sys_time;
	}
	deltaT = sys_time - prevtime;
	if (deltaT == 0.0) {
		deltaT = deltaTdefault;
	}
	deltaTvec[0] = deltaTvec[1]; deltaTvec[1] = deltaTvec[2]; deltaTvec[2] = deltaT;

	// If the trial hasn't started
	if (trial_flag == 0.0) {
		xprev[0] = CurrentPosition[PosX]; xprev[1] = CurrentPosition[PosX]; xprev[2] = CurrentPosition[PosX];
		yprev[0] = CurrentPosition[PosY]; yprev[1] = CurrentPosition[PosY]; yprev[2] = CurrentPosition[PosY];
		thetaXprev[0] = sys.Xcurr[0]; thetaXprev[1] = sys.Xcurr[0]; thetaXprev[2] = sys.Xcurr[0];
		thetaYprev[0] = sys.Xcurr[2]; thetaYprev[1] = sys.Xcurr[2]; thetaYprev[2] = sys.Xcurr[2];

		ballXprev[0] = BallPosition[PosX]; ballXprev[1] = BallPosition[PosX]; ballXprev[2] = BallPosition[PosX];
		ballYprev[0] = BallPosition[PosY]; ballYprev[1] = BallPosition[PosY]; ballYprev[2] = BallPosition[PosY];
	}

	// Update simulation once trial has started
	double xacc, yacc;
	xprev[0] = xprev[1]; xprev[1] = xprev[2]; xprev[2] = CurrentPosition[PosX];
	xacc = -(2.0 * xprev[1] - xprev[2] - xprev[0]) / (deltaTvec[1] * deltaTvec[2]);
	yprev[0] = yprev[1]; yprev[1] = yprev[2]; yprev[2] = CurrentPosition[PosY];
	yacc = -(2.0 * yprev[1] - yprev[2] - yprev[0]) / (deltaTvec[1] * deltaTvec[2]);

	if ((CurrentPosition[PosZ] > z_tolerance) && (ball_energy >= max_energy)) { // && (ball_energy < max_energy)
		sys.Ucurr = { xacc * factor,yacc * factor };
	} else {
		sys.Ucurr = { 0.0, 0.0 };
	}
	
	sys.simulate();
	double Z_potential = RR_vis - RR_vis * cos(sys.Xpotential[0]) * cos(sys.Xpotential[2]);
	if ((BallPosition[PosZ] > 0.8*RR_vis) && (Z_potential>BallPosition[PosZ])) { 
		//printf("ball z position: %f\n", BallPosition[PosZ]);
		sys.Xcurr[1] = 0.0* sys.Xcurr[1]; //0.1
		sys.Xcurr[3] = 0.0* sys.Xcurr[3]; //0.1
	}
	if (trial_flag == 1.0) {
		sys.dt = deltaT;
		sys.step();
	}
	//printf("u current: %f %f.\n", sys.Ucurr[0], sys.Ucurr[1]);
	//printf("x current: %f %f %f %f.\n", sys.Xcurr[0], sys.Xcurr[1], sys.Xcurr[2], sys.Xcurr[3]);

	if (mode == 0) // if using HapticMASTER
	{
		haSendCommand(dev, "get measforce", response);
		ParseFloatVec(response, CurrentForce[PosX], CurrentForce[PosY], CurrentForce[PosZ]);
	}
	else
	{
		CurrentForce[PosX] = 0;  CurrentForce[PosY] = 0;  CurrentForce[PosZ] = 0;
	}

	CheckFlags();

	//// Update forces from ball movement (option 1)
	//double ballXacc, ballYacc;
	//double thetaXacc, thetaYacc;
	//thetaXprev[0] = thetaXprev[1]; thetaXprev[1] = thetaXprev[2]; thetaXprev[2] = sys.Xcurr[0];
	//thetaXacc = -(2.0 * thetaXprev[1] - thetaXprev[2] - thetaXprev[0]) / (deltaTvec[1] * deltaTvec[2]);
	//thetaYprev[0] = thetaYprev[1]; thetaYprev[1] = thetaYprev[2]; thetaYprev[2] = sys.Xcurr[2];
	//thetaYacc = -(2.0 * thetaYprev[1] - thetaYprev[2] - thetaYprev[0]) / (deltaTvec[1] * deltaTvec[2]);
	//ballYacc = yacc - RR_vis * sin(sys.Xcurr[2]) * (cos(sys.Xcurr[0]) * (pow(sys.Xcurr[1], 2) + pow(sys.Xcurr[3], 2)) + sin(sys.Xcurr[0]) * thetaXacc) + RR_vis * cos(sys.Xcurr[2]) * (-2 * sin(sys.Xcurr[0]) * sys.Xcurr[1] * sys.Xcurr[3] + cos(sys.Xcurr[0]) * thetaYacc);
	//ballXacc = xacc - RR_vis * cos(sys.Xcurr[0]) * (-2 * sin(sys.Xcurr[2]) * sys.Xcurr[1] * sys.Xcurr[3] + cos(sys.Xcurr[2]) * thetaXacc) - RR_vis * sin(sys.Xcurr[0]) * (cos(sys.Xcurr[2]) * (pow(sys.Xcurr[1], 2) + pow(sys.Xcurr[3], 2)) + sin(sys.Xcurr[2]) * thetaYacc);
	
	// Update forces from ball movement (option 2)
	double ballXacc, ballYacc;
	ballXprev[0] = ballXprev[1]; ballXprev[1] = ballXprev[2]; ballXprev[2] = BallPosition[PosX];
	ballYprev[0] = ballYprev[1]; ballYprev[1] = ballYprev[2]; ballYprev[2] = BallPosition[PosY];
	ballXacc = -(2.0 * ballXprev[1] - ballXprev[2] - ballXprev[0]) / (deltaTvec[1] * deltaTvec[2]);
	ballYacc = -(2.0 * ballYprev[1] - ballYprev[2] - ballYprev[0]) / (deltaTvec[1] * deltaTvec[2]);
	//ballXacc = -2.0 * (deltaTvec[1] * (ballXprev[2] - ballXprev[1]) - deltaTvec[2] * (ballXprev[1] - ballXprev[0])) / (deltaTvec[1] * deltaTvec[2] * (deltaTvec[1] + deltaTvec[2]));
	//ballYacc = -2.0 * (deltaTvec[1] * (ballYprev[2] - ballYprev[1]) - deltaTvec[2] * (ballYprev[1] - ballYprev[0])) / (deltaTvec[1] * deltaTvec[2] * (deltaTvec[1] + deltaTvec[2]));

	//// Add damping in the x direction
	//float BaccX = 0.03;//0.03; // 0.05; // 0.015;
	//ballXacc = ballXacc - BaccX * (ballXacc - ballXaccPREV) / deltaTdefault;
	//float BaccY = 0.005;//0.005; // 0.05; // 0.015;
	//ballYacc = ballYacc - BaccY * (ballYacc - ballYaccPREV) / deltaTdefault;

	//// Add low-pass filter (option 1)
	//double ballXaccF, ballYaccF;
	//double RC = 1.0/(5.0*2*M_PI); // use cutoff freq here
	//float alpha = deltaT/(RC+deltaT);
	//ballXaccF = ballXaccPREV + alpha * (ballXacc - ballXaccPREV);
	//ballYaccF = ballYaccPREV + alpha * (ballYacc - ballYaccPREV);
	//ballXaccPREV = ballXaccF;
	//ballYaccPREV = ballYaccF;
	//ballXaccPREV = ballXacc;
	//ballYaccPREV = ballYacc;

	//// Add low-pass filter (option 2)
	//double ballXaccF, ballYaccF;
	//ballXaccVEC[0] = ballXaccVEC[1]; ballXaccVEC[1] = ballXaccVEC[2]; ballXaccVEC[2] = ballXacc;
	//ballYaccVEC[0] = ballYaccVEC[1]; ballYaccVEC[1] = ballYaccVEC[2]; ballYaccVEC[2] = ballYacc;
	//ballXaccF = (ballXaccVEC[0] + ballXaccVEC[1] + ballXaccVEC[2]) / 3.0;
	//ballYaccF = (ballYaccVEC[0] + ballYaccVEC[1] + ballYaccVEC[2]) / 3.0;

	//// Add scaling and cap to the forces
	//double feedback_cap = 4.0; double feedback_scaling_x = 0.4 * pow(R / R_vis, 1); double feedback_scaling_y = 0.6 * pow(R / R_vis, 1); //0.5
	//double feedback_cap = 4.0; double feedback_scaling_x = 0.3 * pow(R/ R_vis,0.8); double feedback_scaling_y = 0.5 * pow(R/R_vis,0.8); //0.5
	//ballXaccF = -ballXaccF * feedback_scaling_x; ballYaccF = ballYaccF * feedback_scaling_y;
	//// option 1: if filtering off
	//if (ballXacc > feedback_cap) { ballXacc = feedback_cap; }
	//else if (ballXacc < -feedback_cap) { ballXacc = -feedback_cap; }
	//if (ballYacc > feedback_cap) { ballYacc = feedback_cap; }
	//else if (ballYacc < -feedback_cap) { ballYacc = -feedback_cap; }
	//// option 2: if filtering on
	//if (ballXaccF > feedback_cap) { ballXaccF = feedback_cap; }
	//else if (ballXaccF < -feedback_cap) { ballXaccF = -feedback_cap; }
	//if (ballYaccF > feedback_cap) { ballYaccF = feedback_cap; }
	//else if (ballYaccF < -feedback_cap) { ballYaccF = -feedback_cap; }

	
	// Lower forces in the x direction to avoid oscillations
	double ballXaccFactual;
	double ballYaccFactual;
	float ball_mag = sqrt(pow(ballXacc, 2) + pow(ballYacc, 2));
	if ((abs(ballXacc) > 0.005) && (abs(ballYacc) > 0.005)) {
		ballXaccFactual = -2 * ballXacc / ball_mag;
		ballYaccFactual = 2 * ballYacc / ball_mag;

		if (abs(ballXaccFactual) > 1.75) {
			ballXaccFactual = -.5 * ballXacc / ball_mag;
			ballYaccFactual = .5 * ballYacc / ball_mag;
		}

		ball_mag = sqrt(pow(ballXaccFactual, 2) + pow(ballYaccFactual, 2));
		//printf("ball_mag %f \n", ball_mag);
	} else {
		ballXaccFactual = 0.0;
		ballYaccFactual = 0.0;
	}

	if ((ball_energy < 1.5* mid_energy) && (mode == 0)) { //|| (ball_energy > max_energy*4)) 
		ballXaccFactual = -.5 * ballXacc / ball_mag;
		ballYaccFactual = .5 * ballYacc / ball_mag;
	} 
	
	if ((mode == 0) && (feedback_forces == 1)) // if using HapticMASTER
	{
		//printf("forces %f %f \n", ballXaccFactual, ballYaccFactual);
		haSendCommand(dev, "set BallFeedbackForce force", ballXaccFactual, ballYaccFactual, 0.0, response);
	}

	// Write data to .txt file
	if (trial_flag == 1.0) {
		logfile << sys_time << "," << CurrentPosition[PosX] << "," << CurrentPosition[PosY] << "," << CurrentPosition[PosZ] << ",";
		logfile << BallPosition[PosX] << "," << BallPosition[PosY] << ",";
		logfile << BallPosition[PosZ] << "," << score << ",";
		logfile << CurrentForce[PosX] << "," << CurrentForce[PosY] << "," << CurrentForce[PosZ] << "," << scoring_enabled << "," << ballXacc << "," << ballYacc << ",";
		logfile << ballXaccFactual << "," << ballYaccFactual << "," << ball_energy << "\n";
	}

	// Set time of trial 
	if (sys_time > TF) 
	{
		logfile.close();

		if (score > pb_vec[freq_num][condition] && personalbests_flag == 1) {
			pb_vec[freq_num][condition] = score;

			ofstream myfile;
			myfile.open(personalbestspath);
			for (int i = 0; i < num_freqs_tested; i++) {
				myfile << pb_vec[i][0] << ',';
			}
			myfile << '\n';
			for (int i = 0; i < num_freqs_tested; i++) {
				myfile << pb_vec[i][1] << ',';
			}
			// myfile << '\n';
			// for (int i = 0; i < num_freqs_tested; i++) {
			// 	myfile << pb_vec[i][2] << ',';
			// }
			myfile.close();

			printf("-------------Participant beat their personal best!----------------\n");

			glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
			glPushMatrix();
			glutPostRedisplay();
			renderer->DrawTrophy();
			DrawPersonalBestBoard;
			DrawPersonalBestLabel;
			glPopMatrix();
			glutSwapBuffers();
			double curr_time = clock() / (float)CLOCKS_PER_SEC;
			while (clock() / (float)CLOCKS_PER_SEC - curr_time < 2.0) {}
		}

		if (mode == 0) // if using HapticMASTER
		{
			// Move the person back to their home position by enabling springs
			printf("Moving to home position\n");
			haDeviceSendString(dev, "set mySpring enable", response);

			// Disable springs when ready 
			/*printf("Press ENTER to disable spring\n");
			cin.get();*/
			float pos_tol = 0.01;
			while ((abs(CurrentPosition[PosX] - springPos[PosX]) > pos_tol) || (abs(CurrentPosition[PosY] - springPos[PosY]) > pos_tol) || (abs(CurrentPosition[PosZ] - springPos[PosZ] - 0.011) > pos_tol)) {
				haSendCommand(dev, "get measpos", response);
				ParseFloatVec(response, CurrentPosition[PosX], CurrentPosition[PosY], CurrentPosition[PosZ]);
			}
			haDeviceSendString(dev, "set mySpring disable", response);

			// Closing routine
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
	}

	// sleeper function
	timer_val = clock() / (float)CLOCKS_PER_SEC - timer_beg;
	while (timer_val < deltaTdefault) {	
		timer_val = clock() / (float)CLOCKS_PER_SEC - timer_beg;
	}
	timer_beg = clock() / (float)CLOCKS_PER_SEC;

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
	int returnValue;

	srand(static_cast <unsigned> (time(0))); // this is important for random number generation

	FILE* setupfile;
	printf(setuppath);
	setupfile = fopen(setuppath, "r");
	int val2 = fscanf(setupfile, "%f,%f,%f,%f,%f,%f\n", &armWeight, &ymin, &ymax, &xmin, &xmax, &maxForce);

	FILE* pbfile;
	pbfile = fopen(personalbestspath, "r");
	fscanf(pbfile, "%d,%d,%d,%d\n,%d,%d,%d,%d\n", &pb_vec[0][0], &pb_vec[1][0], &pb_vec[2][0], &pb_vec[3][0], &pb_vec[0][1], &pb_vec[1][1], &pb_vec[2][1], &pb_vec[3][1]);
	printf("Personal best for this frequency is: %d.\n", pb_vec[freq_num][condition]);


	// Option to hardcode values
	// armWeight = 20.0; ymin = -0.3; ymax = 0.3; xmin = -0.1; xmax = 0.19;

	if (armWeight == 10.0 || ymin == 10.0 || ymax == 10.0 || xmin == 10.0 || xmax == 10.0 || maxForce == 10.0)
	{
		printf("error reading setupfile, hardcoded values used - max workspace\n");
		armWeight = 60.0; ymin = -0.25; ymax = 0.25; xmin = -0.1; xmax = 0.15; maxForce = 30.0;
		fclose(setupfile);
	}

	fclose(setupfile);

	springPos[0] = (xmin + xmax) / 2; springPos[1] = (ymin + ymax) / 2; springPos[2] = table_z; //-.17
	if (support_num == 0.0) {
		springPos[2] = table_z + 0.045; //-.17
	}
	else {
		springPos[2] = table_z; //-.17
	}
	double springStiffness = 4000.0, springDamping = 10.0, springMaxForce = 100.0;

	DefineTask();

	// Open file to log state data
	logfile.open(logpath);
	logfile << "system_time, x, y, z, xball, yball, zball, flags,fx,fy,fz, scoring_enabled, ball_fx, ball_fy, fil_ball_fx, fil_ball_fy, ball energy \n";

	if (mode == 0) // if using HapticMASTER
	{

		// Open the HapticMaster device
		dev = haDeviceOpen(IPADDRESS);

		// Check if device connected, then initialize
		if (dev == HARET_ERROR) {
			printf("--- ERROR: Unable to connect to device: %s\n", IPADDRESS);
			//return HARET_ERROR; // comment this out to test visualization without the HapticMaster
		}
		else {
			InitializeDevice(dev);
			
		}

		haSendCommand(dev, "set inertia 4.0", response);

		// Fix arm in place using springs at the start position
		printf("Moving...\n");
		returnValue = haDeviceSendString(dev, "create spring mySpring", response);
		haSendCommand(dev, "set mySpring pos", springPos[PosX], springPos[PosY], springPos[PosZ] + 0.011, response);
		haSendCommand(dev, "set mySpring stiffness", springStiffness, response);
		haSendCommand(dev, "set mySpring dampfactor", springDamping, response);
		haSendCommand(dev, "set mySpring maxforce", springMaxForce, response);
		returnValue = haDeviceSendString(dev, "set mySpring enable", response);

		// Set bias force
		printf("Setting bias force to %f of the max abduction force.\n", support_level[support_num]);
		float bias = - maxForce * support_level[support_num] + armWeight; // changed to subtract arm weight
		returnValue = haDeviceSendString(dev, "create biasforce myBiasForce", response);
		printf("bias: %f \n", bias);
		printf("Weight: %f \n", armWeight);
		printf("Max: %f \n", maxForce);
		returnValue = haSendCommand(dev, "set myBiasForce force", 0.0, 0.0, bias, response);
		printf("Bias force: %s\n", response);

		// Wait to start visualization/trial
		// printf("Press ENTER once arm in position.\n");
		// cin.get();
		// Disable springs when ready 
		double pos_tol = 0.01;
		while ((abs(CurrentPosition[PosX] - springPos[PosX]) > pos_tol) || (abs(CurrentPosition[PosY] - springPos[PosY]) > pos_tol) || (abs(CurrentPosition[PosZ] - springPos[PosZ] - 0.011) > pos_tol)) {
			haSendCommand(dev, "get measpos", response);
			ParseFloatVec(response, CurrentPosition[PosX], CurrentPosition[PosY], CurrentPosition[PosZ]);
		}

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

		// Initialize start position for ball simulation
		haSendCommand(dev, "get measpos", response);
		ParseFloatVec(response, CurrentPosition[PosX], CurrentPosition[PosY], CurrentPosition[PosZ]);
		CurrentPosition[PosX] = 0; CurrentPosition[PosY] = 0; CurrentPosition[PosZ] = 0; // should this be commented out?
		
		// Add damping in the x direction
		// haSendCommand(dev, "create damper myDamper", response);
		// haSendCommand(dev, "set myDamper dampcoef", 50.0, 0.0, 0.0, response);
		// haSendCommand(dev, "set myDamper enable", response);
	}
	else
	{
		// Set initial conditions
		CurrentPosition[PosX] = 0; CurrentPosition[PosY] = 0; CurrentPosition[PosZ] = 0;
	}

	float rand_num = get_random() - 0.5;
	float sign = rand_num / abs(rand_num);
	float x_ang = sign * get_random()* angle_disturbance; // change to be both positive and negative
	rand_num = get_random() - 0.5;
	sign = rand_num / abs(rand_num);
	float y_ang = sign * asin(sqrt(pow(sin(angle_disturbance), 2) - pow(sin(x_ang), 2)));
	//float y_ang = sign*sqrt(pow(angle_disturbance, 2) - pow(x_ang, 2));
	printf("sign: %f x_angle: %f y_angle: %f \n", sign, x_ang, y_ang);

	sys.Xcurr = { x_ang,0.0,y_ang,0.0,CurrentPosition[PosX],0.0,CurrentPosition[PosY],0.0 }; // M_PI / 3.0

	printf("Visualization starting...\n");

	// OpenGL Initialization Calls
	glutInit(&argc, argv);
	glfwInit();
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(SCREEN_WIDTH, SCREEN_HEIGHT);	
	
	// Create The OpenGlWindow and Initialize Things
	glutCreateWindow("Ball in Bowl");
	InitOpenGl();
	glutReshapeFunc(Reshape);
	glutDisplayFunc(Display);
	glutKeyboardFunc(Keyboard);

	#ifdef TEXTURES_ON
	renderer = new Renderer;
	renderer->init();
	#endif

	if (mode == 0) // if using HapticMASTER
	{
		// Start Ball Feedback
		returnValue = haDeviceSendString(dev, "create biasforce BallFeedbackForce", response);
		returnValue = haSendCommand(dev, "set BallFeedbackForce force [0.0, 0.0, 0.0]", response);
		haDeviceSendString(dev, "set BallFeedbackForce enable", response);
	}
	else
	{
		// Check for joystick
		if ( !glfwJoystickPresent( GLFW_JOYSTICK_1 ) )
		{
			std::cout << "ERROR: No Joystick Connected " << std::endl;	
			exit(-1);
		}
	}

	timer_beg = clock() / (float)CLOCKS_PER_SEC;
	glutTimerFunc(10, TimerCB, 1);
	glutMainLoop();
	glfwTerminate();
	
	return 0;
}