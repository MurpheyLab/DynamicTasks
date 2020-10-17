//---------------------------------------------------------------------
// ballinbowl_full.cpp
// Interactive bowl visualization (with textures and water) moved using 
// either the HapticMASTER or a game controller.
//---------------------------------------------------------------------
#include "ballinbowl.hpp"
#include "parameters.hpp"
#include "texture.h"
#include "HapticAPI2.h"
#include "HapticUtility.h"
#include <glut.h>
#include <glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h> // sine, cosine
#include "cmath"
#include <cstdlib>
#include <iostream> // cin
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <ctime>
using namespace std;

//Screen dimension constants
const int SCREEN_WIDTH = 6000;
const int SCREEN_HEIGHT = 3000;

//this is a static pointer to a Renderer used in the glut callback functions
static Renderer *renderer;

//---------------------------------------------------------------------
//---------------------------------------------------------------------

float x_scale = 5.0; //used only for visualization

// Define possible support levels
float support_level[3] = { -0.0, -0.2, -0.5 }; // fraction of max shoulder abduction loading (0-1)

// Define possible tasks & setup variables
const char* tasklist[6] = { taskpath0, taskpath1, taskpath2, taskpath3, taskpath4, taskpath5};

// Define file for logging data
const char* taskpath = tasklist[task_num];
string logpath = "OutputData_S" + to_string(subject_num) + "_Task" + to_string(task_num) + "_SL" + to_string(support_num) + "_Trial" + to_string(trial_num) + ".csv";

ofstream logfile;

#define IPADDRESS "10.30.203.36" //"10.30.203.26"

#define PosX 0
#define PosY 1
#define PosZ 2

#define deltaT 0.01
#define R 0.07
#define factor 20.0 // 60.0 without feedback forces
#define z_tolerance -0.16 // if lower, person can eat flags while slacking
#define tol_delta 0.008
#define size_of_flag 0.005

static GLUquadricObj * q; // water object

//GLfloat xmin = -0.08, xmax = 0.15, ymin = -0.2, ymax = 0.25;
GLfloat xmin = 10.0; GLfloat xmax = 10.0; GLfloat ymin = 10.0; GLfloat ymax = 10.0;
GLfloat armWeight = 10.0;

long dev = 0;
char response[100];
double CurrentPosition[3];
double CurrentJoystick[2];
double CurrentForce[3];
GLfloat BallPosition[3];
char ball_direction, ball_direction_temp;
int ball_intensity, ball_intensity_temp;
GLfloat goals[num_csv_flags][3], goals_rand[num_rand_flags][3];
int goal_status[num_csv_flags], goal_rand_status[num_rand_flags];
int score = 0;

float time_current = 0.0;
float trial_flag = 0.0;
float scoring_enabled = 0.0;

float new_g = 2.76;
//BallBowl sys(1.0, -0.005, -9.81, R, deltaT); // m, B, g, radius, dt in seconds (Original 1.88 Hz)
//BallBowl sys(1.0, -0.0005, -44.22, R, deltaT); // 4 Hz
BallBowl sys(1.0, -0.015, -2.76, R, deltaT); // 1 Hz
//BallBowl sys(1.0, -0.001, -24.87, R, deltaT); // 3 Hz
//BallBowl sys(1.0, -0.03, -0.69, R, deltaT); // 0.5 Hz
double xprev[3] = {0.0,0.0,0.0}; double yprev[3] = {0.0,0.0,0.0};
double thetaXprev[3] = {0.0,0.0,0.0}; double thetaYprev[3] = {0.0,0.0,0.0};

// define ball
#define radius 0.01
#define RR (R - radius)
bool penalty = false;
clock_t penalty_timer, drop_timer;
bool reset_penalty_timer = true, reset_drop_timer = true;
double penalty_duration, drop_duration;

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
//                     D R A W   B A L L
//
// This Function Is Called To Draw The Ball
//---------------------------------------------------------------------
void DrawBall(void)
{
	BallPosition[PosX] = CurrentPosition[PosX] + RR*sin(sys.Xcurr[0]) * cos(sys.Xcurr[2]);
	BallPosition[PosY] = CurrentPosition[PosY] + RR*sin(sys.Xcurr[2]) * cos(sys.Xcurr[0]);
	BallPosition[PosZ] = RR - RR * cos(sys.Xcurr[0]) * cos(sys.Xcurr[2]); //atan(pow(pow(tan(sys.Xcurr[0]), 2) + pow(tan(sys.Xcurr[2]), 2), 0.5
	
	if (BallPosition[PosZ] > RR/5 && BallPosition[PosZ] <= 0.4*RR)
	{
		glColor3f(1.0f, 0.5f, 0.0f); // orange
		ball_intensity = 1; // DON'T CHANGE   ball_intensity is on a range from 0 - 2
	} 
	else if (BallPosition[PosZ] > 0.4*RR)
	{
		glColor3f(1.0f, 0.0f, 0.0f); // red
		ball_intensity = 2; // DON'T CHANGE

		// penalty flag is true if timer is reset
		if (reset_penalty_timer)
		{
			penalty_timer = clock();
			penalty = true;	
			reset_penalty_timer = false;		
		}
		else
		{
			// reset timer if given number of seconds has passed
			penalty_duration = (clock() - penalty_timer) / (double) CLOCKS_PER_SEC;
			if (penalty_duration > 1.0) // ADJUSTABLE duration of penalty forgiveness timer in seconds 
			{
				reset_penalty_timer = true;
			}
		}
	}
	else
	{
		glColor3f(1.0f, 1.0f, 0.0f); // yellow
		ball_intensity = 0; // DON'T CHANGE
	}

	glPushMatrix();
	glTranslatef(BallPosition[PosX], BallPosition[PosY], BallPosition[PosZ] + radius);
	glutSolidSphere(radius, 20, 20);
	glPopMatrix();
}

//---------------------------------------------------------------------
//                     D R A W   B O W L
//
// This Function Is Called To Draw The Bowl
//---------------------------------------------------------------------
void DrawBowl(void) {
	int i, j;
	GLfloat r = R;
	// change these to render more or less circles
	const int scalex = 30; 
	const int scaley = 30;
	GLfloat v[scalex * scaley][3];

	for (i = 0; i < scalex; ++i) {
		for (j = 0; j < scaley; ++j) {
			v[i * scaley + j][0] = CurrentPosition[PosX] + r * cos(j * 2 * M_PI / scaley) * cos(i * M_PI / (2 * scalex));
			v[i * scaley + j][1] = CurrentPosition[PosY] + r * sin(j * 2 * M_PI / scaley) * cos(i * M_PI / (2 * scalex));
			v[i * scaley + j][2] = r - r * sin(i * M_PI / (2 * scalex));
		}
	}

	if (scoring_enabled==1.0) {
		glColor3f(0.3f, 0.0f, 0.51f);
	}
	else {
		glColor3f(0.0f, 0.0f, 0.0f);
	}

	glBegin(GL_LINES);
	for (i = 10; i < scalex - 1; ++i) { // change the penalty_timer condition for i to make the bowl shorter or taller (was i=1 for full bowl)
		for (j = 0; j < scaley; ++j) {
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
	GLfloat r = 0.06; // ADJUSTABLE radius of water in bowl
	const GLint numberOfSides = 360; 
	const GLint numberOfVertices = numberOfSides + 2;	

	GLfloat circleVerticesX[numberOfVertices];
	GLfloat circleVerticesY[numberOfVertices];
	GLfloat circleVerticesZ[numberOfVertices];

	// DON'T CHANGE must all be 0
	circleVerticesX[0] = 0; 
	circleVerticesY[0] = 0; 
	circleVerticesZ[0] = 0;

	for (int i = 1; i < numberOfVertices; i++)
	{
		circleVerticesX[i] = 0.05 * cos(i * 2 * M_PI / numberOfSides);
		circleVerticesY[i] = 0.05 * sin(i * 2 * M_PI / numberOfSides);
		circleVerticesZ[i] = 0;
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
	glTranslatef(CurrentPosition[PosX], CurrentPosition[PosY], r);
	glRotatef(ball_angle_y, 0.0, 1.0, 0.0);
	glRotatef(ball_angle_x, 1.0, 0.0, 0.0);
	glEnable(GL_CLIP_PLANE0);
	double clipEq[4] = { 0.0, 0.0, -1.0, -0.03 }; // ADJUSTABLE plane that removes the top half of the sphere of water
   	glClipPlane(GL_CLIP_PLANE0, clipEq);
	gluSphere(q, r, 16, 16);
	glTranslatef(0.0, 0.0, -r/2);
	glDisable(GL_CLIP_PLANE0);

	// draw water surface
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
			//printf("%f,%f,%f\n", x_temp, y_temp, z_temp);
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
		//generates x and y positions from 0.0 to 1.0 inclusively
		GLfloat rand_x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
		GLfloat rand_y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

		//scale based on reachable workspace
		goals_rand[0][PosX] = xmax * rand_x + xmin * (1 - rand_x);
		goals_rand[0][PosY] = ymax * rand_y + ymin * (1 - rand_y);
		goals_rand[0][PosZ] = 0.0;

		for (int i = 1; i < num_rand_flags; i++)
		{
			GLfloat min_dist = 0.0;

			//flags cannot be placed less than 0.08 units from another flag
			while (min_dist < 0.08) // ADJUSTABLE: minimum distance between two randomly placed flags
			{
				//generates x and y positions from 0.0 to 1.0 inclusively
				rand_x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
				rand_y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

				//scale based on reachable workspace
				goals_rand[i][PosX] = xmax * rand_x + xmin * (1 - rand_x);
				goals_rand[i][PosY] = ymax * rand_y + ymin * (1 - rand_y);
				goals_rand[i][PosZ] = 0.0;

				//check that flag is not too close to the previously created flags
				GLfloat distances[num_rand_flags];
				for (int j = 0; j < i; j++)
				{
					distances[j] = sqrt(pow(goals_rand[i][PosX] - goals_rand[j][PosX], 2) + pow(goals_rand[i][PosY] - goals_rand[j][PosY], 2));
				}

				min_dist = distances[0];

				for (int k = 0; k < i; k++)
				{
					if (min_dist > distances[k])
					{
						min_dist = distances[k];
					}
				}
			} 
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

	float zBall = R - R * cos(sys.Xcurr[0]) * cos(sys.Xcurr[2]); // only allows flags to disappear if ball within bowl

	// here we check for: 1. ball is in bowl, 2. trial is running, 3. person is above the haptic table
	if ((zBall < 0.4 * RR) && (trial_flag == 1.0) && (CurrentPosition[PosZ] > z_tolerance)) {
		scoring_enabled = 1.0;
	}
	else {
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
			if ((abs(CurrentPosition[PosX] - goals[i][PosX] + 0.015) < 3.0*tol_delta) && (abs(CurrentPosition[PosY] - goals[i][PosY]) < tol_delta) && (goal_status[i]==0) && (scoring_enabled==1.0))
			
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
			if ((abs(CurrentPosition[PosX] - goals_rand[i][PosX] + 0.015) < 3.0*tol_delta) && (abs(CurrentPosition[PosY] - goals_rand[i][PosY]) < tol_delta) && (goal_rand_status[i]==0) && (scoring_enabled==1.0))
			{
				// replace eaten flag with new, randomly placed flag
				GLfloat min_dist = 0.0;

				while (min_dist < 0.08) //ADJUSTABLE minimum distance between two randomly placed flags
				{
					//generates x and y positions from 0.0 to 1.0 inclusively
					rand_x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
					rand_y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

					//scale based on reachable workspace
					goals_rand[i][PosX] = xmax * rand_x + xmin * (1 - rand_x);
					goals_rand[i][PosY] = ymax * rand_y + ymin * (1 - rand_y);

					//check that flag is not too close to the previously created flags
					GLfloat distances[num_rand_flags];
					for (int j = 0; j < num_rand_flags; j++)
					{
						distances[j] = sqrt(pow(goals_rand[i][PosX] - goals_rand[j][PosX], 2) + pow(goals_rand[i][PosY] - goals_rand[j][PosY], 2));
					}

					min_dist = 1000; // DON'T CHANGE  arbitrarily high number

					for (int k = 0; k < num_rand_flags; k++)
					{
						if (k == i)
						{
							continue;
						}
						else if (min_dist > distances[k])
						{
							min_dist = distances[k];
						}
					}
				} 
				score++;
			}
		}
	}
	else if (game_version == 'p')
	{
		for (int i = 0; i < num_rand_flags; i++) 
		{
			// because of the perspective of the camera, the x position is given a greater tolerance
			if ((abs(CurrentPosition[PosX] - goals_rand[i][PosX] + 0.015) < 3.0*tol_delta) && (abs(CurrentPosition[PosY] - goals_rand[i][PosY]) < tol_delta) && (goal_rand_status[i]==0) && (scoring_enabled==1.0))
			{
				// eats flags when only x and y overlap
				{
					goal_rand_status[i] = 1; // set achieved goal to disappear (turn black)
					score++;
				}
			} 
		}

		if (penalty)
		{   
			for (int i = 0; i < num_rand_flags; i++)
			{
				if (goal_rand_status[i] == 1)
				{
					goal_rand_status[i] = 0;

					// replace eaten flag with new, randomly placed flag
					GLfloat min_dist = 0.0;

					while (min_dist < 0.08) //ADJUSTABLE minimum distance between two randomly placed flags
					{
						//generates x and y positions from 0.0 to 1.0 inclusively
						rand_x = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
						rand_y = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);

						//scale based on reachable workspace
						goals_rand[i][PosX] = xmax * rand_x + xmin * (1 - rand_x);
						goals_rand[i][PosY] = ymax * rand_y + ymin * (1 - rand_y);

						//check that flag is not too close to the previously created flags
						GLfloat distances[num_rand_flags];
						for (int j = 0; j < num_rand_flags; j++)
						{
							distances[j] = sqrt(pow(goals_rand[i][PosX] - goals_rand[j][PosX], 2) + pow(goals_rand[i][PosY] - goals_rand[j][PosY], 2));
						}

						min_dist = 1000; // DON'T CHANGE  arbitrarily high number

						for (int k = 0; k < num_rand_flags; k++)
						{
							if (k == i)
							{
								continue;
							}
							else if (min_dist > distances[k])
							{
								min_dist = distances[k];
							}
						}
					} 

					score--;
					penalty = false;
					break;
				}
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
	char *cscore = const_cast<char*>(score_string.c_str());

	int time_current_rounded = round(TF - time_current);
	string time_current_string = to_string(time_current_rounded);
	char *ctime = const_cast<char*>(time_current_string.c_str());

    //glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	// thicken the text
	for( float i = 0.0; i <= 0.002; i = i + 0.0001) 
    {
		glColor3f(0.0f, 0.7f, 1.0f);
        DrawText(-0.44 + i, 0.19, "SCORE");

		if (score < 10) {
			DrawText(-0.385 + i, 0.14, cscore);
		}
		else {
			DrawText(-0.395 + i, 0.14, cscore);
		}

		glColor3f(1.0f, 0.0f, 0.0f);
		DrawText(0.325 + i, 0.19, "TIME");

		if (time_current_rounded < 10) {
			DrawText(0.355 + i, 0.14, ctime);
		}
		else {
			DrawText(0.345 + i, 0.14, ctime);
		}
    }
}

//---------------------------------------------------------------------
//                     D R A W   T I M E R   B A R
//
// This Function Is Called To Draw A Timer Bar
//---------------------------------------------------------------------
void DrawTimerBar(void) {
	int i, j;
	GLfloat v[4][3];
	GLfloat bottom = 0.07, width = 0.053, x_offset = 0.3, length = 0.89 * (1 - time_current / TF);

	for (i = 0; i < 4; ++i) { // parallel with x-y plane
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
void DrawScoreBar(void) {
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
	if (scoring_enabled==1.0) {
		glColor3f(0.01f, 0.5f, 1.0f);
	}
	else {
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
	renderer->DrawTable();
	DrawBall();
	DrawBowl();
	
	DrawTimerBar();

	// no score bar for infinite game version
	if (game_version == 'f' || game_version == 'p')
	{
		DrawScoreBar();
	}

	DrawBowlBottom();

	if (game_version == 'f')
	{
		renderer->DrawFlags(goals, goal_status);
	}
	else
	{
		renderer->DrawRandomFlags(goals_rand, goal_rand_status);
	}
	
	DrawWater();
	DrawLabels();

	// water drop flag is true if timer is reset
	if (reset_drop_timer)
	{
		drop_timer = clock();
		reset_drop_timer = false;
		ball_direction_temp = ball_direction;
		ball_intensity_temp	 = ball_intensity;		
	}
	else
	{
		// reset timer after a set time has passed
		drop_duration = (clock() - drop_timer) / (double) CLOCKS_PER_SEC;
		if (drop_duration > 0.13)
		{
			reset_drop_timer = true;
		}
	}

	renderer->DrawWaterDrops(CurrentPosition, ball_direction_temp, ball_intensity_temp, drop_duration, R);

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
	switch (ucKey)
	{
	case 's': // penalty_timer trial
		// Allow movement and penalty_timer timer function
		printf("S pressed: Trial and timer starting...\n");
		time_current = 0.0;
		trial_flag = 1.0;
		if (mode == 0) // if using HapticMASTER
		{
			haSendString(dev, "set myBiasForce enable", response);
			haSendString(dev, "set mySpring disable", response);
		}
		break;
	case 'b': // enable bias force
		printf("B pressed \n");
		haSendString(dev, "set myBiasForce enable", response);
		printf("%s\n", response);
		break;
	case 'd': // disable bias force
		printf("D pressed \n");
		haSendString(dev, "set myBiasForce disable", response);
		printf("%s\n", response);
		break;
	case 27: // ESC key pressed - stop motion and simulation
		printf("ESC pressed: exit routine starting \n");
		logfile.close();
		if (mode == 0) // if using HapticMASTER
		{
			haSendCommand(dev, "remove all", response);
			printf("remove all ==> %s\n", response);
			haSendCommand(dev, "set state stop", response);
			printf("set state stop ==> %s\n", response);
			if (haCloseDevice(dev)) {
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
  	float offset = 1.0f;
  	renderer->t += offset;

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

	// If the trial hasn't started
	if (trial_flag == 0.0) 
	{
		xprev[0] = CurrentPosition[PosX]; xprev[1] = CurrentPosition[PosX]; xprev[2] = CurrentPosition[PosX];
		yprev[0] = CurrentPosition[PosY]; yprev[1] = CurrentPosition[PosY]; yprev[2] = CurrentPosition[PosY];
		thetaXprev[0] = sys.Xcurr[0]; thetaXprev[1] = sys.Xcurr[0]; thetaXprev[2] = sys.Xcurr[0];
		thetaYprev[0] = sys.Xcurr[2]; thetaYprev[1] = sys.Xcurr[2]; thetaYprev[2] = sys.Xcurr[2];
	}

	// Update simulation once trial has started
	if (trial_flag == 1.0) {
		time_current = time_current + deltaT;
	}

	double xacc, yacc;
	xprev[0] = xprev[1]; xprev[1] = xprev[2]; xprev[2] = CurrentPosition[PosX];
	xacc = -(2.0 * xprev[1] - xprev[2] - xprev[0]) / (deltaT * deltaT);
	yprev[0] = yprev[1]; yprev[1] = yprev[2]; yprev[2] = CurrentPosition[PosY];
	yacc = -(2.0 * yprev[1] - yprev[2] - yprev[0]) / (deltaT * deltaT);
	sys.Ucurr = { xacc * factor,yacc * factor };  
	sys.step();

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

	// Update forces from ball movement
	double ballXacc, ballYacc;
	double thetaXacc, thetaYacc;
	thetaXprev[0] = thetaXprev[1]; thetaXprev[1] = thetaXprev[2]; thetaXprev[2] = sys.Xcurr[0];
	thetaXacc = -(2.0 * thetaXprev[1] - thetaXprev[2] - thetaXprev[0]) / (deltaT * deltaT);
	thetaYprev[0] = thetaYprev[1]; thetaYprev[1] = thetaYprev[2]; thetaYprev[2] = sys.Xcurr[2];
	thetaYacc = -(2.0 * thetaYprev[1] - thetaYprev[2] - thetaYprev[0]) / (deltaT * deltaT);
	ballYacc = yacc - sys.h * sin(sys.Xcurr[2]) * (cos(sys.Xcurr[0]) * (pow(sys.Xcurr[1], 2) + pow(sys.Xcurr[3], 2)) + sin(sys.Xcurr[0]) * thetaXacc) + sys.h * cos(sys.Xcurr[2]) * (-2 * sin(sys.Xcurr[0]) * sys.Xcurr[1] * sys.Xcurr[3] + cos(sys.Xcurr[0]) * thetaYacc);
	ballXacc = yacc - sys.h * cos(sys.Xcurr[0]) * (-2 * sin(sys.Xcurr[2]) * sys.Xcurr[1] * sys.Xcurr[3] + cos(sys.Xcurr[2]) * thetaXacc) - sys.h * sin(sys.Xcurr[0]) * (cos(sys.Xcurr[2]) * (pow(sys.Xcurr[1], 2) + pow(sys.Xcurr[3], 2)) + sin(sys.Xcurr[2]) * thetaYacc);

	float feedback_cap = 4.0; float feedback_scaling_x = 0.5 * (9.81/new_g); float feedback_scaling_y = 0.5 * (9.81 / new_g);
	ballXacc = -ballXacc * feedback_scaling_x; ballYacc = ballYacc * feedback_scaling_y;

	if (ballXacc > feedback_cap) { ballXacc = feedback_cap; }
	else if (ballXacc < -feedback_cap) { ballXacc = -feedback_cap; }
	if (ballYacc > feedback_cap) { ballYacc = feedback_cap; }
	else if (ballYacc < -feedback_cap) { ballYacc = -feedback_cap; }
	
	if (mode == 0) // if using HapticMASTER
	{
		haSendCommand(dev, "set BallFeedbackForce force", ballXacc, ballYacc, 0.0, response);
	}

	// Write data to .txt file
	logfile << sys.tcurr << "," << CurrentPosition[PosX] << "," << CurrentPosition[PosY] << "," << CurrentPosition[PosZ] << ",";
	logfile << CurrentPosition[PosX] + RR * sin(sys.Xcurr[0]) * cos(sys.Xcurr[2]) << "," << CurrentPosition[PosY] + RR * sin(sys.Xcurr[2]) * cos(sys.Xcurr[0]) << ",";
	logfile << R - R * cos(sys.Xcurr[0]) * cos(sys.Xcurr[2]) << "," << score <<",";
	logfile << CurrentForce[PosX] << "," << CurrentForce[PosY] << "," << CurrentForce[PosZ] << "," << scoring_enabled << "," << ballXacc << "," << ballYacc << "\n";

	// Set time of trial 
	if (time_current > TF) {
		logfile.close();
		if (mode == 0) // if using HapticMASTER
		{
			haSendCommand(dev, "remove all", response);
			printf("remove all ==> %s\n", response);
			haSendCommand(dev, "set state stop", response);
			printf("set state stop ==> %s\n", response);
			if (haCloseDevice(dev)) {
				printf("--- ERROR: closing device\n");
			}
		}
		exit(0);
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
	int returnValue;
	double springPos[3] = { -0.05,0.05,-0.17 }; // x was 0.0
	double springStiffness = 4000.0, springDamping = 10.0, springMaxForce = 100.0;

	srand(static_cast <unsigned> (time(0)));

	FILE* setupfile;
	setupfile = fopen(setuppath, "r");
	int val2 = fscanf(setupfile, "%f,%f,%f,%f,%f\n", &armWeight, &ymin, &ymax, &xmin, &xmax);


	// Option to hardcode values
	//armWeight = 60.0; ymin = -0.3; ymax = 0.3; xmin = -0.1; xmax = 0.19;

	if (armWeight == 10.0 || ymin == 10.0 || ymax == 10.0 || xmin == 10.0 || xmax == 10.0)
	{
		printf("error reading setupfile, hardcoded values used - max workspace\n");
		armWeight = 60.0; ymin = -0.25; ymax = 0.25; xmin = -0.1; xmax = 0.15;
		fclose(setupfile);
	}

	fclose(setupfile);

	// Read in task from file and scale
	DefineTask();

	// Open file to log state data
	logfile.open(logpath);
	logfile << "time, x, y, z, xball, yball, zball, flags,fx,fy,fz, scoring_enabled, ball_fx, ball_fy\n";

	if (mode == 0) // if using HapticMASTER
	{

		// Open the HapticMaster device
		dev = haOpenDevice(IPADDRESS);

		// Check if device connected, then initialize
		if (dev == HARET_ERROR) {
			printf("--- ERROR: Unable to connect to device: %s\n", IPADDRESS);
			//return HARET_ERROR; // comment this out to test visualization without the HapticMaster
		}
		else {
			InitializeDevice(dev);
		}

		haSendCommand(dev, "set inertia 4.0", response);

		// Fix arm in place using springs 
		printf("Moving...\n");
		returnValue = haSendString(dev, "create spring mySpring", response);
		haSendCommand(dev, "set mySpring pos", springPos[PosX], springPos[PosY], springPos[PosZ] + 0.011, response);
		haSendCommand(dev, "set mySpring stiffness", springStiffness, response);
		haSendCommand(dev, "set mySpring dampfactor", springDamping, response);
		haSendCommand(dev, "set mySpring maxforce", springMaxForce, response);
		returnValue = haSendString(dev, "set mySpring enable", response);

		// Set bias force
		printf("Setting bias force to %f of the max abduction force.\n", support_level[support_num]);
		float bias = maxForce * support_level[support_num]; //+ 20;
		returnValue = haSendString(dev, "create biasforce myBiasForce", response);
		printf("bias: %f \n", bias);
		returnValue = haSendCommand(dev, "set myBiasForce force", 0.0, 0.0, bias, response);
		printf("Bias force: %s\n", response);

		// Wait to start visualization/trial
		printf("Press ENTER once arm in position.\n");
		cin.get();

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
		sys.Xcurr = { 0.0,0.0,0.0,0.0,CurrentPosition[PosX],0.0,CurrentPosition[PosY],0.0 };
	}
	else
	{
		// Set initial conditions
		CurrentPosition[PosX] = 0; CurrentPosition[PosY] = 0; CurrentPosition[PosZ] = 0;
		sys.Xcurr = { 0.0,0.0,0.0,0.0,CurrentPosition[PosX],0.0,CurrentPosition[PosY],0.0 };
	}

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

	renderer = new Renderer;
	renderer->init();

	if (mode == 0) // if using HapticMASTER
	{
		// Start Ball Feedback
		returnValue = haSendString(dev, "create biasforce BallFeedbackForce", response);
		returnValue = haSendCommand(dev, "set BallFeedbackForce force [0.0, 0.0, 0.0]", response);
		haSendString(dev, "set BallFeedbackForce enable", response);
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

	glutTimerFunc(10, TimerCB, 1);
	
	glutMainLoop();

	glfwTerminate( );
	
	return 0;

}