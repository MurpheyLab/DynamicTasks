//---------------------------------------------------------------------
// ballinbowl_joystick.cpp
// Interactive bowl visualization moved using a PlayStation 4 Controller
//---------------------------------------------------------------------
#define _USE_MATH_DEFINES
#include "HapticAPI.h"
#include "HapticUtility.h"
#include <glfw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <glut.h>
#include <math.h> // sine, cosine
#include "cmath"
#include <cstdlib>
#include <iostream> // cin
#include <fstream>
#include <sstream>
#include <ctime>
#include "paths.hpp"
#include "quaternions.hpp"
using namespace std;

//Screen dimension constants
const int SCREEN_WIDTH = 2500;
const int SCREEN_HEIGHT = 1200;

//---------------------------------------------------------------------
//    Variables for adjustment during experiments
//---------------------------------------------------------------------

// at the beginning only
int subject_num = 3;
float maxForce = 71.465; //input value from isometric protocol (Fx)
#define TF 20.0 // length of trial

// define game version
const int mode = 0; // '0' connected to robot, '1' joystick
char game_version = 'i'; // 2 different game versions: 'a' for nails on an arc, 'i' for infinite nails (two visible at a time)
const int nail_forces = 1; // option to enable wall feedback: '0' if off, '1' if on
const int rightleft = 0; // '0' for nails to appear on the right and '1' for nails to appear on the left and '2' for nails to appear all around

int trial_num = 15;
int support_num = 0; // 0: 0%, 1: 30%

//---------------------------------------------------------------------
//---------------------------------------------------------------------

// Define possible support levels
float support_level[3] = { -0.0, -0.3 }; // fraction of max shoulder abduction loading (0-1)

// Files paths for defining task, reading in arm weight and workspace, and logging data
//string logpath = dirpath + "NailHammerData_S" + to_string(subject_num) + "_Trial" + to_string(trial_num) + "_Task" + to_string(rightleft) + "_SL" + to_string(support_num) + ".csv";
string logpath = "NailHammerData_S" + to_string(subject_num) + "_Trial" + to_string(trial_num) + "_Task" + to_string(rightleft) + "_SL" + to_string(support_num) + ".csv";
ofstream logfile;

#define IPADDRESS "10.30.203.26" //"10.30.203.26"

#define PosX 0
#define PosY 1
#define PosZ 2

//#define deltaT 0.01
float sys_time = 0.0;
float beg_sys_time = 0.0;
float deltaT = 0.01;
#define table_z -0.14 // -0.17
double z_tolerance = table_z + 0.01; // -0.16 // if lower, person can eat flags while slacking
float size_of_section = 0.018;
float tol_delta = size_of_section * 0.6; // for "nailing in" a section
#define angle_tol 20 // 15 // alignment of movement direction with angle of nail
#define minvelx 0.1
#define minvely 0.1

//GLfloat xmin = -0.08, xmax = 0.15, ymin = -0.2, ymax = 0.25;
GLfloat xmin = 10.0; GLfloat xmax = 10.0; GLfloat ymin = 10.0; GLfloat ymax = 10.0;
GLfloat armWeight = 10.0;

long dev = 0;
char response[100];
//char response[3] = { '0','0','0' };
double CurrentPosition[3];
double CurrentVisPosition[3];
double CurrentJoystick[2];
double CurrentForce[3];

// initial starting position
double springPos[3] = { -0.05,0.05,table_z }; // x was 0.0

// variables for nail task
#define numsections 7
#define numnails 9
const int numgoals = numnails * numsections; // 7 nails x 5 sections each
GLfloat goals[numnails][numsections][5]; // x, y, z, angle, color, nail#
GLfloat wallnails[numnails][numsections][5];
int goal_status[numnails][numsections][2];

int currentnail = 0;
int current_block = (int) (numsections - 1);
int anglevector[numnails];
float quat[4] = { 0.0, 0.0, 0.0, 1.0 };
float currentangle;
double velx, vely;
int score = 0;
float trial_flag = 0.0;
float prevtime = 0.0;
float scoring_enabled = 0.0;
double xprev[3] = {0.0,0.0,0.0}; 
double yprev[3] = {0.0,0.0,0.0};

// timer
bool resume_flag = true;
float penalty_timer;
bool timer_initialized = false;
float penalty_duration;
#define time_of_inactivity 0.4

// define arc
float arcrad = 0.29;
GLfloat XCenter = 0.1; GLfloat YCenter = 0.0;
#define cursorsize 0.02

//movement scaling
float maxXrange = arcrad + 0.05;
float maxYrange = 2*arcrad + 0.05;
float actualXrange, actualYrange; // read in from a file in main
float workspaceXstart, workspaceYstart;

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
//               D R A W   T A B L E   &   A R C
//
// This Function Is Called To Draw A Plane
//---------------------------------------------------------------------
void DrawTable(void) {
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

	// draw arc on table
	float wid = size_of_section * 0.9;
	float radius1 = arcrad - wid / 2;	float radius2 = arcrad + wid / 2;
	int startang, endang;
	if (rightleft == 1) {
		startang = 180;	endang = 270;
	} else if (rightleft == 0) {
		startang = 90; endang = 180;
	} else {
		startang = 90; endang = 270;
	}

	//glColor3f(1.0f, 1.1f, 0.0f); // yellow
	if (scoring_enabled == 1.0) {
		glColor3f(0.01f, 0.5f, 1.0f); // blue
	} else {
		glColor3f(1.0f, 1.0f, 1.0f); // white
	}
	glBegin(GL_TRIANGLE_STRIP);
	for (int i = startang; i <= endang; i += 10)
	{
		GLfloat ang = (GLfloat)(i / 180.0 * M_PI);
		glVertex3f((float)(XCenter + cos(ang) * (radius1)), (float)(YCenter + sin(ang) * (radius1)), 0.005);
		glVertex3f((float)(XCenter + cos(ang) * (radius2)), (float)(YCenter + sin(ang) * (radius2)), 0.005);
	}
	glEnd();
}

//---------------------------------------------------------------------
//                     D R A W   N A I L S
//
// This Function Is Called To Draw Nails
//---------------------------------------------------------------------
void DrawNails(void)
{
	int draw_flag = 0;

	if (game_version == 'a') {
		// draw all nails at once
		for (int n = 0; n < numnails; n += 1) {
			for (int nn = 0; nn < numnails; nn += 1) {
				if (goal_status[n][nn][0] == 1) {
					glColor3f(0.0f, 0.0f, 0.0f); // black
				}
				else if (goals[n][nn][4] == 0) {
					glColor3f(1.0f, 0.0f, 1.0f); // purple
				}
				else {
					glColor3f(1.0f, 0.5f, 0.0f); // orange
				}
				glPushMatrix();
				glTranslatef(goals[n][nn][PosX], goals[n][nn][PosY], goals[n][nn][PosZ]);
				glRotatef(goals[n][nn][3], 0.0, 0.0, 1.0);
				glutSolidCube(size_of_section); // change this to make the cubes bigger or smaller
				glPopMatrix();
			}
		}
	}
	else if (game_version == 'i') {
		for (int n = 0; n < numgoals; n += 1) {
			// draw one nail at a time
			for (int nn = 0; nn < numsections; nn += 1) {
				if ((goals[n][nn][4] == 0) && (n == currentnail) && (goal_status[n][nn][0] == 0)) {
					glColor3f(1.0f, 0.0f, 1.0f); // purple
					draw_flag = 1;
				}
				else if ((goals[n][nn][4] == 1) && (n == currentnail) && (goal_status[n][nn][0] == 0)) {
					glColor3f(1.0f, 0.5f, 0.0f); // orange
					draw_flag = 1;
				}
				else if ((goals[n][nn][4] == 0) && (n == currentnail + 1)) {
					glColor4f(1.0f, 0.0f, 1.0f, 0.2f); // transparent purple
					draw_flag = 1;
				}
				else if ((goals[n][nn][4] == 1) && (n == currentnail + 1)) {
					glColor4f(1.0f, 0.5f, 0.0f, 0.2f); // transparent orange
					draw_flag = 1;
				}
				else {
					glColor4f(0.0f, 0.0f, 0.0f, 0.0f); // black
					draw_flag = 0;
				}

				if (draw_flag == 1)
				{
					glPushMatrix();
					glTranslatef(goals[n][nn][PosX], goals[n][nn][PosY], goals[n][nn][PosZ]);
					glRotatef(goals[n][nn][3], 0.0, 0.0, 1.0);
					glutSolidCube(size_of_section); // change this to make the cubes bigger or smaller
					glPopMatrix();
				}

				// add nail inside wall
				if (wallnails[n][nn][4] == 1) {
					glColor4f(0.5f, 0.5f, 0.5f, 1.0f); // grey
					glPushMatrix();
					glTranslatef(wallnails[n][nn][PosX], wallnails[n][nn][PosY], wallnails[n][nn][PosZ]);
					glRotatef(wallnails[n][nn][3], 0.0, 0.0, 1.0);
					glutSolidCube(size_of_section); // change this to make the cubes bigger or smaller
					glPopMatrix();
				}
			}
		}
	}
}

//---------------------------------------------------------------------
//                    D E F I N E   T A S K
//
// This Function Is Called To Define the Location of Goal Flags
//---------------------------------------------------------------------
void DefineTask(void)
{
	if (rightleft == 1) {
		YCenter = YCenter + 0.1;
	} else if (rightleft == 0) {
		YCenter = YCenter - 0.1;
	}

	int anglemin, anglemax;
	if (rightleft == 1) {
		anglemin = 180; anglemax = 270;
	}
	else if (rightleft == 0) {
		anglemin = 100; anglemax = 190;
	}
	else {
		anglemin = 90; anglemax = 280;
	}

	if (game_version == 'i') {
		int ii = 0;
		for (int i = 0; i < numnails; i++) {
			int randangle = 10 * (rand() % (anglemax - anglemin) / 10 + anglemin / 10);  // 10 * (rand() % 18 + 9);
			//printf("i %d randangle %d.\n", i, randangle);
			while (ii<i) {
				//printf("diff %d.\n", abs(randangle - anglevector[ii]));
				if (abs(randangle - anglevector[ii])<1) {
					randangle = 10 * (rand() % (anglemax - anglemin) / 10 + anglemin / 10);
					ii = 0;
					//printf("same detected: ii %d randangle %d", ii, randangle);
				} else {
					ii++;
				}
			}
			ii = 0;
			anglevector[i] = randangle;
			//printf("angle %d\n", randangle);
		}
	}
	else {
		int ii = 0;
		for (int i = 90; i <= 270; i += 30) {
			anglevector[ii] = i;
			ii++;
		}
	}
	//draw "nails"
	float blockradius = arcrad;
	float blockradiusWALL = arcrad - size_of_section;
	int val = numsections - 1;
	for (int nn = 0; nn < numsections; nn += 1) {
		blockradius = blockradius - size_of_section;
		blockradiusWALL = blockradiusWALL + size_of_section;
		for (int n = 0; n < numnails; n += 1) {
			//GLfloat ang = (GLfloat)(i / 180.0 * M_PI);
			//goals[n][3] = i;
			GLfloat ang = (GLfloat)(anglevector[n % numnails] / 180.0 * M_PI);
			goals[n][nn][3] = anglevector[n % numnails];
			goals[n][nn][0] = XCenter + cos(ang) * (blockradius);
			goals[n][nn][1] = YCenter + sin(ang) * (blockradius);
			goals[n][nn][2] = 0.0;
			if ((nn % 2) == 0) {
				goals[n][nn][4] = 0;
			}
			else {
				goals[n][nn][4] = 1;
			}
			goal_status[n][nn][1] = val;
			goal_status[n][nn][0] = 0;
			//printf("goals:     x %f, y %f, angle %f\n", goals[n][nn][0], goals[n][nn][1], goals[n][nn][3]);

			wallnails[n][nn][0] = XCenter + cos(ang) * (blockradiusWALL);
			wallnails[n][nn][1] = YCenter + sin(ang) * (blockradiusWALL);
			wallnails[n][nn][2] = 0.0;
			wallnails[n][nn][3] = anglevector[n % numnails];
			wallnails[n][nn][4] = 0;
			//printf("wallnails: x %f, y %f, angle %f\n", wallnails[n][nn][0], wallnails[n][nn][1], wallnails[n][nn][3]);
		}
		val = val - 1;
	}

}

//---------------------------------------------------------------------
//                     C H E C K   N A I L S
//
// This Function Is Called By TimerCB to check if a flag has been
// reached
//---------------------------------------------------------------------
void CheckNails(void)
{
    GLfloat rand_x, rand_y;

	// here we check for: 1. ball is in bowl, 2. trial is running, 3. person is above the haptic table
	//printf("Nail angle is %f.\n", goals[currentnail][0][3]);
	if ((trial_flag == 1.0) && (CurrentPosition[PosZ] > z_tolerance)) {
		scoring_enabled = 1.0;
	}
	else {
		scoring_enabled = 0.0;
	}

	// game version determines how flags are scored
	if (game_version == 'a') {
		// UPDATE THIS FOR 'a' version to work
		////printf("current pos x%f, y%f\n", CurrentPosition[PosX], CurrentPosition[PosY]);
		//for (int i = 0; i < numgoals; i++) { 
		//	//printf("goal #%d x%f, y%f\n", i, goals[i][PosX], goals[i][PosY]);
		//	// eats flags when only x and y overlap
		//	if (resume_flag && (abs(CurrentPosition[PosX] - goals[i][PosX]) < tol_delta) && (abs(CurrentPosition[PosY] - goals[i][PosY]) < tol_delta) && (goal_status[i][0]==0) && (goal_status[i][1] == 0) && (scoring_enabled==1.0)) {
		//		goal_status[i][0] = 1; // set achieved goal to disappear (turn transparent)
		//		score++;
		//		resume_flag = false;
		//		// enable scoring of the next block in the nail 
		//		for (int nn = 0; nn < numgoals; nn = nn + 1) {
		//			if (goals[nn][5] == goals[i][5]) { // if the same nail
		//				goal_status[nn][1] = goal_status[nn][1] - 1; // enable next block
		//			}
		//		}
		//	}
		//	// penalty flag is true if timer is reset
		//	if ((resume_flag == false) && (timer_initialized == false)) {
		//		penalty_timer = clock();
		//		timer_initialized = true;
		//	}
		//	else if (resume_flag == false) {
		//		// reset timer if given number of seconds has passed
		//		penalty_duration = (clock() - penalty_timer) / (double)CLOCKS_PER_SEC;
		//		if (penalty_duration > time_of_inactivity) { // ADJUSTABLE duration of penalty forgiveness timer in seconds 
		//			timer_initialized = false;
		//			resume_flag = true;
		//		}
		//	}
		//}
	} else if (game_version == 'i') {
		//printf("current pos x%f, y%f\n", CurrentPosition[PosX], CurrentPosition[PosY]);
		//printf("current_block %d.\n", current_block);
		for (int i = 0; i < numgoals; i++) {
			//printf("expected angle %f, actual angle %f.\n", currentangle, goals[currentnail][current_block][3]);
			//printf("expected x pos %f, actual x pos %f.\n", CurrentVisPosition[PosX],goals[currentnail][current_block][PosX]);
			if ((abs(currentangle - goals[currentnail][current_block][3]) < angle_tol) && resume_flag && (abs(CurrentVisPosition[PosX] - goals[currentnail][current_block][PosX]) < tol_delta) && (abs(CurrentVisPosition[PosY] - goals[currentnail][current_block][PosY]) < tol_delta) && (scoring_enabled == 1.0)) {
				goal_status[currentnail][current_block][0] = 1; // set achieved goal to disappear (turn transparent)
				wallnails[currentnail][numsections - current_block - 1][4] = 1;
				score++;
				resume_flag = false;
				current_block = current_block - 1;
				// enable scoring of the next block in the nail 
				for (int nn = 0; nn < numsections; nn++) {
					goal_status[currentnail][nn][1] = goal_status[currentnail][nn][1] - 1; // enable next block
					// change colors to make it seem like nail is going into wall
					if (goals[currentnail][nn][4] == 0.0) { goals[currentnail][nn][4] = 1.0; }
					else { goals[currentnail][nn][4] = 0.0; }
				}
				if (goal_status[currentnail][0][1] < 0) {
					printf("RESET\n");
					currentnail++;
					printf("current nail %d\n", currentnail);
					current_block = numsections - 1;
				}
			}
			// penalty flag is true if timer is reset
			if ((resume_flag==false) && (timer_initialized==false)) {
				penalty_timer = clock();
				timer_initialized = true;
			} else if (resume_flag==false) {
				// reset timer if given number of seconds has passed
				// printf("penalty duration: %f.\n", penalty_duration);
				penalty_duration = (clock() - penalty_timer) / (float) CLOCKS_PER_SEC;
				if (penalty_duration > time_of_inactivity) { // ADJUSTABLE duration of penalty forgiveness timer in seconds 
					timer_initialized = false;
					resume_flag = true;
					if ((mode == 0) && (nail_forces == 1)) {
						haSendCommand(dev, "set myNail disable", response);
						printf("current nail: %d, current angle: %f, current block: %d.\n", currentnail, goals[currentnail][0][3], current_block);
						float wallX, wallY;
						if ((current_block - 2) >= 0) {
							wallX = (goals[currentnail][current_block - 2][0] - XCenter) * actualXrange / maxXrange + springPos[PosX];
							wallY = (goals[currentnail][current_block - 2][1] - YCenter) * actualYrange / maxYrange + springPos[PosY];
						}
						else if (current_block == 1) {
							wallX = (wallnails[currentnail][0][0] - XCenter) * actualXrange / maxXrange + springPos[PosX];
							wallY = (wallnails[currentnail][0][1] - YCenter) * actualYrange / maxYrange + springPos[PosY];
						}
						else if (current_block == 0) {
							wallX = (wallnails[currentnail][1][0] - XCenter) * actualXrange / maxXrange + springPos[PosX];
							wallY = (wallnails[currentnail][1][1] - YCenter) * actualYrange / maxYrange + springPos[PosY];
						}
						haSendCommand(dev, "set myNail pos", wallX, wallY, goals[currentnail][0][PosZ], response); //- size_of_section / 4
						//haSendCommand(dev, "set myNail pos", goals[currentnail][current_block][0], goals[currentnail][current_block][1], goals[currentnail][current_block][2], response);
						haSendCommand(dev, "set myNail att", QuaternionLookup(goals[currentnail][0][3], 0), QuaternionLookup(goals[currentnail][0][3], 1), QuaternionLookup(goals[currentnail][0][3], 2), QuaternionLookup(goals[currentnail][0][3], 3), response);
						haSendCommand(dev, "set myNail enable", response);
						printf("Nail wall enabled...\n");
					}
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

	int time_current_rounded = round(TF - sys_time);
	string time_current_string = to_string(time_current_rounded);
	char *ctime = const_cast<char*>(time_current_string.c_str());

    //glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

	// thicken the text
	for( float i = 0.0; i <= 0.002; i = i + 0.0001) {
		for (float j = 0.0; j <= 0.002; j = j + 0.0001) {
			glColor3f(0.0f, 0.7f, 1.0f);
			DrawText(-0.45 + i, 0.19 - j, "SCORE");

			if (score < 10) {
				DrawText(-0.395 + i, 0.14 - j, cscore); //DrawText(-0.385 + i, 0.14, cscore);
			}
			else {
				DrawText(-0.405 + i, 0.14 - j, cscore); //DrawText(-0.395 + i, 0.14, cscore);
			}

			glColor3f(1.0f, 0.0f, 0.0f);
			DrawText(0.355 + i, 0.19 - j, "TIME");

			if (time_current_rounded < 10) {
				DrawText(0.385 + i, 0.14 - j, ctime);
			}
			else {
				DrawText(0.375 + i, 0.14 - j, ctime);
			}
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
	GLfloat bottom = 0.07, width = 0.053, x_offset = 0.3, length = 0.6 * (1 - sys_time / TF);

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
	GLfloat bottom = 0.07, width = 0.053, x_offset = 0.3, length = 0.89 * score; // / num_rand_flags;

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
//                       D R A W   C U R S O R
//
// This Function Is Called To Draw Labels for the time and score
//---------------------------------------------------------------------
void DrawCursor(void)
{
	if (scoring_enabled == 1.0) {
		glColor3f(0.01f, 0.5f, 1.0f);
	}
	else {
		glColor3f(1.0f, 1.0f, 1.0f); // white 
		//glColor3f(0.388f, 0.675f, 0.745f); // always show blue square for taking screenshots
		//glColor3f(0.0f, 0.0f, 0.0f); //black
	}
	glPushMatrix();
	glTranslatef(CurrentVisPosition[PosX], CurrentVisPosition[PosY], 0.0);
	glutSolidSphere(cursorsize, 20, 20); 
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
	//glPushMatrix();
	//gluLookAt(0.68, 0.0, 0.4, -0.015, 0.0, 0.0, 0.0, 0.0, 1.0);
	//gluLookAt(0.0, 0.0, 1.5, 0.0, 0.0, 0.0, -1.0, 0.0, 0.0); // for top down view 
	//gluLookAt(0.5, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0); // for participant view point
	//gluLookAt(0.5, 0.0, 0.5, -0.2, 0.0, 0.0, 0.0, 0.0, 1.0); // for participant view point
	gluLookAt(0.645, 0.0, 0.4, 0.05, 0.0, 0.0, 0.0, 0.0, 1.0); // works well
	//gluLookAt(0.68, 0.0, 0.4, -0.015, 0.0, 0.0, 0.0, 0.0, 1.0);
	
	glutPostRedisplay();
	DrawTable();
	DrawNails();
	
	DrawTimerBar();

	// no score bar for infinite game version
	//if (game_version == 'a') {
	//	DrawScoreBar();
	//}

	DrawCursor();
	DrawLabels();

	//glPopMatrix();
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
	trial_flag = 1.0;
	if (mode == 0) // if using HapticMASTER
	{
		haDeviceSendString(dev, "set myBiasForce enable", response);
		haDeviceSendString(dev, "set mySpring disable", response);
	}
	//Set up haptic nail
	if ((mode == 0) && (nail_forces == 1)) {
		haSendCommand(dev, "create block myNail", response);
		float wallX, wallY;
		if ((current_block - 2) >= 0) {
			wallX = (goals[currentnail][current_block - 2][0] - XCenter) * actualXrange / maxXrange + springPos[PosX];
			wallY = (goals[currentnail][current_block - 2][1] - YCenter) * actualYrange / maxYrange + springPos[PosY];
		} else if (current_block == 1) {
			wallX = (wallnails[currentnail][0][0] - XCenter) * actualXrange / maxXrange + springPos[PosX];
			wallY = (wallnails[currentnail][0][1] - YCenter) * actualYrange / maxYrange + springPos[PosY];
		} else if (current_block == 0) {
			wallX = (wallnails[currentnail][1][0] - XCenter) * actualXrange / maxXrange + springPos[PosX];
			wallY = (wallnails[currentnail][1][1] - YCenter) * actualYrange / maxYrange + springPos[PosY];
		}
		haSendCommand(dev, "set myNail pos", wallX, wallY, goals[currentnail][0][PosZ], response); //-size_of_section/4
		//printf("this is the nail location %f, %f, %f.\n", goals[currentnail][current_block][0], goals[currentnail][current_block][1], -0.16);
		//haSendCommand(dev, "set myNail pos", -0.05, 0.05, -0.17, response);
		//printf("this is the nail location %f, %f, %f.\n", goals[4 * numnails][0], goals[4 * numnails][1], 0.0);
		haSendCommand(dev, "set myNail size", 0.045, 0.02, 0.5, response); //set to be a little larger than workspace
		//size_of_section * (numsections - 1)
		haSendCommand(dev, "set myNail stiffness", 20000.0, response);
		haSendCommand(dev, "set myNail att", QuaternionLookup(goals[currentnail][0][3], 0), QuaternionLookup(goals[currentnail][0][3], 1), QuaternionLookup(goals[currentnail][0][3], 2), QuaternionLookup(goals[currentnail][0][3], 3), response);
		//printf("quat: %f, %f, %f, %f.\n", QuaternionLookup(goals[currentnail][0][3], 0), QuaternionLookup(goals[currentnail][0][3], 1), QuaternionLookup(goals[currentnail][0][3], 2), QuaternionLookup(goals[currentnail][0][3], 3));
		//haSendCommand(dev, "set myNail att", 0.0, 0.0, 0.0, 1.0, response);
		haSendCommand(dev, "set myNail enable", response);
		printf("Nail wall enabled...\n");
	}
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

	// Start normal shutdown routine
	logfile.close();

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
		if (haDeviceClose(dev)) {
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
// New EndEffector Position and Updates The Bowl and Ball Position
//---------------------------------------------------------------------
void TimerCB(int iTimer)
{
	// Get current position
	if (mode == 0) { // if using HapticMASTER
		haSendCommand(dev, "get measpos", response);
		ParseFloatVec(response, CurrentPosition[PosX], CurrentPosition[PosY], CurrentPosition[PosZ]);
		CurrentVisPosition[PosX] = XCenter + (CurrentPosition[PosX] - springPos[PosX]) * maxXrange / actualXrange;
		CurrentVisPosition[PosY] = YCenter + (CurrentPosition[PosY] - springPos[PosY]) * maxYrange / actualYrange;
		//printf("current position is: %f, %f, %f.\n", CurrentPosition[PosX], CurrentPosition[PosY], CurrentPosition[PosZ]);
	}
	else { // if using a joystick
		int axesCount;
		const float* axes = glfwGetJoystickAxes(GLFW_JOYSTICK_1, &axesCount);
		CurrentJoystick[PosX] = axes[1];
		CurrentJoystick[PosY] = axes[0];
		CurrentPosition[PosX] = CurrentPosition[PosX] + CurrentJoystick[PosX] / 300;
		CurrentPosition[PosY] = CurrentPosition[PosY] + CurrentJoystick[PosY] / 300;
		CurrentVisPosition[PosX] = CurrentPosition[PosX];
		CurrentVisPosition[PosY] = CurrentPosition[PosY];
		CurrentPosition[PosZ] = 1; // always above table
	}

	// If the trial hasn't started
	if (trial_flag == 0.0) {
		xprev[0] = CurrentPosition[PosX]; xprev[1] = CurrentPosition[PosX]; xprev[2] = CurrentPosition[PosX];
		yprev[0] = CurrentPosition[PosY]; yprev[1] = CurrentPosition[PosY]; yprev[2] = CurrentPosition[PosY];
	}

	// Update simulation once trial has started
	if (trial_flag == 1.0) {
		prevtime = sys_time;
		sys_time = clock() / (float)CLOCKS_PER_SEC - beg_sys_time;
	}
	//printf("sys time: %f.\n", sys_time);
	deltaT = sys_time - prevtime;
	if (deltaT == 0.0) {
		deltaT = 0.03;
	}

	// Calculate current angle
	velx = (xprev[0] - CurrentPosition[PosX]) / (3 * deltaT);
	vely = (yprev[0] - CurrentPosition[PosY]) / (3 * deltaT);
	
	if (((abs(velx) < minvelx) && (abs(vely) < minvely)) || velx < 0.0) {
		velx = 0;
		vely = 0;
		currentangle = 0;
	} else {
		double total = pow(velx, 2.0) + pow(vely, 2.0);
		currentangle = asin(vely / sqrt(total)) * 180.0 / M_PI + 180.0;
	}

	if (trial_flag == 1.0) {
		// printf("Current velx %f and vely %f.\n", velx, vely);
		// printf("Current angle is %f.\n", currentangle);
	}
	// Update pos history for future angle calculations
	xprev[0] = xprev[1]; xprev[1] = xprev[2]; xprev[2] = CurrentPosition[PosX];
	yprev[0] = yprev[1]; yprev[1] = yprev[2]; yprev[2] = CurrentPosition[PosY];

	if (mode == 0) {// if using HapticMASTER
		haSendCommand(dev, "get measforce", response);
		ParseFloatVec(response, CurrentForce[PosX], CurrentForce[PosY], CurrentForce[PosZ]);
	} else {
		CurrentForce[PosX] = 0;  CurrentForce[PosY] = 0;  CurrentForce[PosZ] = 0;
	}

	CheckNails();

	// Write data to .txt file
	if (trial_flag == 1.0) {
		logfile << sys_time << "," << CurrentPosition[PosX] << "," << CurrentPosition[PosY] << "," << CurrentPosition[PosZ] << ",";
		logfile << CurrentVisPosition[PosX] << "," << CurrentVisPosition[PosY] << "," << CurrentVisPosition[PosZ] << ",";
		logfile << currentangle << "," << goals[currentnail][0][3] << ",";
		logfile << CurrentForce[PosX] << "," << CurrentForce[PosY] << "," << CurrentForce[PosZ] << "," << score << "," << scoring_enabled << "\n";
	}

	// Set time of trial 
	if (sys_time > TF) {
		logfile.close(); 
		if (mode == 0) { // if using HapticMASTER

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
			if (haDeviceClose(dev)) {
				printf("--- ERROR: closing device\n");
			}
		}
		exit(0);
	}

	// Poll for and process events
	glfwPollEvents();

	// Set The Timer For This Function Again
	glutTimerFunc(10, TimerCB, 1);
}

//---------------------------------------------------------------------
//                            M A I N
//---------------------------------------------------------------------
int main(int argc, char** argv)
{
	int returnValue;
	double springStiffness = 4000.0, springDamping = 10.0, springMaxForce = 100.0;

	srand(static_cast <unsigned> (time(0)));

	FILE* setupfile;
	#pragma warning (disable : 4996)
	setupfile = fopen(setuppath, "r");
	#pragma warning (disable : 4996)
	int val2 = fscanf(setupfile, "%f,%f,%f,%f,%f\n", &armWeight, &ymin, &ymax, &xmin, &xmax);

	// Option to hardcode values
	//armWeight = 60.0; ymin = -0.26; ymax = 0.26; xmin = -0.08; xmax = 0.15;

	if (armWeight == 10.0 || ymin == 10.0 || ymax == 10.0 || xmin == 10.0 || xmax == 10.0) {
		printf("error reading setupfile, hardcoded values used - max workspace\n");
		armWeight = 60.0; ymin = -0.25; ymax = 0.25; xmin = -0.1; xmax = 0.15;
		fclose(setupfile);
	}

	actualXrange = xmax - xmin;
	actualYrange = ymax - ymin;
	springPos[0] = xmax - 0.02;
	springPos[1] = ymin + actualYrange / 2;

	fclose(setupfile);

	// Read in task from file and scale
	DefineTask();

	// Open file to log state data
	logfile.open(logpath);
	logfile << "time, x, y, z, visx, visy, visz, mov_angle, nail_angle, fx, fy, fz, score, scoring_enabled\n";

	if (mode == 0) { // if using HapticMASTER
		// Open the HapticMaster device
		dev = haDeviceOpen(IPADDRESS);

		// Check if device connected, then initialize
		if (dev == HARET_ERROR) {
			printf("--- ERROR: Unable to connect to device: %s\n", IPADDRESS);
			//return HARET_ERROR; // comment this out to test visualization without the HapticMaster
		} else {
			InitializeDevice(dev);
		}

		haSendCommand(dev, "set inertia 4.0", response);

		// Fix arm in place using springs 
		printf("Moving...\n");
		returnValue = haDeviceSendString(dev, "create spring mySpring", response);
		haSendCommand(dev, "set mySpring pos", springPos[PosX], springPos[PosY], springPos[PosZ] + 0.011, response);
		haSendCommand(dev, "set mySpring stiffness", springStiffness, response);
		haSendCommand(dev, "set mySpring dampfactor", springDamping, response);
		haSendCommand(dev, "set mySpring maxforce", springMaxForce, response);
		returnValue = haDeviceSendString(dev, "set mySpring enable", response);

		// Set bias force
		printf("Setting bias force to %f of the max abduction force.\n", support_level[support_num]);
		float bias = maxForce * support_level[support_num] + armWeight; //+ 20;
		returnValue = haDeviceSendString(dev, "create biasforce myBiasForce", response);
		printf("bias: %f \n", bias);
		returnValue = haSendCommand(dev, "set myBiasForce force", 0.0, 0.0, bias, response);
		printf("Bias force: %s\n", response);

		// Wait to start visualization/trial
		printf("Press ENTER once arm in position.\n");
		cin.get();

		//// Disable springs when ready 
		//double pos_tol = 0.01;
		//while ((abs(CurrentPosition[PosX] - springPos[PosX]) > pos_tol) || (abs(CurrentPosition[PosY] - springPos[PosY]) > pos_tol) || (abs(CurrentPosition[PosZ] - springPos[PosZ] - 0.011) > pos_tol)) {
		//	haSendCommand(dev, "get measpos", response);
		//	ParseFloatVec(response, CurrentPosition[PosX], CurrentPosition[PosY], CurrentPosition[PosZ]);
		//}
		//haDeviceSendString(dev, "set mySpring disable", response);


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
		CurrentVisPosition[PosX] = XCenter; 
		CurrentVisPosition[PosY] = YCenter;
		CurrentVisPosition[PosZ] = CurrentPosition[PosZ];
		// CurrentPosition[PosX] = 0; CurrentPosition[PosY] = 0; CurrentPosition[PosZ] = 0; // should this be commented out?
	} else {
		// Set initial conditions
		CurrentPosition[PosX] = XCenter; CurrentPosition[PosY] = YCenter; CurrentPosition[PosZ] = 0;
		CurrentVisPosition[PosX] = XCenter; CurrentVisPosition[PosY] = YCenter; CurrentVisPosition[PosZ] = 0;
	}

	printf("Visualization starting...\n");

	// OpenGL Initialization Calls
	glutInit(&argc, argv);
	glfwInit();
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
	glutInitWindowSize(SCREEN_WIDTH, SCREEN_HEIGHT);

	// Create The OpenGlWindow and Initialize Things
	glutCreateWindow("Nail and hammer Game");
	InitOpenGl();
	glutReshapeFunc(Reshape);
	glutDisplayFunc(Display);
	glutKeyboardFunc(Keyboard);

	// Check for joystick
	if ((mode == 1) && (!glfwJoystickPresent(GLFW_JOYSTICK_1))) {
		std::cout << "ERROR: No Joystick Connected " << std::endl;
		exit(-1);
	}

	beg_sys_time = clock() / (float)CLOCKS_PER_SEC;
	sys_time = clock() / (float)CLOCKS_PER_SEC - beg_sys_time;

	glutTimerFunc(10, TimerCB, 1);

	glutMainLoop();

	glfwTerminate();

	return 0;

}