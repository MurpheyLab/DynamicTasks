//---------------------------------------------------------------------
// parameters.hpp
// This file is used to set adjustable parameters including the number
// of flags and the subject-speficic parameters (subject number, trial
// number, loading level, etc.). *This should be the only file changed
// while running experiments.*
//---------------------------------------------------------------------

// Trial parameters (update every iteration)
int trial_num = 1;
int support_num = 0; // 0: 0%, 1: 30%
int freq_num = 1; // 0: 0.5Hz 1: 1Hz 2: 1.5Hz 3:2Hz 4:2.5Hz
int feedback_forces = 1; // 0: off, 1: on

// Subject parameters (update once at the beginning)
int subject_num = 3;

// Input mode
const int mode = 2; // 0: HapticMASTER control 1: joystick control 2: no input

// Visualization mode
#define TEXTURES_ON // uncomment this line to use textures

// Game parameters
#define TF 30.0 // length of trial
double percent_height = 0.3; // Allowed ball height
char game_version = 'i'; // 3 different game versions: 'f' for reading flag coords from csv file, 'i' for infinite flags, and 'p' for penalty for dropping the ball
// NOTE: version 'p' is not currently supported
int task_num = 0; // 0 is an example task, 1-5 are actual tasks

// Number of flags for task
const int num_csv_flags = 20; //change according to the number of flags in csv file
const int num_rand_flags = 5; //change according to the number of random flags to appear on table

// Define possible frequency R/damping combinations -- corresponds to 0.5,1,1.5,2,2.5Hz
const int num_freqs_tested = 4;
double radius_options[5] = { 0.995, 0.249, 0.111, 0.062, 0.04 };
double damping_options[5] = { 0.23838, 0.0149, 0.00294, 0.00093, 0.00038 }; // 10s settling
//float damping_options[5] = { 0.15892, 0.00993, 0.00196, 0.00062, 0.00025 }; // 15s settling
double R = radius_options[freq_num];
double damping = damping_options[freq_num];

// Define possible support levels
float support_level[2] = { -0.0, -0.3 }; // fraction of max shoulder abduction loading (0-1)

//float damping = 0.0;
// Difficulty/frequency parameters (R represents pendulum length)
// original (1.88 Hz)
//float R = 0.07;
//float damping = 0.00313; // 0.005;

// 0.5 Hz
//int freq = 0;
//float R = 0.995;
//float damping = 0.15892; // 0.23838;
// 1 Hz
//int freq = 1;
//float R = 0.249;
//float damping = 0.00993; // 0.0149;
// 1.5 Hz
//int freq = 4;
//float R = 0.111;
//float damping = 0.00196; // 0.00294;
// 2 Hz
//int freq = 2;
//float R = 0.062;
//float damping = 0.00062; // 0.00093;
// 2.5 Hz
//int freq = 5;
//float R = 0.04;
//float damping = 0.00025; // 0.00038;
// 3 Hz
//int freq = 3;
//float R = 0.028;
//float damping = 0.00012; // 0.00018;
// 4 Hz
//float R = 0.016;
//float damping = 0.00004; // 0.00006;

