//---------------------------------------------------------------------
// parameters.hpp
// This file is used to set adjustable parameters including the number
// of flags and the subject-speficic parameters (subject number, trial
// number, loading level, etc.). *This should be the only file changed
// while running experiments.*
//---------------------------------------------------------------------

// Subject parameters (update once at the beginning)
int subject_num = 103;

// Trial parameters (update every iteration)
int arm = 0; // 0: paretic; 1: non-paretic
int trial_num = 80;
int support_num = 2; // 0: haptic table (always), 1: 0%, 2: 30%, 3: 50%
int freq_num = 0; // 0: 0.5Hz 1: 1Hz 2: 1.5Hz 3: 2.5Hz
// 0: 0.25Hz 1: 0.5Hz 2: 0.75Hz 3: 1Hz 4: 1.5Hz 5:2Hz 6:2.5Hz
int personalbests_flag = 0; // 0: off for learning, 1: on for data collection

int condition = 0;
const int feedback_forces_vec[4] = { 1, 1, 0, 0 }; // 0: off, 1: on
int ball_moving_vec[4] = { 1, 0, 1, 0 }; // 0: not moving, 1: moving 
const int feedback_forces = feedback_forces_vec[condition];
int ball_moving = ball_moving_vec[condition];
	
// Input mode
const int mode = 0; // 0: HapticMASTER control 1: joystick control 2: no input

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
const int num_rand_flags = 3; //change according to the number of random flags to appear on table
float angle_disturbance = 3.14 / 3.5;

// Define possible frequency R/damping combinations -- corresponds to 0.5,1,1.5,2,2.5Hz
const int num_freqs_tested = 4;
//double radius_options[7] = { 3.98, 0.995, 0.442, 0.249, 0.111, 0.062, 0.04 };
double radius_options[4] = { 0.995, 0.249, 0.111, 0.04 };

//double damping_options[5] = { 0.23838, 0.0149, 0.00294, 0.00093, 0.00038 }; // 10s settling
//double damping_options[5] = { 0.15892, 0.00993, 0.00196, 0.00062, 0.00025 }; // 15s settling
//double damping_options[5] = { 0.15892/2, 0.00993/2, 0.00196/2, 0.00062/2, 0.00025/2 }; // 15s settling
//double damping_options[5] = { 0.137238/2, 0.008577/2, 0.001694/2, 0.000536/2, 0.000220/2 }; // 10s settling to 30 degrees
//double damping_options[5] = { 0.03531, 0.002207, 0.000436, 0.000138, 0.000056 }; // 15s settling

//double damping_options[7] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
//double damping_options[4] = { 0.0, 0.0, 0.0, 0.0 };

double R = radius_options[freq_num];
double damping = 0.0; // damping_options[freq_num];

// Define possible support levels
float support_level[5] = { 0.5, 0.0, 0.3, 0.5, 0.5 }; // fraction of max shoulder abduction loading (0-1)
const int max_support = 2; // support level index that represents the most load

//// Difficulty/frequency parameters (R represents pendulum length)
//// original (1.88 Hz)
//float R = 0.07;
//float damping = 0.005;