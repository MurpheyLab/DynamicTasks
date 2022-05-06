//---------------------------------------------------------------------
// taskpaths.hpp
// This file contains paths for the setup file, task files, and texture
// files. These paths should be updated when switching to a new computer
// or reorganizing the repo.
//---------------------------------------------------------------------

#include "parameters.hpp"
using namespace std;

// File paths
const char* taskpath1 = "C:\\Users\\numur\\Documents\\DynamicTasks\\Tasks20\\spiral_task2.csv";
const char* taskpath2 = "C:\\Users\\numur\\Documents\\DynamicTasks\\Tasks20\\bowtie_task2_20.csv";
const char* taskpath3 = "C:\\Users\\numur\\Documents\\DynamicTasks\\Tasks20\\clusters_task1.csv";
const char* taskpath4 = "C:\\Users\\numur\\Documents\\DynamicTasks\\Tasks20\\clusters_task2.csv";
const char* taskpath5 = "C:\\Users\\numur\\Documents\\DynamicTasks\\Tasks20\\clusters_task5.csv";
const char* taskpath0 = "C:\\Users\\numur\\Documents\\DynamicTasks\\Tasks20\\ex_task.csv";

string flagspathstring = "C:\\Users\\numur\\Documents\\DynamicTasks\\WorkspaceCSVs\\flags_example.csv";
setuppathstring = "C:\\Users\\numur\\Documents\\DynamicTasks\\WorkspaceCSVs\\setup_example.csv";
// string flagspathstring = "C:\\Users\\numur\\Documents\\DynamicTasks\\SubjectData\\Flags_S" + to_string(subject_num) + "_SL2_A" + to_string(arm) + ".csv";
// string setuppathstring = "C:\\Users\\numur\\Documents\DynamicTasks\\SubjectData\\setup_S" + to_string(subject_num) + "_A" + to_string(arm) + ".csv";
const char* setuppath = setuppathstring.c_str();
const char* flagspath = flagspathstring.c_str();

//const char* tablepath = "C:\\Users\\numur\\Documents\\DynamicTasks\\Textures\\wood_table.png"; // old table
const char* tablepath = "C:\\Users\\numur\\Documents\\DynamicTasks\\Textures\\wood_2.png";
const char* flagpath = "C:\\Users\\numur\\Documents\\DynamicTasks\\Textures\\mushroom.png";
const char* easyLpath = "C:\\Users\\numur\\Documents\\DynamicTasks\\Textures\\easy_L.png";
const char* easyRpath = "C:\\Users\\numur\\Documents\\DynamicTasks\\Textures\\easy_R.png";
const char* mediumLpath = "C:\\Users\\numur\\Documents\\DynamicTasks\\Textures\\medium_L.png";
const char* mediumRpath = "C:\\Users\\numur\\Documents\\DynamicTasks\\Textures\\medium_R.png";
const char* intenseLpath = "C:\\Users\\numur\\Documents\\DynamicTasks\\Textures\\intense_L.png";
const char* intenseRpath = "C:\\Users\\numur\\Documents\\DynamicTasks\\Textures\\intense_R.png";

string dirpath = "C:\\Users\\numur\\Documents\\DynamicTasks\\SubjectData\\";
// string personalbestspathstring = "C:\\Users\\numur\\Documents\\DynamicTasks\\SubjectData\\personalbests_A" + to_string(arm) + ".csv";
string personalbestspathstring = "C:\\Users\\numur\\Documents\\DynamicTasks\\WorkspaceCSVs\\personalbests_example.csv";
const char* personalbestspath = personalbestspathstring.c_str();
const char* trophypath = "C:\\Users\\numur\\Documents\\DynamicTasks\\Textures\\personal-best-trophy2.png"; // trophy.png";
const char* liftpath = "C:\\Users\\numur\\Documents\\DynamicTasks\\Textures\\liftred.png";
