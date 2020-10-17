//*****************************************************************
// Copyright © 2010 Moog-FCS B.V.
//*****************************************************************
//
// Name        : HapticAPI2Def.h
// Class       : 
// Description : See *.cpp file
// History     : See *.cpp file
// ----------------------------------------------------------------
/*
   HapticAPI Definitions
*/

/* haptic device states */
#define HA_STATE_OFF                   "off"                // not powered
#define HA_STATE_INIT                  "init"               // initialing... (calibrating force sensor, searching endstop...)
#define HA_STATE_STOP                  "stop"               // device is powered and not rendering forces and not responding to the force sensor
#define HA_STATE_POSITION              "position"           // device is powered and is rendering forces, but not responding to the force sensor
#define HA_STATE_FORCE                 "force"              // force sensitive and applying newtons law

/* haptic device/scope parameters */
#define HA_INERTIA                     "inertia"            // [kg]
#define HA_MODEL_POS                   "modelpos"           // [m]
#define HA_MODEL_VEL                   "modelvel"           // [m/s]
#define HA_MODEL_ACC                   "modelacc"           // [m/s^2]
#define HA_MEAS_FORCE                  "measforce"          // [N]
#define HA_MEAS_ORIENT                 "measorient"        	// [rad]
#define HA_MEAS_POS                    "measpos"            // [m]
#define HA_MEAS_VEL                    "measvel"            // [m/s]
#define HA_EXTL_FORCE                  "extlforce"          // [N]
#define HA_MODEL_POS_JOINT             "modelposjoint"      // [m,rad,m]
#define HA_MODEL_VEL_JOINT             "modelveljoint"      // [m/s,rad/s,m/s]
#define HA_MODEL_ACC_JOINT             "modelaccjoint"      // [m/s^2,rad/s^2,m/s^2]
#define HA_MEAS_FORCE_JOINT            "measforcejoint"     // [N]
#define HA_MEAS_POS_JOINT              "measposjoint"       // [m,rad,m]

// Haptic effects supported by the renderer, for use in
// the haCreateEffect function. After a successfull 
// creation, the effect parameters can be set using 
// haSetEffect calls.
#define HA_EFFECT_SPRING               "spring"
#define HA_EFFECT_DAMPER               "damper"
#define HA_EFFECT_SHAKER               "shaker"
#define HA_EFFECT_SPHERE               "sphere"
#define HA_EFFECT_BLOCK                "block"
#define HA_EFFECT_BIASFORCE            "biasforce"

/* generic effects properties */
#define HA_EFFECT_POS                  "pos"               // [m]
#define HA_EFFECT_VEL                  "vel"               // [m/s]
#define HA_EFFECT_ACC                  "acc"               // [m/s^2]

/* parameters for the spring effect */
#define HA_SPRING_STIFFNESS            "stiffness"         // [Nm]
#define HA_SPRING_DAMPFACTOR           "dampfactor"        // [-]
#define HA_SPRING_DEADBAND             "deadband"          // [m]
#define HA_SPRING_DIRECTION            "direction"
#define HA_SPRING_MAXFORCE             "maxforce"          // [N]
#define HA_SPRING_DAMPGLOBAL           "dampglobal"

/* parameters for the damper effect */
#define HA_DAMPER_DAMPCOEF             "dampcoef"          // [N s/m]

/* parameters for the sphere effect */
#define HA_SPHERE_STIFFNESS            "stiffness"         // [Nm]
#define HA_SPHERE_DAMPFACTOR           "dampfactor"       
#define HA_SPHERE_RADIUS               "radius"            

/* parameters for the block shape */
#define HA_BLOCK_STIFFNESS             "stiffness"          // [Nm]
#define HA_BLOCK_DAMPFACTOR            "dampfactor"
#define HA_BLOCK_SIZE                  "size"               // [m]

/* parameters for the bias force effect */
#define HA_BIASFORCE_FORCE             "force"              // [N]


/* parameter for the shaker effect */
#define HA_SHAKER_FREQUENCY1           "frequency1"         // [hz]
#define HA_SHAKER_FREQUENCY2           "frequency2"         // [hz]
#define HA_SHAKER_DIRECTION            "direction"
#define HA_SHAKER_POSMAX               "posmax"             // [m]
#define HA_SHAKER_VELMAX               "velmax"             // [m/s]
#define HA_SHAKER_ACCMAX               "accmax"             // [m/s^2]
#define HA_SHAKER_STIFFNESS            "stiffness"          // [N/m]
#define HA_SHAKER_DAMPFACTOR           "dampfactor"         // [N s/m]
#define HA_SHAKER_DEADBAND             "deadband"           // [m]
#define HA_SHAKER_MAXFORCE             "maxforce"           // [N]

/* extra parameters supported by the scope only */
#define HA_SCOPE_SAMPLENR              "samplenr"
#define HA_SCOPE_TIMESTAMP             "sampletime"

/* bubblesheet */
#define HA_BS_TOOLTYPE_DRILLSPHERE     0
#define HA_BS_TOOLTYPE_DRILLCYLINDER   1

/* grid3d */
#define HA_GRID3D_TOOLTYPE_ELECTRICAL   0
#define HA_GRID3D_TOOLTYPE_AIRROTOR     1
#define HA_GRID3D_TOOLSHAPE_CYLINDRICAL 0
#define HA_GRID3D_TOOLSHAPE_SPHERICAL   1

/* types and return codes */
typedef char HASTATE [32];
typedef char HAPRM   [32];
typedef char HATYPE  [32];

#define HA_STATE                       "state"
#define HA_ENABLE                      "enable"
#define HA_MODE                        "mode"

#define HARET_ERROR    -1  // we got an error, use the haGetError for the human readable version
#define HARET_SUCCESS   0  // call ended successfully
#define HARET_BUFFERED  1  // call ended successfully, message is buffered

