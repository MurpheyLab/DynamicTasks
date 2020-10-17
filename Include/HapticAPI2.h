//*****************************************************************
// Copyright © 2010 Moog-FCS B.V.
//*****************************************************************
//
// file    : HapticAPI2.h
// purpose : See *.cpp file ( HapticAPI C interface )
// history : See *.cpp file
// ----------------------------------------------------------------

#ifndef HAPTICAPI_H_
#define HAPTICAPI_H_

#include "HapticAPI2Def.h"

#include <stdio.h>

#ifdef WIN32
   #ifdef HAPTICAPI_EXPORTS
      #define HAPTIC_API __declspec(dllexport) 
   #else
      #define HAPTIC_API __declspec(dllimport)
   #endif
#else
   #define HAPTIC_API
#endif

typedef struct {
   char address [32];
   char device  [32];
   char name    [32];
   char version [32];
   char date    [32];
   char os      [32];
} deviceInfo_t;

// HapticAPI
extern "C" HAPTIC_API int  haGetDevices               ( deviceInfo_t inDeviceInfo[], int inSize );
extern "C" HAPTIC_API int  haBeginBuffer              ( long inDev );
extern "C" HAPTIC_API int  haEndBuffer                ( long inDev );

// Haptic Device
extern "C" HAPTIC_API long haOpenDevice               ( const char* inAddress );
extern "C" HAPTIC_API int  haCloseDevice              ( long inDev );

extern "C" HAPTIC_API int  haSendString               ( long inDev, const char* inCommand, char* outCommand );

extern "C" HAPTIC_API int  haSetState                 ( long inDev, const char* inState );
extern "C" HAPTIC_API int  haGetState                 ( long inDev, char* outState );

extern "C" HAPTIC_API int  haSetParameter             ( long inDev, const char* inParameter, const double inValue );
extern "C" HAPTIC_API int  haGetParameter             ( long inDev, const char* inParameter, double* outValue );
extern "C" HAPTIC_API int  haSetParameter3            ( long inDev, const char* inParameter, const double inValue[3] );
extern "C" HAPTIC_API int  haGetParameter3            ( long inDev, const char* inParameter, double outValue[3] );
extern "C" HAPTIC_API int  haSetParameter4            ( long inDev, const char* inParameter, const double inValue[4] );
extern "C" HAPTIC_API int  haGetParameter4            ( long inDev, const char* inParameter, double outValue[4] );

extern "C" HAPTIC_API int  haCalibrateForceSensor     ( long inDev );

extern "C" HAPTIC_API int  haCreateEffect             ( long inDev, const char* inType, const char* inName );
extern "C" HAPTIC_API int  haRemoveEffect             ( long inDev, const char* inType, const char* inName );
extern "C" HAPTIC_API int  haRemoveAllEffects         ( long inDev );

extern "C" HAPTIC_API int  haSetEffectParameter       ( long inDev, const char* inName, const char* inParameter, const double inValue );
extern "C" HAPTIC_API int  haGetEffectParameter       ( long inDev, const char* inName, const char* inParameter, double* outValue );
extern "C" HAPTIC_API int  haSetEffectParameter3      ( long inDev, const char* inName, const char* inParameter, const double inValue[3] );
extern "C" HAPTIC_API int  haGetEffectParameter3      ( long inDev, const char* inName, const char* inParameter, double outValue[3] );
extern "C" HAPTIC_API int  haEnableEffect             ( long inDev, const char* inName );
extern "C" HAPTIC_API int  haDisableEffect            ( long inDev, const char* inName );
extern "C" HAPTIC_API int  haIsEffectEnabled          ( long inDev, const char* inName, bool* outEnabled );

// Realtime Scope
typedef double** matrix;
extern "C" HAPTIC_API long haOpenScope                ( const char* inAddress );
extern "C" HAPTIC_API int  haCloseScope               ( long inScope );

extern "C" HAPTIC_API int  haScopeConfigure           ( long inScope, const int inFlushSize = 128, const int inSubSampleRatio = 1);
extern "C" HAPTIC_API void haScopeCreateMatrix        ( const int inFlushSize, const int inSamplesCount, matrix& outMatrix );
extern "C" HAPTIC_API void haScopeDeleteMatrix        ( const int inFlushSize, const int inSamplesCount, matrix& outMatrix );

extern "C" HAPTIC_API int  haAddScopeParameter        ( long inScope, const char* inParam, int& outSampleCount );
extern "C" HAPTIC_API int  haRemoveScopeParameter     ( long inScope, const char* inParam, int& outSampleCount );
extern "C" HAPTIC_API int  haRemoveAllScopeParameters ( long inScope );

extern "C" HAPTIC_API int  haStartScope               ( long inScope );
extern "C" HAPTIC_API int  haStopScope                ( long inScope );
extern "C" HAPTIC_API int  haGetScopeStatus           ( long inScope, bool& outScoping, int& outBufferSize );

extern "C" HAPTIC_API int  haScopeFlushMatrix         ( long inScope, const int inFlushSize, const int inSamplesCount, matrix& outMatrix );
extern "C" HAPTIC_API int  haScopeFlushVector         ( long inScope, const int inFlushSize, const int inSamplesCount, double outData[] );

#endif

