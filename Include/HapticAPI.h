//*****************************************************************
// Copyright © 2010 Moog-FCS B.V.
//*****************************************************************
//
// file    : HapticAPI.h
// purpose : See *.cpp file ( HapticAPI C interface )
// history : See *.cpp file
// ----------------------------------------------------------------

#ifndef HAPTICAPI_H_
#define HAPTICAPI_H_

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
extern "C" HAPTIC_API int  haDevicesGetList           ( deviceInfo_t inDeviceInfo[], int inSize );

#define HARET_SUCCESS 0
#define HARET_ERROR -1

// Haptic Device
extern "C" HAPTIC_API long haDeviceOpen               ( const char* inAddress );
extern "C" HAPTIC_API int  haDeviceClose              ( long inDev );
extern "C" HAPTIC_API int  haDeviceSendString         ( long inDev, const char* inCommand, char* outCommand );
extern "C" HAPTIC_API int  haGetAPIDllVersion         ( char* outVersion );


// Realtime DataLogger
typedef float** matrix;
#define DEF_FLUSHSIZE_MAX     2048
#define DEF_SUBSAMPLE_RATIO   1

extern "C" HAPTIC_API long haDataLoggerOpen           ( const char* inAddress );
extern "C" HAPTIC_API int  haDataLoggerClose          ( long inDataLogger );

extern "C" HAPTIC_API int  haDataLoggerConfigure      ( long inDataLogger, const int inFlushSizeMax, const int inSubSampleRatio );
extern "C" HAPTIC_API void haDataLoggerAllocMatrix    ( const int inFlushSizeMax, const int inColumnCount, matrix& ioMatrix );
extern "C" HAPTIC_API void haDataLoggerFreeMatrix     ( const int inColumnCount, matrix& ioMatrix );

extern "C" HAPTIC_API int  haDataLoggerAddParameter   ( long inDataLogger, const char* inParameter, int* outColumnCount );
extern "C" HAPTIC_API int  haDataLoggerRemoveParameter( long inDataLogger, const char* inParameter, int* outColumnCount );
extern "C" HAPTIC_API int  haDataLoggerRemoveAllParameters ( long inDataLogger );

extern "C" HAPTIC_API int  haDataLoggerStart          ( long inDataLogger );
extern "C" HAPTIC_API int  haDataLoggerStop           ( long inDataLogger );
extern "C" HAPTIC_API int  haDataLoggerGetStatus      ( long inDataLogger, bool* outLogging, int* outBufferSize, int* outFlushSizeMax, int* outColumnCount );

extern "C" HAPTIC_API int  haDataLoggerFlushMatrix    ( long inDataLogger, matrix& outMatrix );
extern "C" HAPTIC_API int  haDataLoggerFlushVector    ( long inDataLogger, float* outData );

#ifdef WIN32
extern "C" HAPTIC_API int  haDataLoggerStartThread    ( long inDataLogger, const char* inFileName, long& outProcessHandle );
extern "C" HAPTIC_API int  haDataLoggerStopThread     ( long inDataLogger, long inProcessHandle );
#endif

#endif

