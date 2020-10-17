
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
//#define _USE_MATH_DEFINES
#define STB_IMAGE_IMPLEMENTATION

#include <glut.h> 
#include <math.h>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include "stb/stb_image.h" 
#include <string.h>
#include <cstring>
#include "paths.hpp"

//using namespace std;

//const int num_csv_flags = 20; //change according to the number of flags in csv file
//const int num_rand_flags = 5; //change according to the number of random flags to appear on table


class Renderer {

public:
  float t, drop_size; // 0.015
private:
  GLuint texID_table, texID_flag, texID_easyL, texID_easyR, texID_mediumL, texID_mediumR, texID_intenseL, texID_intenseR;

public:
  // constructor
  Renderer() : t(0.0), texID_table(0), texID_flag(0) {}

  // destructor
  ~Renderer() {
    if(texID_table !=0 || texID_flag != 0 || texID_easyL != 0 || texID_easyR != 0 || texID_mediumL != 0 || texID_mediumR != 0 || texID_intenseL != 0 || texID_intenseR != 0) {
      glDeleteTextures( 1, &texID_table);
      glDeleteTextures( 1, &texID_flag);
      glDeleteTextures( 1, &texID_easyL);
      glDeleteTextures( 1, &texID_easyR);
      glDeleteTextures( 1, &texID_mediumL);
      glDeleteTextures( 1, &texID_mediumR);
      glDeleteTextures( 1, &texID_intenseL);
      glDeleteTextures( 1, &texID_intenseR);
    }
  }

public:
  void init() {
    glEnable(GL_DEPTH_TEST);
    glEnable(GLenum(GL_POINT_SIZE));
    glEnable(GLenum(GL_POINT_SMOOTH));

    // assign ID# to each loaded texture
    texID_table = loadTexture(tablepath);
    texID_flag = loadTexture(flagpath);
    texID_easyL = loadTexture(easyLpath);
    texID_easyR = loadTexture(easyRpath);
    texID_mediumL = loadTexture(mediumLpath);
    texID_mediumR = loadTexture(mediumRpath);
    texID_intenseL = loadTexture(intenseLpath);
    texID_intenseR = loadTexture(intenseRpath);
  }

  void DrawTable() 
  {
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    // set camera
    //gluLookAt(0.3, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0); // water zoom in for testing water level
    gluLookAt(0.65, 0.0, 0.35, -0.015, 0.0, 0.0, 0.0, 0.0, 1.0);

    glEnable(GL_TEXTURE_2D);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBindTexture(GL_TEXTURE_2D, texID_table);
    glColor3f(1.0f,1.0f,0.0f);
    glBegin(GL_POLYGON);
    glTexCoord2f(0.0f,0.0f); glVertex3f(-0.85f, 0.3f, 0.0f);
    glTexCoord2f(1.0f,0.0f); glVertex3f(0.3f,0.3f, 0.0f);
    glTexCoord2f(1.0f,1.0f); glVertex3f( 0.3f,-0.3f, 0.0f);
    glTexCoord2f(0.0f,1.0f); glVertex3f( -0.85f, -0.3f, 0.0f);
    glEnd();

    glDisable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
  }

  void DrawFlags(GLfloat goals[num_csv_flags][3], int goal_status[num_csv_flags], float flag_size) 
  {
    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  	for (int i = 0; i < num_csv_flags; i++)
    {
      if (goal_status[i] == 0) 
      {
        GLfloat x = goals[i][0];
        GLfloat y = goals[i][1];
        GLfloat z = goals[i][2];

        glBindTexture(GL_TEXTURE_2D, texID_flag);
        glColor3f(1.0f,0.0f,0.0f);
        glBegin(GL_POLYGON);
        glTexCoord2f(0.0f,0.0f); glVertex3f(x + flag_size, y - flag_size, z);
        glTexCoord2f(1.0f,0.0f); glVertex3f(x + flag_size, y + flag_size, z);
        glTexCoord2f(1.0f,1.0f); glVertex3f(x - flag_size, y + flag_size, z + 2*flag_size);
        glTexCoord2f(0.0f,1.0f); glVertex3f(x - flag_size, y - flag_size, z + 2*flag_size);
        glEnd();

      }
    }
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
  }

  void DrawRandomFlags(GLfloat goals[num_rand_flags][3], float flag_size) 
  {
    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  	for (int i = 0; i < num_rand_flags; i++)
    {
        GLfloat x = goals[i][0];
        GLfloat y = goals[i][1];
        GLfloat z = goals[i][2];

        glBindTexture(GL_TEXTURE_2D, texID_flag);
        glColor3f(1.0f,0.0f,0.0f);
        glBegin(GL_POLYGON);
        glTexCoord2f(0.0f,0.0f); glVertex3f(x + flag_size, y - flag_size, z);
        glTexCoord2f(1.0f,0.0f); glVertex3f(x + flag_size, y + flag_size, z);
        glTexCoord2f(1.0f,1.0f); glVertex3f(x - flag_size, y + flag_size, z + 2*flag_size);
        glTexCoord2f(0.0f,1.0f); glVertex3f(x - flag_size, y - flag_size, z + 2*flag_size);
        glEnd();
	  }

    glDisable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
  }

   void DrawWaterDrops(double CurrentPosition[3], char direction, int intensity, double t, double bowl_rad) 
  {
    gluLookAt(0.68, 0.0, 0.4, -0.015, 0.0, 0.0, 0.0, 0.0, 1.0);

    if (direction == 'l' && intensity == 0)
    {
      glBindTexture(GL_TEXTURE_2D, texID_easyL);
      drop_size = 0.01;
    }
    else if (direction == 'r' && intensity == 0)
    {
      glBindTexture(GL_TEXTURE_2D, texID_easyR);
      drop_size = 0.01;
    }
    else if (direction == 'l' && intensity == 1)
    {
      glBindTexture(GL_TEXTURE_2D, texID_mediumL);
      drop_size = 0.015;
    }
    else if (direction == 'r' && intensity == 1)
    {
      glBindTexture(GL_TEXTURE_2D, texID_mediumR);
      drop_size = 0.015;
    }
    else if (direction == 'l' && intensity == 2)
    {
      glBindTexture(GL_TEXTURE_2D, texID_intenseL);
      drop_size = 0.025;
    }
    else if (direction == 'r' && intensity == 2)
    {
      glBindTexture(GL_TEXTURE_2D, texID_intenseR);
      drop_size = 0.025;
    }
    else
    {
      return;
    }

    GLfloat x, y, z; 
    x = CurrentPosition[0];

    if (direction == 'l')
    {
      y = CurrentPosition[1] - bowl_rad - t*(1.03)/3; // t* (intensity + 0.03) / 3;
    } 
    else
    {
      y = CurrentPosition[1] + bowl_rad + t*(1.03)/3; // t*(intensity+0.03)/3;
    }
    
    z = bowl_rad +0.005 - 9.81*t/100; //63 // - 0.025

    glMatrixMode(GL_MODELVIEW);
    glEnable(GL_TEXTURE_2D);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBegin(GL_POLYGON);
    glColor3f(0.0f,0.0f,1.0f);
    glTexCoord2f(0.0f,0.0f); glVertex3f(x, y-drop_size, z-drop_size);
    glTexCoord2f(1.0f,0.0f); glVertex3f(x, y+drop_size, z-drop_size);
    glTexCoord2f(1.0f,1.0f); glVertex3f(x, y+drop_size, z+drop_size);
    glTexCoord2f(0.0f,1.0f); glVertex3f(x, y-drop_size, z+drop_size);

    glEnd();
    glDisable(GL_TEXTURE_2D);
    glDisable(GL_BLEND);
  } 


private:
  // returns a valid textureID on success, otherwise 0
  GLuint loadTexture(const char* filename) {

    int width, height, numComponents;
    int level = 0;
    int border = 0;

    // loads image in the correct orientation
    stbi_set_flip_vertically_on_load(1);

    // load image data
    unsigned char* image = stbi_load(filename, &width, &height, &numComponents, STBI_rgb_alpha);

    if (image == NULL)
    {
      cout << "Cannot load texture" << endl;
    }
    
    // data is aligned in byte order
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    //request textureID
    GLuint textureID;
    glGenTextures(1, &textureID);

    // bind texture
    glBindTexture(GL_TEXTURE_2D, textureID);

    //define how to filter the texture (important but ignore for now)
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

    if(numComponents == 3)
    { 
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, image);
    }
    else if (numComponents == 4) 
    {
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
    }
    
    //texture colors should replace the original color values
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE); //GL_MODULATE

    // specify the 2D texture map
    glTexImage2D(GL_TEXTURE_2D, level, GL_RGBA, width, height, border, GL_RGBA, GL_UNSIGNED_BYTE, image);

    glBindTexture(GL_TEXTURE_2D, 0);

    stbi_image_free(image);

    // return unique texture identifier
    return textureID;
  }
};
