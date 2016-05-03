#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "platformgl.h"

GLfloat xRotated, yRotated, zRotated;
GLdouble radius=1;
GLenum doubleBuffer;

static void
Reshape(int x, int y)
{

  if (y == 0 || x == 0) return;   
    glMatrixMode(GL_PROJECTION);  
    glLoadIdentity(); 
    gluPerspective(39.0,(GLdouble)x/(GLdouble)y,0.6,21.0);
    glMatrixMode(GL_MODELVIEW);
    glViewport(0,0,x,y);
}

static void
Draw(void)
{

  glMatrixMode(GL_MODELVIEW);
  // clear the drawing buffer.
  glClear(GL_COLOR_BUFFER_BIT);
  // clear the identity matrix.
  glLoadIdentity();
  // traslate the draw by z = -4.0
  // Note this when you decrease z like -8.0 the drawing will looks far , or smaller.
  glTranslatef(0.0,0.0,-5.0);
  // Red color used to draw.
  glColor3f(0.9, 0.3, 0.2); 
  // changing in transformation matrix.
  // rotation about X axis
  glRotatef(xRotated,1.0,0.0,0.0);
  // rotation about Y axis
  glRotatef(yRotated,0.0,1.0,0.0);
  // rotation about Z axis
  glRotatef(zRotated,0.0,0.0,1.0);
  // scaling transfomation 
  glScalef(1.0,1.0,1.0);
  // built-in (glut library) function , draw you a sphere.
  glutSolidSphere(radius,20,20);
  // Flush buffers to screen
   
  glFlush();        
 
}

int
main(int argc, char **argv)
{
  GLenum type;

  glutInit(&argc, argv);
  //Args(argc, argv);

  type = GLUT_RGB | GLUT_DEPTH;
  type |= (doubleBuffer) ? GLUT_DOUBLE : GLUT_SINGLE;
  glutInitDisplayMode(type);
  glutInitWindowSize(350, 350);
  glutCreateWindow("High Resolution Fluid Simulation");
  xRotated = yRotated = zRotated = 30.0;
  xRotated=43;
  yRotated=50;

  glutReshapeFunc(Reshape);
  // glutKeyboardFunc(Key);
  // glutSpecialFunc(SpecialKey);
  glutDisplayFunc(Draw);
  glutMainLoop();
  return 0;
}