/*
* OGL02Animation.cpp: 3D Shapes with animation
*/

#include "platformgl.h" 
#include <math.h>
#include <omp.h>
#include "nv_seq.h"
#include <stdio.h>

/* Global variables */
char title[] = "3D Shapes with animation";
GLfloat angleCube = 0.0f;     // Rotational angle for cube [NEW]
int refreshMills = 15;        // refresh interval in milliseconds [NEW]
GLfloat camera_angle = 0.0f;

FluidBox *box = FluidBoxCreate(100, 100, 100, refreshMills/1000.0);

/* Initialize OpenGL Graphics */
void initGL() {
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Set background color to black and opaque
	glClearDepth(1.0f);                   // Set background depth to farthest
	glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_LIGHT1);
	glEnable(GL_NORMALIZE);
	glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
	glShadeModel(GL_SMOOTH);   // Enable smooth shading
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections
}

// //https://gist.github.com/strife25/803118
// void drawFilledCircle(GLdouble x, GLdouble y, GLdouble z, GLfloat radius){
// 	int i;
// 	int triangleAmount = 24; //# of triangles used to draw circle

// 	//GLfloat radius = 0.8f; //radius
// 	GLfloat twicePi = 2.0f * M_PI;

// 	glBegin(GL_TRIANGLE_FAN);
// 	glColor3f(0.0f, 0.7f, 0.8f);
// 	glVertex3f(x, y, z); // center of circle
// 	for(i = 0; i <= triangleAmount;i++) { 
// 		glVertex2f(
// 							x + (radius * cos(i *  twicePi / triangleAmount)), 
// 				y + (radius * sin(i * twicePi / triangleAmount))
// 		);
// 	}
// 	glColor3f(0.0f, 0.0f, 0.0f);
// 	glEnd();
// }

// void draw() {

// 	glEnable(GL_TEXTURE_2D);
// 	glBindTexture(GL_TEXTURE_2D, textureId);
// 	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
// 	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

// 	glBegin(GL_QUADS);
// 	for(unsigned int i = 0; i < ps.size(); i++) {

// 		Particle* p = ps[i];
// 		glColor4f(p->color[0], p->color[1], p->color[2],
// 				  (1 - p->timeAlive / p->lifespan));
// 		float size = PARTICLE_SIZE / 2;
	
// 		Vec3f pos = adjParticlePos(p->pos);
	
// 		glTexCoord2f(0, 0);
// 		glVertex3f(pos[0] - size, pos[1] - size, pos[2]);
// 		glTexCoord2f(0, 1);
// 		glVertex3f(pos[0] - size, pos[1] + size, pos[2]);
// 		glTexCoord2f(1, 1);
// 		glVertex3f(pos[0] + size, pos[1] + size, pos[2]);
// 		glTexCoord2f(1, 0);
// 		glVertex3f(pos[0] + size, pos[1] - size, pos[2]);
// 	}
// 	glEnd();
// }


void drawParticle(int i, int j, int k){

 	float half_length = 1.0/(box->length);

 	GLfloat x = -1.0 + 2.0*i/(box->depth);
 	GLfloat y = -1.0 + 2.0*j/(box->width);
 	GLfloat z = -1.0 + 2.0*k/(box->length);

 	glBegin(GL_QUADS);

 		glColor3f(0.4f,0.84f,0.91f);
 		glVertex3f(x-half_length,y,z-half_length);
 		glVertex3f(x-half_length,y,z+half_length);
 		glVertex3f(x+half_length,y,z-half_length);
 		glVertex3f(x+half_length,y,z+half_length);

 	glEnd();


}

/* Handler for window-repaint event. Called back when the window first appears and
whenever the window needs to be re-painted. */
void display() {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); // Clear color and depth buffers
	glMatrixMode(GL_MODELVIEW);     // To operate on model-view matrix

	// Render a color-cube consisting of 6 quads with different colors
	glLoadIdentity();                 // Reset the model-view matrix


	GLfloat ambientColor[] = {0.2f, 0.2f, 0.2f, 0.2f};
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambientColor);

	//Add positioned light
	GLfloat lightColor0[] = {0.5f, 0.5f, 0.5f, 1.0f}; //Color (0.5, 0.5, 0.5)
	GLfloat lightPos0[] = {4.0f, 0.0f, 8.0f, 1.0f}; //Positioned at (4, 0, 8)
	glLightfv(GL_LIGHT0, GL_DIFFUSE, lightColor0);
	glLightfv(GL_LIGHT0, GL_POSITION, lightPos0);

	//Add directed light
	GLfloat lightColor1[] = {0.5f, 0.2f, 0.2f, 1.0f}; //Color (0.5, 0.2, 0.2)
	//Coming from the direction (-1, 0.5, 0.5)
	GLfloat lightPos1[] = {-1.0f, 0.5f, 0.5f, 0.0f};
	glLightfv(GL_LIGHT1, GL_DIFFUSE, lightColor1);
	glLightfv(GL_LIGHT1, GL_POSITION, lightPos1);

	glRotatef(-camera_angle, 0.0f, 0.0f, 1.0f);

	glTranslatef(0.0f, 0.0f, -7.0f);  // Move right and into the screen
	glRotatef(angleCube, 0.0f, 1.0f, 0.0f); 
	glPushMatrix(); 
	glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

	glColor3f(1.0f, 1.0f, 1.0f);
	glBegin(GL_QUADS);                // Begin drawing the color cube with 6 quads
		// Top face (y = 1.0f)
		// Define vertices in counter-clockwise (CCW) order with normal pointing out
		// glColor3f(0.0f, 1.0f, 0.0f);     // Green
		glNormal3f(0.0f, 1.0f, 0.0f);
		glVertex3f( 1.0f, 1.0f, -1.0f);
		glVertex3f(-1.0f, 1.0f, -1.0f);
		glVertex3f(-1.0f, 1.0f,  1.0f);
		glVertex3f( 1.0f, 1.0f,  1.0f);

		// Bottom face (y = -1.0f)
		// glColor3f(1.0f, 0.5f, 0.0f);     // Orange
		glNormal3f(0.0f, -1.0f, 0.0f);
		glVertex3f( 1.0f, -1.0f,  1.0f);
		glVertex3f(-1.0f, -1.0f,  1.0f);
		glVertex3f(-1.0f, -1.0f, -1.0f);
		glVertex3f( 1.0f, -1.0f, -1.0f);

		// Front face  (z = 1.0f)
		// glColor3f(1.0f, 0.0f, 0.0f);     // Red
		glNormal3f(0.0f, 0.0f, 1.0f);
		glVertex3f( 1.0f,  1.0f, 1.0f);
		glVertex3f(-1.0f,  1.0f, 1.0f);
		glVertex3f(-1.0f, -1.0f, 1.0f);
		glVertex3f( 1.0f, -1.0f, 1.0f);

		// Back face (z = -1.0f)
		// glColor3f(1.0f, 1.0f, 0.0f);     // Yellow
		glNormal3f(0.0f, 0.0f, -1.0f);
		glVertex3f( 1.0f, -1.0f, -1.0f);
		glVertex3f(-1.0f, -1.0f, -1.0f);
		glVertex3f(-1.0f,  1.0f, -1.0f);
		glVertex3f( 1.0f,  1.0f, -1.0f);

		// Left face (x = -1.0f)
		// glColor3f(0.0f, 0.0f, 1.0f);     // Blue
		glNormal3f(-1.0f, 0.0f, 0.0f);
		glVertex3f(-1.0f,  1.0f,  1.0f);
		glVertex3f(-1.0f,  1.0f, -1.0f);
		glVertex3f(-1.0f, -1.0f, -1.0f);
		glVertex3f(-1.0f, -1.0f,  1.0f);

		// Right face (x = 1.0f)
		// glColor3f(1.0f, 0.0f, 1.0f);     // Magenta
		glNormal3f(1.0f, 0.0f, 0.0f);
		glVertex3f(1.0f,  1.0f, -1.0f);
		glVertex3f(1.0f,  1.0f,  1.0f);
		glVertex3f(1.0f, -1.0f,  1.0f);
		glVertex3f(1.0f, -1.0f, -1.0f);
	glEnd();
	glPopMatrix();

	//#pragma omp parallel for
	for(int i = 0; i < box->depth; i++){

		//#pragma omp parallel for
		for(int j = 0; j < box->width; j++){

			//#pragma omp parallel for
			for(int k = 0; k < box->length; k++){

				if(box->particle[i][j][k]){
					drawParticle(i,j,k);
				}
			}
		}

	}

	glutSwapBuffers();  // Swap the front and back frame buffers (double buffering)

}

void wrapper() {

	display();

}

/* Called back when timer expired [NEW] */
void timer(int value) {

	glutPostRedisplay();      // Post re-paint request to activate display()

	timeStep(box);

	// camera_angle += 0.5;
	// if(camera_angle > 360)
	// 	camera_angle -= 360;

	glutTimerFunc(refreshMills, timer, 0); // next timer call milliseconds later

}

/* Handler for window re-size event. Called back when the window first appears and
whenever the window is re-sized with its new width and height */
void reshape(GLsizei width, GLsizei height) {  // GLsizei for non-negative integer
// Compute aspect ratio of the new window
	if (height == 0) height = 1;                // To prevent divide by 0
		GLfloat aspect = (GLfloat)width / (GLfloat)height;

	// Set the viewport to cover the new window
	glViewport(0, 0, width, height);

	// Set the aspect ratio of the clipping volume to match the viewport
	glMatrixMode(GL_PROJECTION);  // To operate on the Projection matrix
	glLoadIdentity();             // Reset
	// Enable perspective projection with fovy, aspect, zNear and zFar
	gluPerspective(45.0f, aspect, 1.0f, 100.0f);
}


/* Main function: GLUT runs as a console application starting at main() */
int main(int argc, char** argv) {

	glutInit(&argc, argv);            // Initialize GLUT
	glutInitDisplayMode(GLUT_DOUBLE); // Enable double buffered mode
	glutInitWindowSize(1000, 800);   // Set the window's initial width & height
	glutInitWindowPosition(50, 50); // Position the window's initial top-left corner
	glutCreateWindow(title);          // Create window with the given title
	glutDisplayFunc(wrapper);       // Register callback handler for window re-paint event
	glutReshapeFunc(reshape);       // Register callback handler for window re-size event
	initGL();                       // Our own OpenGL initialization
	glutTimerFunc(0, timer, 0);     // First timer call immediately [NEW]
	glutMainLoop();                 // Enter the infinite event-processing loop
	return 0;

}