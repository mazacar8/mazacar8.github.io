/*
* OGL02Animation.cpp: 3D Shapes with animation
*/

#include "platformgl.h" 
#include <math.h>
#include <omp.h>
//#include "nv_seq.h"
#include "nv_seq2d.h"
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
#include <getopt.h>
#include "cudaRenderer.h"

#define WINDOW_HEIGHT 350
#define WINDOW_WIDTH 350

void startRendererWithDisplay(CudaRenderer* renderer);

void usage(const char* progname) {
    printf("Usage: %s [options]\n", progname);
    printf("Program Options:\n");
    printf("  -r  --renderer <ref/cuda>  Select renderer: ref or cuda\n");
    printf("  -?  --help                 This message\n");
}

struct Point
{
    float x, y;
    unsigned char r, g, b, a;
} ;

std::vector< Point > points;

/* Global variables */
char title[] = "Fluid Simulation";
GLfloat angleCube = 0.0f;     // Rotational angle for cube [NEW]
int refreshMills = DT * 1000;        // refresh interval in milliseconds [NEW]
GLfloat camera_angle = 0.0f;
bool mousePressed = false;

FluidBox *box;

/* Initialize OpenGL Graphics */
void initGL() {
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f); // Set background color to white and opaque
	glClearDepth(1.0f);                   // Set background depth to farthest
	glEnable(GL_DEPTH_TEST);   // Enable depth testing for z-culling
	glEnable(GL_COLOR_MATERIAL);

	// glEnable(GL_LIGHTING);
	//glEnable(GL_LIGHT0);
	//glEnable(GL_LIGHT1);
	//glEnable(GL_NORMALIZE);

	glDepthFunc(GL_LEQUAL);    // Set the type of depth-test
	glShadeModel(GL_SMOOTH);   // Enable smooth shading
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);  // Nice perspective corrections
}

void drawParticle(int i, int j){

 	GLfloat x = WINDOW_WIDTH*j/(box->width);
 	GLfloat y = WINDOW_HEIGHT*i/(box->length);

 	glBegin(GL_POINTS);

 		if(i < box->length/2 && j < box->length/2)
 			glColor3f(0.4f,0.84f,0.91f);

 		else if(i < box->length/2 && j < box->length/2)
 			glColor3f(0.1f,0.5f,0.3f);
 		else if(i < box->length/2 && j < box->length/2)
 			glColor3f(0.8f,0.6f,0.2f);
 		else
 			glColor3f(0.5f,0.3f,0.9f);

 		glVertex2f(x,y);
 		glColor3f(1.0f,1.0f,1.0f);
 	glEnd();


}

void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-50, 50, -50, 50, -1, 1);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    // draw
    glColor3ub( 255, 255, 255 );
    glEnableClientState( GL_VERTEX_ARRAY );
    glEnableClientState( GL_COLOR_ARRAY );
    glVertexPointer( 2, GL_FLOAT, sizeof(Point), &points[0].x );
    glColorPointer( 4, GL_UNSIGNED_BYTE, sizeof(Point), &points[0].r );
    glPointSize( 5.0 );
    glDrawArrays( GL_POINTS, 0, points.size() );
    glDisableClientState( GL_VERTEX_ARRAY );
    glDisableClientState( GL_COLOR_ARRAY );

    glFlush();
    glutSwapBuffers();
}

/* Called back when timer expired [NEW] */
void timer(int value) {

	glutPostRedisplay();      // Post re-paint request to activate display()

	timeStep2D(box);

	points.clear();

	for(int i = 0; i < box->length; i++){
		for(int j = 0; j < box->width; j++){

			if(box->particle[i][j]){

				Point pt;
				pt.x = -50+WINDOW_WIDTH*j/(box->width);
				pt.y = -50+WINDOW_HEIGHT*i/(box->length);
				pt.r = 100;
				pt.g = 200;
				pt.b = 230;
				pt.a = 255;
				points.push_back(pt);
			}

		}
	}

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

void handleMousePress(int button, int state, int x, int y){

	if(button == GLUT_LEFT_BUTTON && state == GLUT_DOWN){

		box->mousePressed = true;
		box->mouse_i = round(((float)x)/WINDOW_WIDTH * box->width);
		box->mouse_j = box->length - round(((float)y)/WINDOW_HEIGHT * box->length);

	}

}

int main(int argc, char** argv) {

	bool cuda = false;

	 // parse commandline options ////////////////////////////////////////////
    int opt;
    static struct option long_options[] = {
        {"help",     0, 0,  '?'},
        {"renderer", 1, 0,  'r'},
        {0 ,0, 0, 0}
    };

    while ((opt = getopt_long(argc, argv, "r:?", long_options, NULL)) != EOF) {

        switch (opt) {
        case 'r':
            if (std::string(optarg).compare("cuda") == 0)
                cuda = true;

            else if(std::string(optarg).compare("ref") == 0)
            	cuda = false;

            else{
            	usage(argv[0]);
            	return 1;
            }

            break;
        case '?':
            usage(argv[0]);
            break;
        default:
            usage(argv[0]);
            return 1;
        }
    }
    // end parsing of commandline options //////////////////////////////////////

    if(cuda){
        CudaRenderer* renderer = new CudaRenderer();
        renderer->allocOutputImage(WINDOW_WIDTH, WINDOW_HEIGHT);
        SceneName sceneName = WATER_CUBE;
        renderer->loadScene(sceneName);
        renderer->setup();

        glutInit(&argc, argv);
        startRendererWithDisplay(renderer);
    }

    else{

    	box = FluidBoxCreate2D(LENGTH, WIDTH, DT);
		glutInit(&argc, argv);            // Initialize GLUT
		glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH | GLUT_DOUBLE); // Enable double buffered mode
		glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);   // Set the window's initial width & height
		glutCreateWindow(title);          // Create window with the given title
		glutDisplayFunc(display);       // Register callback handler for window re-paint event
		glutReshapeFunc(reshape);       // Register callback handler for window re-size event
		glutMouseFunc(handleMousePress);
		initGL();                       // Our own OpenGL initialization
		glutTimerFunc(0, timer, 0);     // First timer call immediately [NEW]
		glutMainLoop();                 // Enter the infinite event-processing loop
	}


	return 0;

}