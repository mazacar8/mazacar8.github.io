#ifndef __NV_SEQ2D_H__
#define __NV_SEQ2D_H__

#define FILL_LEVEL 0.75
#define DIFF_ITER 50
#define IMPACT_RADIUS 3
#define LENGTH 10
#define WIDTH 10
#define DT 0.015
#define DIFF_CONST 0.0016

typedef struct {

	int length;
	int width;

	float time_step_size;
	float diff_const;

	int size;
	int numParticles;

	float** vel_x;
	float** vel_y;

	float** temp_vel_x;
	float** temp_vel_y;

	float** pre_x;
	float** pre_y;

	float** temp_pre_x;
	float** temp_pre_y;

	float** grad_x;
	float** grad_y;

	float** divergence;

	bool** particle;
	bool** temp_particle;

	int mouse_i, mouse_j;
	bool mousePressed;

} FluidBox;

FluidBox *FluidBoxCreate2D(int length, int width, float ts);
void FluidBoxFree2D(FluidBox *box);
void addVelocity2D(FluidBox *box, int x, int y, float vel_x, float vel_y);
void advectCube2D(FluidBox *box);
void diffuseCube2D(FluidBox *box);
void addForce2D(FluidBox *box);
void computeDivergence2D(FluidBox *box);
void projectBox2D(FluidBox *box);
void accountForGradient2D(FluidBox *box);
void timeStep2D(FluidBox *box);
void copy2dArray(float** dst,float** src, int length, int width);
void setZero2D(float** array, int length, int width);
int countParticles(FluidBox *box);

#endif



