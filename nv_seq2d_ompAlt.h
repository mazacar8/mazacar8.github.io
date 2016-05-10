#ifndef __NV_SEQ2D_OMPALT_H__
#define __NV_SEQ2D_OMPALT_H__

#define FILL_LEVEL 0.5
#define DIFF_ITER 50
#define IMPACT_RADIUS 3
#define LENGTH 1000
#define WIDTH 1000
#define DT 0.015
#define DIFF_CONST 0.0016

typedef struct {

	int length;
	int width;

	float time_step_size;
	float diff_const;

	int size;
	int numParticles;

	float* vel_x;
	float* vel_y;

	float* temp_vel_x;
	float* temp_vel_y;

	float* pre_x;
	float* pre_y;

	float* temp_pre_x;
	float* temp_pre_y;

	float* grad_x;
	float* grad_y;

	float* divergence;

	bool* particle;
	bool* temp_particle;

	int mouse_i, mouse_j;
	bool mousePressed;

} FluidBoxAlt;

FluidBoxAlt *FluidBoxCreate2D_ompAlt(int length, int width, float ts);
void FluidBoxFree2D_ompAlt(FluidBoxAlt *box);
void addVelocity2D_ompAlt(FluidBoxAlt *box, int x, int y, float vel_x, float vel_y);
void advectCube2D_ompAlt(FluidBoxAlt *box);
void diffuseCube2D_ompAlt(FluidBoxAlt *box);
void addForce2D_ompAlt(FluidBoxAlt *box);
void computeDivergence2D_ompAlt(FluidBoxAlt *box);
void projectBox2D_ompAlt(FluidBoxAlt *box);
void accountForGradient2D_ompAlt(FluidBoxAlt *box);
void timeStep2D_ompAlt(FluidBoxAlt *box);
void copy2dArray_ompAlt(float* dst,float* src, int length, int width);
void setZero2D_ompAlt(float* array, int length, int width);
int countParticles_ompAlt(FluidBoxAlt *box);

#endif



