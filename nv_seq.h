#ifndef __NV_SEQ_H__
#define __NV_SEQ_H__

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

	int numParticles;
	int size;

	float* vel_x;
	float* vel_y;

	float* temp_vel_x;
	float* temp_vel_y;
	float* temp_vel_z;

	float* pre_x;
	float* pre_y;

	float* temp_pre_x;
	float* temp_pre_y;

	float* grad_x;
	float* grad_y;

	float* divergence;

	bool* particle;
	bool* temp_particle;

} FluidBox1D;

FluidBox1D *FluidBoxCreate(int length, int width, float ts, float diff);

#endif


