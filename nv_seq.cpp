#include <stdio.h>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "nv_seq.h"

FluidBox1D *FluidBoxCreate(int length, int width,float ts,float diff) {

	FluidBox1D *box = new FluidBox1D;

	box->length = length;
	box->width = width;
	box->size = length*width;
	box->numParticles = ceil(length * width * FILL_LEVEL);
	box->diff_const = diff;

	box->time_step_size = ts;


	box->vel_x = new float[box->size];
	box->temp_vel_x = new float[box->size];
	box->vel_y = new float[box->size];
	box->temp_vel_y = new float[box->size];
	box->pre_x = new float[box->size];
	box->temp_pre_x = new float[box->size];
	box->pre_y = new float[box->size];
	box->temp_pre_y = new float[box->size];
	box->divergence = new float[box->size];
	box->grad_x = new float[box->size];
	box->grad_y = new float[box->size];
	box->particle = new bool[box->size];
	box->temp_particle = new bool[box->size];

	for(int i = 0; i < length/2; ++i)
		box->particle[i] = true;

	return box;

}