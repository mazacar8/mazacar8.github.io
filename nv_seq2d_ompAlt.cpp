#include <stdio.h>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "nv_seq2d_ompAlt.h"

FluidBox *FluidBoxCreate2D_ompAlt(int length, int width, float ts) {

	FluidBox *box = new FluidBox;

	box->length = LENGTH;
	box->width = WIDTH;
	box->size = length*width;
	box->numParticles = ceil(length * width * FILL_LEVEL);
	box->time_step_size = DT;
	box->diff_const = DIFF_CONST;

	 
	box->vel_x = new float[length * width]();
	box->temp_vel_x = new float[length * width]();
	box->vel_y = new float[length * width]();
	box->temp_vel_y = new float[length * width]();
	box->pre_x = new float[length * width]();
	box->temp_pre_x = new float[length * width]();
	box->pre_y = new float[length * width]();
	box->temp_pre_y = new float[length * width]();
	box->divergence = new float[length * width]();
	box->grad_x = new float[length * width]();
	box->grad_y = new float[length * width]();
	box->particle = new bool[length * width]();
	box->temp_particle = new bool[length * width]();

	int ctr = 0;

	while(ctr < box->numParticles * FILL_LEVEL){

		box->particle[ctr] = true;
		ctr++;
	}

	return box;

}

void FluidBoxFree2D_ompAlt(FluidBox *box) {


	delete [] box->vel_x;
	delete [] box->vel_y;
	delete [] box->temp_vel_x;
	delete [] box->temp_vel_y;
	delete [] box->pre_x;
	delete [] box->pre_y;
	delete [] box->temp_pre_x;
	delete [] box->temp_pre_y;
	delete [] box->particle;
	delete [] box->divergence;
	delete [] box->grad_x;
	delete [] box->grad_y;
	delete [] box->temp_particle;
	delete box;
}

// void addDensity(FluidBox *box, int x, int y, int z, float density) {

// 	int index = x + (box->length * y) + (box->length * box->width * z);
// 	box->density[index] += density;
// }

void addVelocity2D_ompAlt(FluidBox *box, int x, int y, float vel_x, float vel_y) {

	int length = box->length;
	int width = box->width;

	int index = (y * width) + length;

	box->vel_x[index] += vel_x;
	box->vel_y[index] += vel_y;
}

void advectCube2D_ompAlt(FluidBox *box) {

	float dt = box->time_step_size;
	int old_i,old_j,old_index;
	int length = box->length;
	int width = box->width;

	#pragma omp parallel for
	for(int i = 0; i < length * width; i++){

		int index = i;

		int x = i%length;
		int y = i/length;

		old_i = round(x - box->temp_vel_x[i]*dt);
		old_j = round(y - box->temp_vel_y[i]*dt);

		old_index = old_j*length + old_i;

		if(old_index < box->size and old_index > 0){

			box->temp_particle[index] = box->particle[index];

			if(box->temp_particle[index]){

				box->temp_vel_x[index] = box->vel_x[index];
				box->vel_x[index] = box->vel_x[old_index];
				box->temp_vel_y[index] = box->vel_y[index];
				box->vel_y[index] = box->vel_y[old_index];

			}
		}
	}

	for (int i = 0; i < length*width; ++i)
		box->particle[i] = box->temp_particle[i];
}

void diffuseCube2D_ompAlt(FluidBox *box) {

	int length = box->length;

	float diff = box->diff_const;
	float alpha = (box->length) * (box->width)/(box->time_step_size * diff);
	float beta = alpha + 4;

	float sumx, sumy;

	for(int iter = 0; iter < DIFF_ITER; iter++) {

		copy2dArray_ompAlt(box->temp_vel_x,box->vel_x,box->length,box->width);
		copy2dArray_ompAlt(box->temp_vel_y,box->vel_y,box->length,box->width);

		#pragma omp parallel for
		for (int i = 1; i < box->length * box-> width; ++i){

			int x = i%box->length;
			int y = i/box->length;

			if(x == 0 or x == box->length - 1 or y == 0 or y == box->width - 1){

			}

			else{
				sumx = box->temp_vel_x[((y-1)*length) + x] + box->temp_vel_x[((y+1)*length) + x] +
				   box->temp_vel_x[((y)*length) + x - 1] + box->temp_vel_x[((y)*length) + x + 1];

				sumy = box->temp_vel_y[((y-1)*length) + x] + box->temp_vel_y[((y+1)*length) + x] +
				   box->temp_vel_y[((y)*length) + x - 1] + box->temp_vel_y[((y)*length) + x + 1];

				box->vel_x[i] = (sumx + alpha*box->temp_vel_x[i])/beta;
				box->vel_y[i] = (sumy + alpha*box->temp_vel_y[i])/beta;

			}
		}
	}
}

void addForce2D_ompAlt(FluidBox *box){

	int length = box->length;

	printf("Add Force called at particle %d, %d\n",box->mouse_i,box->mouse_j);
	box->mousePressed = false;

	int i = box->mouse_i;
	int j = box->mouse_j;

	float vel_x = box->vel_x[j*length + i];
	float vel_y = box->vel_y[j*length + i];

	addVelocity2D_ompAlt(box,i,j,vel_x,vel_y);

	for(int m = 0; m < IMPACT_RADIUS; m++) {
		for(int n = 0; n < IMPACT_RADIUS; n++) {

			if(i-m > 0 && j-n > 0 )
				addVelocity2D_ompAlt(box,i-m,j-n,box->vel_x[((j-n)*length) + i - m],box->vel_y[((j-n)*length) + i - m]);

			if(i+m < box->length && j+n < box->width )
				addVelocity2D_ompAlt(box,i+m,j+n,box->vel_x[((j+n)*length) + i + m],box->vel_y[((j+n)*length) + i + m]);
			
		}
	}
}

void computeDivergence2D_ompAlt(FluidBox *box) {

	int length = box->length;

	#pragma omp parallel for
	for (int i = 1; i < box->length * box-> width; ++i){

		int x = i%box->length;
		int y = i/box->length;

		if(x == 0 or x == box->length - 1 or y == 0 or y == box->width - 1) {

		}
		
		else {
			box->divergence[i] = ((box->vel_x[((y+1)*length) + x] - box->vel_x[((y-1)*length) + x]) +
			   (box->vel_y[((y)*length) + x + 1] - box->vel_y[((y)*length) + x - 1]))/2;
		}
	}
}

void projectBox2D_ompAlt(FluidBox *box){

	float alpha = (box->length) * (box->width);
	float beta = 4;
	computeDivergence2D_ompAlt(box);
	float sumx, sumy;
	int length = box->length;

	setZero2D_ompAlt(box->pre_x,box->length,box->width);
	setZero2D_ompAlt(box->pre_y,box->length,box->width);

	for(int iter = 0; iter < DIFF_ITER; iter++) {

		copy2dArray_ompAlt(box->temp_pre_x,box->pre_x,box->length,box->width);
		copy2dArray_ompAlt(box->temp_pre_y,box->pre_y,box->length,box->width);

		#pragma omp parallel for
		for (int i = 1; i < box->length * box-> width; ++i){

			int x = i%box->length;
			int y = i/box->length;

			if(x == 0 or x == box->length - 1 or y == 0 or y == box->width - 1) {

			}
				
			else {	
			   sumx = box->temp_pre_x[((y-1)*length) + x] + box->temp_pre_x[((y+1)*length) + x] +
			   box->temp_pre_x[((y)*length) + x - 1] + box->temp_pre_x[((y)*length) + x + 1];

			   sumy = box->temp_pre_y[((y-1)*length) + x] + box->temp_pre_y[((y+1)*length) + x] +
			   box->temp_pre_y[((y)*length) + x - 1] + box->temp_pre_y[((y)*length) + x + 1];

				box->pre_x[i] = (sumx + alpha*box->divergence[i])/beta;
				box->pre_y[i] = (sumy + alpha*box->divergence[i])/beta;
			}
		}
	}

}

void setZero2D_ompAlt(float* array, int length, int width){

	#pragma omp parallel for
	for(int i = 0; i < length * width; i++)
		array[i] = 0.0f;

}

void copy2dArray_ompAlt(float* dst,float* src, int length, int width){

	#pragma omp parallel for
	for(int i = 0; i < length * width; i++)
		dst[i] = src[i];

}

void accountForGradient2D_ompAlt(FluidBox *box) {

	int length = box->length;
	#pragma omp parallel for
	for(int i = 1; i < box->length * box->width; i++){

		int x = i%box->length;
		int y = i/box->length;

		if(x == 0 or x == box->length - 1 or y == 0 or y == box->width - 1) {

		}
		
		else {
			box->grad_x[i] = (box->pre_x[((y+1)*length) + x] - box->pre_x[((y-1)*length) + x])/2;
			box->grad_y[i] = (box->pre_y[((y)*length) + x + 1] - box->pre_y[((y)*length) + x - 1])/2;

			box->temp_vel_x[i] -= box->grad_x[i];
			box->temp_vel_y[i] -= box->grad_y[i];

			box->vel_x[i] -= box->grad_x[i];
			box->vel_y[i] -= box->grad_y[i];

		}
	}
}

int countParticles_ompAlt(FluidBox *box){

	int numParticles = 0;

	#pragma omp parallel for
	for (int i = 1; i < box->length * box->width; i++){

		if(box->particle[i]) {
			#pragma omp critical 
			{
				numParticles ++;
			}
		}
	}

	return numParticles;

}

void timeStep2D_ompAlt(FluidBox *box){


	advectCube2D_ompAlt(box);
	// printf("Advected\n");
	diffuseCube2D_ompAlt(box);
	// printf("Diffused\n");

	if(box->mousePressed && box->particle[box->mouse_j*box->length + box->mouse_i])
		addForce2D_ompAlt(box);

	projectBox2D_ompAlt(box);
	// printf("Projected\n");
	accountForGradient2D_ompAlt(box);
	// printf("Done\n");

	// int numParticles = countParticles(box);
	// printf("%d\n",numParticles);

}






