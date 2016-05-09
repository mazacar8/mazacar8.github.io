#include <stdio.h>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "nv_seq2d.h"
#include "cycleTimer.h"

FluidBox *FluidBoxCreate2D(int length, int width, float ts) {

	FluidBox *box = new FluidBox;

	box->length = LENGTH;
	box->width = WIDTH;
	box->size = length*width;
	box->numParticles = ceil(length * width * FILL_LEVEL);
	box->time_step_size = DT;
	box->diff_const = DIFF_CONST;


	box->vel_x = new float*[width];
	box->temp_vel_x = new float*[width];
	box->vel_y = new float*[width];
	box->temp_vel_y = new float*[width];
	box->pre_x = new float*[width];
	box->temp_pre_x = new float*[width];
	box->pre_y = new float*[width];
	box->temp_pre_y = new float*[width];
	box->divergence = new float*[width];
	box->grad_x = new float*[width];
	box->grad_y = new float*[width];
	box->particle = new bool*[width];
	box->temp_particle = new bool*[length];

	for (int i = 0; i < width; ++i) {
	 
		box->vel_x[i] = new float[length]();
		box->temp_vel_x[i] = new float[length]();
		box->vel_y[i] = new float[length]();
		box->temp_vel_y[i] = new float[length]();
		box->pre_x[i] = new float[length]();
		box->temp_pre_x[i] = new float[length]();
		box->pre_y[i] = new float[length]();
		box->temp_pre_y[i] = new float[length]();
		box->divergence[i] = new float[length]();
		box->grad_x[i] = new float[length]();
		box->grad_y[i] = new float[length]();
		box->particle[i] = new bool[length];
		box->temp_particle[i] = new bool[length];

	}

	int ctr = 0;

	while(ctr < box->numParticles * FILL_LEVEL){
		box->particle[ctr/length][ctr%length] = true;
		ctr++;
	}

	return box;

}

void FluidBoxFree2D(FluidBox *box) {


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
	delete box;
}

// void addDensity(FluidBox *box, int x, int y, int z, float density) {

// 	int index = x + (box->length * y) + (box->length * box->width * z);
// 	box->density[index] += density;
// }

void addVelocity2D(FluidBox *box, int x, int y, float vel_x, float vel_y) {

	box->vel_x[x][y] += vel_x;
	box->vel_y[x][y] += vel_y;
}

void advectCube2D(FluidBox *box) {

	float dt = box->time_step_size;
	int old_i,old_j,old_index;
	int length = box->length;
	int width = box->width;

	for(int i = 0; i < length; i++){

		//#pragma omp parallel for
		for(int j = 0; j <width; j++){

			old_i = round(i - box->temp_vel_x[i][j]*dt);
			old_j = round(j - box->temp_vel_y[i][j]*dt);

			old_index = old_i*length + old_j;

			if(old_index < box->size and old_index > 0){

				box->particle[i][j] = box->particle[old_i][old_j];

				if(box->particle[i][j]){

					box->temp_vel_x[i][j] = box->vel_x[i][j];
					box->vel_x[i][j] = box->vel_x[old_i][old_j];
					box->temp_vel_y[i][j] = box->vel_y[i][j];
					box->vel_y[i][j] = box->vel_y[old_i][old_j];

				}
					
			}

		}
	}

}

void diffuseCube2D(FluidBox *box) {

	float diff = box->diff_const;
	float alpha = (box->length) * (box->width)/(box->time_step_size * diff);
	float beta = alpha + 4;

	float sumx, sumy;

	for(int iter = 0; iter < DIFF_ITER; iter++) {

		copy2dArray(box->temp_vel_x,box->vel_x,box->length,box->width);
		copy2dArray(box->temp_vel_y,box->vel_y,box->length,box->width);

		for (int i = 1; i < box->length - 1; i++){

			for (int j = 1; j < box->width - 1; j++) {

				sumx = box->temp_vel_x[i-1][j] + box->temp_vel_x[i+1][j] +
					   box->temp_vel_x[i][j-1] + box->temp_vel_x[i][j+1];

				sumy = box->temp_vel_y[i-1][j] + box->temp_vel_y[i+1][j] +
					   box->temp_vel_y[i][j-1] + box->temp_vel_y[i][j+1];

				box->vel_x[i][j] = (sumx + alpha*box->temp_vel_x[i][j])/beta;
				box->vel_y[i][j] = (sumy + alpha*box->temp_vel_y[i][j])/beta;

			}
		}

	}

}

void addForce2D(FluidBox *box){

	printf("Add Force called at particle %d, %d\n",box->mouse_i,box->mouse_j);
	box->mousePressed = false;

	int i = box->mouse_i;
	int j = box->mouse_j;

	float vel_x = box->vel_x[i][j];
	float vel_y = box->vel_y[i][j];
	addVelocity2D(box,i,j,vel_x,vel_y);

	for(int m = 0; m < IMPACT_RADIUS; m++) {
		for(int n = 0; n < IMPACT_RADIUS; n++) {

			if(i-m > 0 && j-n > 0 )
				addVelocity2D(box,i-m,j-n,box->vel_x[i-m][j-n],box->vel_y[i-m][j-n]);

			if(i+m < box->length && j+n < box->width )
				addVelocity2D(box,i+m,j+n,box->vel_x[i+m][j+n],box->vel_y[i+m][j+n]);
			
		}
	}

	// for (int i = 1; i < box->length - 1; i++){

	// 	for (int j = 1; j < box->width - 1; j++) {

	// 		box->vel_x[i][j] = 3;
	// 		box->vel_y[i][j] = 3;

	// 	}
	// }


}

void computeDivergence2D(FluidBox *box) {

	for (int i = 1; i < box->length - 1; i++){

		for (int j = 1; j < box->width - 1; j++) {

				box->divergence[i][j] = ((box->vel_x[i+1][j] - box->vel_x[i-1][j]) +
					   (box->vel_y[i][j+1] - box->vel_y[i][j-1]))/2;

			
		}

	}

}

void projectBox2D(FluidBox *box){

	float alpha = (box->length) * (box->width);
	float beta = 4;
	computeDivergence2D(box);
	float sumx, sumy;

	setZero2D(box->pre_x,box->length,box->width);
	setZero2D(box->pre_y,box->length,box->width);

	for(long i = 0; i < 1000000000; i++){
        
        ;

    }


	for(int iter = 0; iter < DIFF_ITER; iter++) {

		copy2dArray(box->temp_pre_x,box->pre_x,box->length,box->width);
		copy2dArray(box->temp_pre_y,box->pre_y,box->length,box->width);
		
		for (int i = 1; i < box->length - 1; i++){

			for (int j = 1; j < box->width - 1; j++) {


				sumx = box->temp_pre_x[i-1][j] + box->temp_pre_x[i+1][j] +
					   box->temp_pre_x[i][j-1] + box->temp_pre_x[i][j+1];

				sumy = box->temp_pre_y[i-1][j] + box->temp_pre_y[i+1][j] +
					   box->temp_pre_y[i][j-1] + box->temp_pre_y[i][j+1];

				box->pre_x[i][j] = (sumx + alpha*box->divergence[i][j])/beta;
				box->pre_y[i][j] = (sumy + alpha*box->divergence[i][j])/beta;

				
			}

		}

	}

}

void setZero2D(float** array, int length, int width){

	for(int i = 0; i < length; i++){

		for(int j = 0; j < width; j++){

			array[i][j] = 0.0f;
			
		}
			
	}

}

void copy2dArray(float** dst,float** src, int length, int width){

	for(int i = 0; i < length; i++){

		for(int j = 0; j < width; j++){
			
			dst[i][j] = src[i][j];

		}
			
	}

}

void accountForGradient2D(FluidBox *box) {

	for (int i = 1; i < box->length - 1; i++){

		for (int j = 1; j < box->width - 1; j++) {

			box->grad_x[i][j] = (box->pre_x[i+1][j] - box->pre_x[i-1][j])/2;
			box->grad_y[i][j] = (box->pre_y[i][j+1] - box->pre_y[i][j-1])/2;

			box->temp_vel_x[i][j] -= box->grad_x[i][j];
			box->temp_vel_y[i][j] -= box->grad_y[i][j];

			box->vel_x[i][j] -= box->grad_x[i][j];
			box->vel_y[i][j] -= box->grad_y[i][j];

			
		}

	}

}

int countParticles(FluidBox *box){

	int numParticles = 0;

	for (int i = 1; i < box->length - 1; i++){

		for (int j = 1; j < box->width - 1; j++) {

			if(box->particle[i][j])
				numParticles++;

		}
	}

	return numParticles;

}

void timeStep2D(FluidBox *box){

	double startTime = CycleTimer::currentSeconds();
	advectCube2D(box);
	// printf("Advected\n");
	diffuseCube2D(box);
	// printf("Diffused\n");

	if(box->mousePressed && box->particle[box->mouse_i][box->mouse_j])
		addForce2D(box);

	projectBox2D(box);
	// printf("Projected\n");
	accountForGradient2D(box);
	double endTime = CycleTimer::currentSeconds();

    printf("Sequential Version takes %f ms\n",(endTime - startTime)*1000);
	// printf("Done\n");

	// int numParticles = countParticles(box);
	// printf("%d\n",numParticles);

}






