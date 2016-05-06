#include <stdio.h>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define FILL_LEVEL 0.75
#define DIFF_ITER 50
#define IMPACT_RADIUS 3

typedef struct {

	int length;
	int width;
	int depth;

	float time_step_size;
	float diff_const;

	int numParticles;
	int size;

	float*** vel_x;
	float*** vel_y;
	float*** vel_z;

	float*** temp_vel_x;
	float*** temp_vel_y;
	float*** temp_vel_z;

	float*** pre_x;
	float*** pre_y;
	float*** pre_z;

	float*** temp_pre_x;
	float*** temp_pre_y;
	float*** temp_pre_z;

	float*** grad_x;
	float*** grad_y;
	float*** grad_z;

	float*** divergence;

	bool *particle;

} FluidBox;


FluidBox *FluidBoxCreate(int length, int width, int depth, float ts) {

	FluidBox *box = new FluidBox;

	box->length = length;
	box->width = width;
	box->depth = depth;
	box->size = length*width*depth;
	box->numParticles = ceil(length * width * depth * FILL_LEVEL);

	box->time_step_size = ts;

	box->particle = new bool [box->size];

	for(int i = 0; i < box->size; i += 1) {

		if(i < box->numParticles)
			box->particle[i] = true;

		else
			box->particle[i] = false;

	}

   box->vel_x = new float**[depth];
   box->temp_vel_x = new float**[depth];
   box->vel_y = new float**[depth];
   box->temp_vel_y = new float**[depth];
   box->vel_z = new float**[depth];
   box->temp_vel_z = new float**[depth];
   box->pre_x = new float**[depth];
   box->temp_pre_x = new float**[depth];
   box->pre_y = new float**[depth];
   box->temp_pre_y = new float**[depth];
   box->pre_z = new float**[depth];
   box->temp_pre_z = new float**[depth];
   box->divergence = new float**[depth];
   box->grad_x = new float**[depth];
   box->grad_y = new float**[depth];
   box->grad_z = new float**[depth];

   for (int i = 0; i < depth; ++i) {
   
   		box->vel_x[i] = new float*[width];
   		box->temp_vel_x[i] = new float*[width];
   		box->vel_y[i] = new float*[width];
   		box->temp_vel_y[i] = new float*[width];
   		box->vel_z[i] = new float*[width];
   		box->temp_vel_z[i] = new float*[width];
   		box->pre_x[i] = new float*[width];
   		box->temp_pre_x[i] = new float*[width];
   		box->pre_y[i] = new float*[width];
   		box->temp_pre_y[i] = new float*[width];
   		box->pre_z[i] = new float*[width];
   		box->temp_pre_z[i] = new float*[width];
   		box->divergence[i] = new float*[width];
   		box->grad_x[i] = new float*[width];
   		box->grad_y[i] = new float*[width];
   		box->grad_z[i] = new float*[width];

   		for (int j = 0; j < width; ++j){
      		box->vel_x[i][j] = new float[length]();
      		box->temp_vel_x[i][j] = new float[length]();
      		box->vel_y[i][j] = new float[length]();
      		box->temp_vel_y[i][j] = new float[length]();
      		box->vel_x[i][j] = new float[length]();
      		box->temp_vel_z[i][j] = new float[length]();
      		box->pre_x[i][j] = new float[length]();
      		box->temp_pre_x[i][j] = new float[length]();
      		box->pre_y[i][j] = new float[length]();
      		box->temp_pre_y[i][j] = new float[length]();
      		box->pre_x[i][j] = new float[length]();
      		box->temp_pre_z[i][j] = new float[length]();
      		box->divergence[i][j] = new float[length]();
      		box->grad_x[i][j] = new float[length]();
   			box->grad_y[i][j] = new float[length]();
   			box->grad_z[i][j] = new float[length]();
      	}

    }

    return box;

}

void FluidBoxFree(FluidBox *box) {


	delete [] box->vel_x;
	delete [] box->vel_y;
	delete [] box->vel_z;
	delete [] box->temp_vel_x;
	delete [] box->temp_vel_y;
	delete [] box->temp_vel_z;
	delete [] box->pre_x;
	delete [] box->pre_y;
	delete [] box->pre_z;
	delete [] box->temp_pre_x;
	delete [] box->temp_pre_y;
	delete [] box->temp_pre_z;
	delete [] box->particle;
	delete [] box->divergence;
	delete [] box->grad_x;
	delete [] box->grad_y;
	delete [] box->grad_z;
	delete box;
}

// void addDensity(FluidBox *box, int x, int y, int z, float density) {

// 	int index = x + (box->length * y) + (box->length * box->width * z);
// 	box->density[index] += density;
// }

void addVelocity(FluidBox *box, int x, int y, int z, float vel_x, float vel_y, float vel_z) {

	box->vel_x[x][y][z] += vel_x;
	box->vel_y[x][y][z] += vel_y;
	box->vel_z[x][y][z] += vel_z;
}

void advectCube(FluidBox *box) {

	float dt = box->time_step_size;
	int old_i,old_j,old_k,index,old_index;
	int length = box->length;
	int width = box->width;

	for(int i = 0; i < box->depth; i++){

		for(int j = 0; j < box->width; j++){

			for(int k = 0; k < box->length; k++){

				index = k + (length * j) + (length * width * i);

				old_i = round(i - box->vel_z[i][j][k]*dt);
				old_j = round(j - box->vel_y[i][j][k]*dt);
				old_k = round(k - box->vel_x[i][j][k]*dt);

				old_index = old_k + (length * old_j) + (length * width * old_i);

				if(index < box->size)
					box->particle[index] = box->particle[old_index];

				if(box->particle[index]){

					box->temp_vel_x[index] = box->vel_x[index];
					box->vel_x[index] = box->vel_x[old_index];
					box->temp_vel_y[index] = box->vel_y[index];
					box->vel_y[index] = box->vel_y[old_index];
					box->temp_vel_z[index] = box->vel_z[index];
					box->vel_z[index] = box->vel_z[old_index];

				}

			}
		}
	}

}

void diffuseCube(FluidBox *box) {

	float diff = box->diff_const;
	float alpha = (box->length) * (box->width)/(box->time_step_size * diff);
	float beta = alpha + 6;

	float sumx, sumy, sumz;

	for(int i = 0; i < DIFF_ITER; i++) {

		memcpy((void *)box->temp_vel_x, (void *) box->vel_x, (std::size_t) box->size);
		memcpy((void *)box->temp_vel_y, (void *) box->vel_y, (std::size_t) box->size);
		memcpy((void *)box->temp_vel_z, (void *) box->vel_z, (std::size_t) box->size);

		for (int i = 1; i < box->depth - 1; i++){

			for (int j = 1; j < box->width - 1; j++) {

				for (int k = 1; j < box->length - 1; k++) {

					sumx = box->temp_vel_x[i-1][j][k] + box->temp_vel_x[i+1][j][k] +
						   box->temp_vel_x[i][j-1][k] + box->temp_vel_x[i][j+1][k] +
						   box->temp_vel_x[i][j][k-1] + box->temp_vel_x[i][j][k+1];

					sumy = box->temp_vel_y[i-1][j][k] + box->temp_vel_y[i+1][j][k] +
						   box->temp_vel_y[i][j-1][k] + box->temp_vel_y[i][j+1][k] +
						   box->temp_vel_y[i][j][k-1] + box->temp_vel_y[i][j][k+1];

					sumz = box->temp_vel_z[i-1][j][k] + box->temp_vel_z[i+1][j][k] +
						   box->temp_vel_z[i][j-1][k] + box->temp_vel_z[i][j+1][k] +
						   box->temp_vel_z[i][j][k-1] + box->temp_vel_z[i][j][k+1];

					box->vel_x[i][j][k] = (sumx + alpha*box->temp_vel_x[i][j][k])/beta;
					box->vel_y[i][j][k] = (sumy + alpha*box->temp_vel_y[i][j][k])/beta;
					box->vel_z[i][j][k] = (sumz + alpha*box->temp_vel_z[i][j][k])/beta;

				}
			}

		}

	}

}

void addForce(FluidBox *box, int x, int y, int z,
	float vel_x, float vel_y, float vel_z){

	addVelocity(box,x,y,z,vel_x,vel_y,vel_z);

	for(int i = 0; i < IMPACT_RADIUS; i++) {
		for(int j = 0; j < IMPACT_RADIUS; j++) {
			for(int k = 0; k < IMPACT_RADIUS; k++) {

				addVelocity(box,x-k,y-j,z-i,vel_x, vel_y, vel_z);
				addVelocity(box,x+k,y+j,z+i,vel_x, vel_y, vel_z);
			}
		}
	}

}

void computeDivergence(FluidBox *box) {

	for (int i = 1; i < box->depth - 1; i++){

		for (int j = 1; j < box->width - 1; j++) {

			for (int k = 1; j < box->length - 1; k++) {

				box->divergence[i][j][k] = (box->vel_z[i+1][j][k] - box->vel_z[i-1][j][k] +
					   box->vel_y[i][j-1][k] - box->vel_y[i][j+1][k] +
					   box->vel_x[i][j][k-1] - box->vel_x[i][j][k+1])/2;

			}
		}

	}

}

void projectBox(FluidBox *box){

	float alpha = (box->length) * (box->width);
	float beta = 6;
	computeDivergence(box);

	float sumx, sumy, sumz;

	memset(box->pre_x,0,box->size);
	memset(box->pre_y,0,box->size);
	memset(box->pre_z,0,box->size);

	for(int i = 0; i < DIFF_ITER; i++) {

		memcpy((void *)box->temp_pre_x, (void *) box->pre_x, (std::size_t) box->size);
		memcpy((void *)box->temp_pre_y, (void *) box->pre_y, (std::size_t) box->size);
		memcpy((void *)box->temp_pre_z, (void *) box->pre_z, (std::size_t) box->size);

		for (int i = 1; i < box->depth - 1; i++){

			for (int j = 1; j < box->width - 1; j++) {

				for (int k = 1; j < box->length - 1; k++) {

					sumx = box->temp_pre_x[i-1][j][k] + box->temp_pre_x[i+1][j][k] +
						   box->temp_pre_x[i][j-1][k] + box->temp_pre_x[i][j+1][k] +
						   box->temp_pre_x[i][j][k-1] + box->temp_pre_x[i][j][k+1];

					sumy = box->temp_pre_y[i-1][j][k] + box->temp_pre_y[i+1][j][k] +
						   box->temp_pre_y[i][j-1][k] + box->temp_pre_y[i][j+1][k] +
						   box->temp_pre_y[i][j][k-1] + box->temp_pre_y[i][j][k+1];

					sumz = box->temp_pre_z[i-1][j][k] + box->temp_pre_z[i+1][j][k] +
						   box->temp_pre_z[i][j-1][k] + box->temp_pre_z[i][j+1][k] +
						   box->temp_pre_z[i][j][k-1] + box->temp_pre_z[i][j][k+1];

					box->pre_x[i][j][k] = (sumx + alpha*box->divergence[i][j][k])/beta;
					box->pre_y[i][j][k] = (sumy + alpha*box->divergence[i][j][k])/beta;
					box->pre_z[i][j][k] = (sumz + alpha*box->divergence[i][j][k])/beta;

				}
			}

		}

	}

}

void accountForGradient(FluidBox *box) {

	for (int i = 1; i < box->depth - 1; i++){

		for (int j = 1; j < box->width - 1; j++) {

			for (int k = 1; j < box->length - 1; k++) {

				box->grad_x[i][j][k] = (box->pre_x[i+1][j][k] - box->pre_x[i-1][j][k])/2;
				box->grad_y[i][j][k] = (box->pre_y[j+1][j][k] - box->pre_y[j-1][j][k])/2;
				box->grad_z[i][j][k] = (box->pre_z[k+1][j][k] - box->pre_z[k-1][j][k])/2;

				box->temp_vel_x[i][j][k] -= box->grad_x[i][j][k];
				box->temp_vel_y[i][j][k] -= box->grad_y[i][j][k];
				box->temp_vel_z[i][j][k] -= box->grad_z[i][j][k];

				box->vel_x[i][j][k] -= box->grad_x[i][j][k];
				box->vel_y[i][j][k] -= box->grad_y[i][j][k];
				box->vel_z[i][j][k] -= box->grad_z[i][j][k];

			}
		}

	}

}

void timeStep(FluidBox *box, int mouseX, int mouseY, int mouseZ, 
	float vel_x, float vel_y, float vel_z){

	advectCube(box);
	diffuseCube(box);
	addForce(box,mouseX,mouseY,mouseZ,vel_x,vel_y,vel_z);
	projectBox(box);
	accountForGradient(box);

}






