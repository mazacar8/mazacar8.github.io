#include <stdio.h>
#include <math.h>
#define FILL_LEVEL 0.75

typedef struct {

	int length;
	int width;
	int depth;

	float time_step_size;

	int numParticles;
	int size;

	float *vel_x;
	float *vel_y;
	float *vel_z;

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

	box->vel_x = new float [box->size]();
	box->vel_y = new float [box->size]();
	box->vel_z = new float [box->size]();	
}

void FluidBoxFree(FluidBox *box) {

	delete [] box->vel_x;
	delete [] box->vel_y;
	delete [] box->vel_z;
	delete [] box->particle;
	delete box;
}

// void addDensity(FluidBox *box, int x, int y, int z, float density) {

// 	int index = x + (box->length * y) + (box->length * box->width * z);
// 	box->density[index] += density;
// }

// void addVelocity(FluidBox *box, int x, int y, int z, float vel_x, float vel_y, float vel_z) {

// 	int index = x + (box->length * y) + (box->length * box->width * z);
// 	box->vel_x[index] = vel_x;
// 	box->vel_y[index] = vel_y;
// 	box->vel_z[index] = vel_z;
// }

void advectCube(FluidBox *box) {

	float dt = box->time_step_size;
	int x,y,z,old_x,old_y,old_z,index;

	int length = box->length;
	int width = box->width;

	for(int i = 0; i < box->size; i += 1) {

		x = i % length;
		y = (i % (length*width))/length;
		z =  i/(length*width);

		old_x = round(x - box->vel_x[i] * dt);
		old_y = round(y - box->vel_y[i] * dt);
		old_z = round(z - box->vel_z[i] * dt);

		index = old_x + (length * old_y) + (length * width * old_z);
		box->particle[i] = box->particle[index];

		if(box->particle[i]){

			box->vel_x[i] = box->vel_x[index];
		}
		
	}
}

void timeStep(FluidBox *box){

	int length = box->length;
	int depth = box->depth;
	int width = box->width;

	advectCube(box);
}






