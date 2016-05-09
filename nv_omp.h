#include "nv_seq2d.h"


FluidBox *FluidBoxCreateOMP(int length, int width, float ts);
void FluidBoxFreeOMP(FluidBox *box);
void addVelocityOMP(FluidBox *box, int x, int y, float vel_x, float vel_y);
void advectCubeOMP(FluidBox *box);
void diffuseCube2D(FluidBox *box);
void addForceOMP(FluidBox *box);
void computeDivergenceOMP(FluidBox *box);
void projectBoxOMP(FluidBox *box);
void setZeroOMP(float** array, int length, int width);
void copy2dArrayOMP(float** dst,float** src, int length, int width);
void accountForGradientOMP(FluidBox *box);
void timeStepOMP(FluidBox *box);