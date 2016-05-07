
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

FluidBox *FluidBoxCreate(int length, int width, int depth, float ts);
void FluidBoxFree(FluidBox *box);
void addVelocity(FluidBox *box, int x, int y, int z, float vel_x, float vel_y, float vel_z);
void advectCube(FluidBox *box);
void diffuseCube(FluidBox *box);
void addForce(FluidBox *box, int x, int y, int z,
	float vel_x, float vel_y, float vel_z);
void computeDivergence(FluidBox *box);
void projectBox(FluidBox *box);
void accountForGradient(FluidBox *box);
void timeStep(FluidBox *box);
void copy3dArray(float*** dst,float*** src, int length, int width, int depth);
void setZero3d(float*** array, int length, int width, int depth);



