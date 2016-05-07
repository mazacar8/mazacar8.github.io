#include <stdio.h>
#include <math.h>
#include <vector>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <string>
#include <cuda.h>
#include <cuda_runtime.h>

#include "cudaRenderer.h"
#include "image.h"
#include "sceneLoader.h"

////////////////////////////////////////////////////////////////////////////////////////
// Putting all the cuda kernels here
////////////////////////////////////////////////////////////////////////////////////


typedef struct {

	SceneName sceneName;

	int length;
	int width;

	float time_step_size;
	float diff_const;

	int numParticles;
	int size;      

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

	int imageWidth;
    int imageHeight;
    float* imageData;

} GlobalConstants;

// Global variable that is in scope, but read-only, for all cuda
// kernels.  The __constant__ modifier designates this variable will
// be stored in special "constant" memory on the GPU. (we didn't talk
// about this type of memory in class, but constant memory is a fast
// place to put read-only variables).
__constant__ GlobalConstants cuConstRendererParams;

__global__ void kernelClearImage(float r, float g, float b, float a) {

    int imageX = blockIdx.x * blockDim.x + threadIdx.x;
    int imageY = blockIdx.y * blockDim.y + threadIdx.y;

    int width = cuConstRendererParams.imageWidth;
    int height = cuConstRendererParams.imageHeight;

    if (imageX >= width || imageY >= height)
        return;

    int offset = 4 * (imageY * width + imageX);
    float4 value = make_float4(r, g, b, a);

    // write to global memory: As an optimization, I use a float4
    // store, that results in more efficient code than if I coded this
    // up as four seperate fp32 stores.
    *(float4*)(&cuConstRendererParams.imageData[offset]) = value;
}

__global__ void kernelAdvanceWaterCube(){


}


////////////////////////////////////////////////////////////////////////////////////////
// Cuda Renderer Class
////////////////////////////////////////////////////////////////////////////////////


CudaRenderer::CudaRenderer() {

	image = NULL;
    box = NULL;
    sceneName = WATER_CUBE;

	// length = 0;
	// width = 0;

	// time_step_size = 0;
	// diff_const = 0;

	// numParticles = 0;
	// size = 0;

	// vel_x = NULL;
	// vel_y = NULL;

	// temp_vel_x = NULL;
	// temp_vel_y = NULL;

	// pre_x = NULL;
	// pre_y = NULL;

	// temp_pre_x = NULL;
	// temp_pre_y = NULL;

	// grad_x = NULL;
	// grad_y = NULL;

	// divergence = NULL;
	// particle = NULL;

	cudaDevice_imageData = NULL;    

	cudaDevice_vel_x = NULL;
	cudaDevice_vel_y = NULL;

	cudaDevice_temp_vel_x = NULL;
	cudaDevice_temp_vel_y = NULL;

	cudaDevice_pre_x = NULL;
	cudaDevice_pre_y = NULL;

	cudaDevice_temp_pre_x = NULL;
	cudaDevice_temp_pre_y = NULL;

	cudaDevice_grad_x = NULL;
	cudaDevice_grad_y = NULL;

	cudaDevice_divergence = NULL;
	cudaDevice_particle = NULL;
    cudaDevice_temp_particle = NULL;

}

CudaRenderer::~CudaRenderer(){

	if(image) {
		delete image;
	}

	delete [] box->vel_x;
	delete [] box->vel_y;
	delete [] box->temp_vel_x;
	delete [] box->temp_vel_y;
	delete [] box->pre_x;
	delete [] box->pre_y;
	delete [] box->temp_pre_x;
	delete [] box->temp_pre_y;
	delete [] box->particle;
    delete [] box->temp_particle;
	delete [] box->divergence;
	delete [] box->grad_x;
	delete [] box->grad_y;
    delete [] box;

	cudaFree(cudaDevice_vel_x);
	cudaFree(cudaDevice_vel_y);
	cudaFree(cudaDevice_temp_vel_x);
	cudaFree(cudaDevice_temp_vel_y);
	cudaFree(cudaDevice_pre_x);
	cudaFree(cudaDevice_pre_y);
	cudaFree(cudaDevice_temp_pre_x);
	cudaFree(cudaDevice_temp_pre_y);
	cudaFree(cudaDevice_particle);
    cudaFree(cudaDevice_temp_particle);
	cudaFree(cudaDevice_divergence);
	cudaFree(cudaDevice_grad_x);
	cudaFree(cudaDevice_grad_y);
	cudaFree(cudaDevice_imageData);
}

// allocOutputImage --
//
// Allocate buffer the renderer will render into.  Check status of
// image first to avoid memory leak.
void
CudaRenderer::allocOutputImage(int width, int height) {

    if (image)
        delete image;
    image = new Image(width, height);
}

const Image*
CudaRenderer::getImage() {

    // need to copy contents of the rendered image from device memory
    // before we expose the Image object to the caller

    printf("Copying image data from device\n");

    cudaMemcpy(image->data,
               cudaDevice_imageData,
               sizeof(float) * 4 * image->width * image->height,
               cudaMemcpyDeviceToHost);

    return image;
}

// clearImage --
//
// Clear's the renderer's target image.  The state of the image after
// the clear depends on the scene being rendered.
void
CudaRenderer::clearImage() {

    // 256 threads per block is a healthy number
    dim3 blockDim(16, 16, 1);
    dim3 gridDim(
        (image->width + blockDim.x - 1) / blockDim.x,
        (image->height + blockDim.y - 1) / blockDim.y);

    kernelClearImage<<<gridDim, blockDim>>>(1.f, 1.f, 1.f, 1.f);
   
    cudaThreadSynchronize();
}

void
CudaRenderer::loadScene(SceneName scene) {
    sceneName = scene;
    // loadWaterScene(length, width, time_step_size, diff_const, numParticles, 
    //     size, sceneName, vel_x, vel_y, temp_vel_x, temp_vel_y, pre_x, pre_y, 
    //     temp_pre_x, temp_pre_y, particle, divergence, grad_x, grad_y);

    box = loadWaterScene();

}

void
CudaRenderer::setup() {

    int deviceCount = 0;
    bool isFastGPU = false;
    std::string name;
    cudaError_t err = cudaGetDeviceCount(&deviceCount);

    printf("---------------------------------------------------------\n");
    printf("Initializing CUDA for CudaRenderer\n");
    printf("Found %d CUDA devices\n", deviceCount);

    for (int i=0; i<deviceCount; i++) {
        cudaDeviceProp deviceProps;
        cudaGetDeviceProperties(&deviceProps, i);
        name = deviceProps.name;
        if (name.compare("GeForce GTX 780") == 0)
        {
            isFastGPU = true;
        }

        printf("Device %d: %s\n", i, deviceProps.name);
        printf("   SMs:        %d\n", deviceProps.multiProcessorCount);
        printf("   Global mem: %.0f MB\n", static_cast<float>(deviceProps.totalGlobalMem) / (1024 * 1024));
        printf("   CUDA Cap:   %d.%d\n", deviceProps.major, deviceProps.minor);
    }
    printf("---------------------------------------------------------\n");
    if (!isFastGPU)
    {
        printf("WARNING: "
               "You're not running on a fast GPU, please consider using "
               "NVIDIA GTX 780.\n");
        printf("---------------------------------------------------------\n");
    }
    
    // By this time the scene should be loaded.  Now copy all the key
    // data structures into device memory so they are accessible to
    // CUDA kernels
    //
    // See the CUDA Programmer's Guide for descriptions of
    // cudaMalloc and cudaMemcpy

    int length = box->length;
    int width = box->width;

    cudaMalloc(&cudaDevice_vel_x, sizeof(float) * length * width);
    cudaMalloc(&cudaDevice_vel_y, sizeof(float) * length * width);
    cudaMalloc(&cudaDevice_temp_vel_x, sizeof(float) * length * width);
    cudaMalloc(&cudaDevice_temp_vel_y, sizeof(float) * length * width);
    cudaMalloc(&cudaDevice_pre_x, sizeof(float) * length * width);
    cudaMalloc(&cudaDevice_pre_y, sizeof(float) * length * width);
    cudaMalloc(&cudaDevice_temp_pre_x, sizeof(float) * length * width);
    cudaMalloc(&cudaDevice_temp_pre_y, sizeof(float) * length * width);
    cudaMalloc(&cudaDevice_grad_x, sizeof(float) * length * width);
    cudaMalloc(&cudaDevice_grad_y, sizeof(float) * length * width);
    cudaMalloc(&cudaDevice_particle, sizeof(float) * length * width);
    cudaMalloc(&cudaDevice_temp_particle, sizeof(float) * length * width);
    cudaMalloc(&cudaDevice_divergence, sizeof(float) * length * width);
    cudaMalloc(&cudaDevice_imageData, sizeof(float) * 4 * length * width);

    cudaMemcpy(cudaDevice_vel_x, box->vel_x, sizeof(float) * length * width, cudaMemcpyHostToDevice);
    cudaMemcpy(cudaDevice_vel_y, box->vel_y, sizeof(float) * length * width, cudaMemcpyHostToDevice);
    cudaMemcpy(cudaDevice_temp_vel_x, box->temp_vel_x, sizeof(float) * length * width, cudaMemcpyHostToDevice);
    cudaMemcpy(cudaDevice_temp_vel_y, box->temp_vel_y, sizeof(float) * length * width, cudaMemcpyHostToDevice);
    cudaMemcpy(cudaDevice_pre_x, box->pre_x, sizeof(float) * length * width, cudaMemcpyHostToDevice);
    cudaMemcpy(cudaDevice_pre_y, box->pre_x, sizeof(float) * length * width, cudaMemcpyHostToDevice);
    cudaMemcpy(cudaDevice_temp_pre_x, box->temp_pre_x, sizeof(float) * length * width, cudaMemcpyHostToDevice);
    cudaMemcpy(cudaDevice_temp_pre_y, box->temp_pre_y, sizeof(float) * length * width, cudaMemcpyHostToDevice);
    cudaMemcpy(cudaDevice_grad_x, box->grad_x, sizeof(float) * length * width, cudaMemcpyHostToDevice);
    cudaMemcpy(cudaDevice_grad_y, box->grad_y, sizeof(float) * length * width, cudaMemcpyHostToDevice);
    cudaMemcpy(cudaDevice_particle, box->particle, sizeof(float) * length * width, cudaMemcpyHostToDevice);
    cudaMemcpy(cudaDevice_temp_particle, box->temp_particle, sizeof(float) * length * width, cudaMemcpyHostToDevice);
    cudaMemcpy(cudaDevice_divergence, box->divergence, sizeof(float) * length * width, cudaMemcpyHostToDevice);
    

    // Initialize parameters in constant memory.  We didn't talk about
    // constant memory in class, but the use of read-only constant
    // memory here is an optimization over just sticking these values
    // in device global memory.  NVIDIA GPUs have a few special tricks
    // for optimizing access to constant memory.  Using global memory
    // here would have worked just as well.  See the Programmer's
    // Guide for more information about constant memory.

    GlobalConstants params;

    params.sceneName = sceneName;

	params.length = length;
	params.width = width;

	params.time_step_size = box->time_step_size;
	params.diff_const = box->diff_const;

	params.numParticles = box->numParticles;
	params.size = box->size;      

	params.vel_x = cudaDevice_vel_x;
	params.vel_y = cudaDevice_vel_y;

	params.temp_vel_x = cudaDevice_temp_vel_x;
	params.temp_vel_y = cudaDevice_temp_vel_y;

	params.pre_x = cudaDevice_pre_x;
	params.pre_y = cudaDevice_pre_y;

	params.temp_pre_x = cudaDevice_temp_pre_x;
	params.temp_pre_y = cudaDevice_temp_pre_y;

	params.grad_x = cudaDevice_grad_x;
	params.grad_y = cudaDevice_grad_y;

	params.divergence = cudaDevice_divergence;
	params.particle = cudaDevice_particle;
    params.temp_particle = cudaDevice_temp_particle;

	params.imageHeight = image->height;
	params.imageWidth = image->width;
	params.imageData = cudaDevice_imageData;

    cudaMemcpyToSymbol(cuConstRendererParams, &params, sizeof(GlobalConstants));
}


void
CudaRenderer::advanceAnimation() {
     // 256 threads per block is a healthy number
    dim3 blockDim(256, 1);
    dim3 gridDim((42 + blockDim.x - 1) / blockDim.x);

    if(sceneName == WATER_CUBE) {
    	kernelAdvanceWaterCube<<<gridDim, blockDim>>>();
    }
    cudaThreadSynchronize();
}

void
CudaRenderer::render() {

    //using a 2d block
    dim3 blockDim(BLOCKDIM,BLOCKDIM);
    //splitting the image into a 2d grid of 2d blocks
    dim3 gridDim(
        (image->width + blockDim.x - 1) / blockDim.x,
        (image->height + blockDim.y - 1) / blockDim.y);
    printf("In Host\n");

    //kernelRender<<<gridDim, blockDim>>>();
    cudaThreadSynchronize();

}