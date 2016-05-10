#ifndef __CUDA_RENDERER_H__
#define __CUDA_RENDERER_H__

#ifndef uint
#define uint unsigned int
#endif

#include "nv_seq.h"

#define BLOCKDIM 32

struct Image;

typedef enum {
    WATER_CUBE
} SceneName;

class CudaRenderer{

private:

    Image* image;
    SceneName sceneName;

    float* cudaDevice_imageData;   

    float* cudaDevice_vel_x;
    float* cudaDevice_vel_y;

    float* cudaDevice_temp_vel_x;
    float* cudaDevice_temp_vel_y;

    float* cudaDevice_pre_x;
    float* cudaDevice_pre_y;

    float* cudaDevice_temp_pre_x;
    float* cudaDevice_temp_pre_y;

    float* cudaDevice_grad_x;
    float* cudaDevice_grad_y;

    float* cudaDevice_divergence;
    bool* cudaDevice_particle;
    bool* cudaDevice_temp_particle;

    FluidBox1D *box;

public:

    CudaRenderer();
    virtual ~CudaRenderer();

    const Image* getImage();

    void setup();

    void loadScene(SceneName name);

    void allocOutputImage(int width, int height);

    void clearImage();

    void advanceAnimation();

    void render();

    void shadePixel(int i, int j);
};


#endif