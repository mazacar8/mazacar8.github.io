#ifndef __CUDA_RENDERER_H__
#define __CUDA_RENDERER_H__

#ifndef uint
#define uint unsigned int
#endif

#include "nv_seq2d.h"

#define BLOCKDIM 48

struct Image;

typedef enum {
    WATER_CUBE
} SceneName;

class CudaRenderer{

private:

    Image* image;
    SceneName sceneName;

    // int length;
    // int width;

    // float time_step_size;
    // float diff_const;

    // int numParticles;
    // int size;      

    // float** vel_x;
    // float** vel_y;

    // float** temp_vel_x;
    // float** temp_vel_y;

    // float** pre_x;
    // float** pre_y;

    // float** temp_pre_x;
    // float** temp_pre_y;

    // float** grad_x;
    // float** grad_y;

    // float** divergence;
    // bool** particle;

    float* cudaDevice_imageData;   

    float** cudaDevice_vel_x;
    float** cudaDevice_vel_y;

    float** cudaDevice_temp_vel_x;
    float** cudaDevice_temp_vel_y;

    float** cudaDevice_pre_x;
    float** cudaDevice_pre_y;

    float** cudaDevice_temp_pre_x;
    float** cudaDevice_temp_pre_y;

    float** cudaDevice_grad_x;
    float** cudaDevice_grad_y;

    float** cudaDevice_divergence;
    bool** cudaDevice_particle;
    bool** cudaDevice_temp_particle;

    FluidBox *box;

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

    // void shadePixel(
    //     int circleIndex,
    //     float pixelCenterX, float pixelCenterY,
    //     float px, float py, float pz,
    //     float* pixelData);
};


#endif