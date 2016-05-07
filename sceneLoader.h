#ifndef __SCENE_LOADER_H__
#define __SCENE_LOADER_H__

#include "cudaRenderer.h"

void
loadWaterScene(int& length, 
				int& width, 
				float& time_step_size, 
				float& diff_const, 
				int& numParticles, 
				int& size, 
				SceneName sceneName, 
	            float**& vel_x, 
	            float**& vel_y, 
	            float**& temp_vel_x, 
	            float**& temp_vel_y, 
	            float**& pre_x, 
	            float**& pre_y, 
	            float**& temp_pre_x, 
	            float**& temp_pre_y, 
	            bool**& particle, 
	            float**& divergence, 
	            float**& grad_x, 
	            float**& grad_y);

FluidBox *loadWaterScene();

#endif
