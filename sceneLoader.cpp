#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include <functional>

#include "sceneLoader.h"

FluidBox *loadWaterScene(){

	return FluidBoxCreate2D(LENGTH,WIDTH,DT);
}