#include "boundary.h"

//Setting up fixed boundary conditions for top & buttom
//returns sawtooth value of input boundary condition

float boundaryval(int i, int m)
{
    float val;

    val = 2.0*((float)(i-1))/((float)(m-1));
    if (i >= m/2+1) 
        val = 2.0-val;

    return val;
}