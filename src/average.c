#include "average.h"

//simply sums up the values of the array and divides by the total count to get the average
//if excluding_halo is equal to 1(DIM_Y), then we will delete the halo values
//start and end points changes by 1 iteration

float average_pixel(int M, int N, float *array, int excluding_halo)
{
    float result, total_pixel = 0;
    int halo_effect;
    int i, j;

    if (excluding_halo == 1) {
        halo_effect = 1;
    } else {
        halo_effect = 0;
    }

    for (i = (0 + halo_effect); i < ( M - halo_effect); i++) {
        for(j= (0 + halo_effect); j < (N - halo_effect); j++) {
            total_pixel += *array;
        }
    }
    result = total_pixel/(M*N);
    return result;
}