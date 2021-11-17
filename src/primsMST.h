/* Image-driven minimum-spanning-tree calcuation using Prim's algorithm.
 Average run-time should be of n*log(n) complexity.
 Uses heap data structure to speed-up finding the next lowest edge-weight.
 Edge-weights 
*/

#ifndef PRIMS_MST_H_INCLUDED
#define PRIMS_MST_H_INCLUDED

void primsGraph(float* im1, int* ordered, int* parents, float* edgemst, int step1, int m2, int n2, int o2);

#endif // PRIMS_MST_H_INCLUDED
