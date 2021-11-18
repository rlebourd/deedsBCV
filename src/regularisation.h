#ifndef REGULARISATION_H_INCLUDED
#define REGULARISATION_H_INCLUDED

/* Incremental diffusion regularisation of parametrised transformation
 using (globally optimal) belief-propagation on minimum spanning tree.
 Fast distance transform uses squared differences.
 Similarity cost for each node and label has to be given as input.
*/
void messageDT(int ind,float* data,short* indout,int len1,float offsetx,float offsety,float offsetz);

void regularisationCL(float* costall,float* u0,float* v0,float* w0,float* u1,float* v1,float* w1,int hw,int step1,float quant,int* ordered,int* parents,float* edgemst);

#endif // REGULARISATION_H_INCLUDED
