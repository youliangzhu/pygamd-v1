/*
PYGAMD - Python GPU-Accelerated Molecular Dynamics Software
VERSION 1
COPYRIGHT
	PYGAMD Copyright (c) (2021) You-Liang Zhu, Zhong-Yuan Lu
LICENSE
	This program is a free software: you can redistribute it and/or 
	modify it under the terms of the GNU General Public License. 
	This program is distributed in the hope that it will be useful, 
	but WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANT ABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
	See the General Public License v3 for more details.
	You should have received a copy of the GNU General Public License
	along with this program. If not, see <http://www.gnu.org/licenses/>.
DISCLAIMER
	The authors of PYGAMD do not guarantee that this program and its 
	derivatives are free from error. In no event shall the copyright 
	holder or contributors be liable for any indirect, incidental, 
	special, exemplary, or consequential loss or damage that results 
	from its use. We also have no responsibility for providing the 
	service of functional extension of this program to general users.
USER OBLIGATION 
	If any results obtained with PYGAMD are published in the scientific 
	literature, the users have an obligation to distribute this program 
	and acknowledge our efforts by citing the paper "Y.-L. Zhu, H. Liu, 
	Z.-W. Li, H.-J. Qian, G. Milano, and Z.-Y. Lu, J. Comput. Chem. 2013,
	34, 2197-2211" in their article.
CORRESPONDENCE
	Dr. You-Liang Zhu
	Email: ylzhu@pygamd.com
*/

#include <cuda_runtime.h>
#include "cuComplex.h"
#include <stdio.h>

void checkCUDAErrors(const char *msg)
	{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
		{
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(-1);
		}                         
	}



#define max_Ntypes 10
__global__ void gpu_compute_domainsize_kernel(float *x,
									 float *y,
                                     float *z,
									 float *w,
                                     unsigned int Num,
									 int Lmax,
									 int Mmax,
									 int Nmax,
									 float Lxinv,
									 float Lyinv,
									 float Lzinv,
									 int Kmax,
									 float TWOPI,
									 unsigned int Ntypes,
									 cuFloatComplex *ss,
									 unsigned int block_size)
    {
    // shared data to store all of the particles we compare against
    extern __shared__ float sdata[];
   
    // load in the particle
    int idx = blockIdx.x * block_size + threadIdx.x;
	
	int n = int (float(idx)/float(Lmax*Mmax));
	int m = int (float(idx - n*Lmax*Mmax)/float(Lmax));
	int l = idx -n*Lmax*Mmax-m*Lmax;
	
	n -= Kmax;	
	m -= Kmax;

//	int kk = l*l + m*m + n*n;


	cuFloatComplex sums[max_Ntypes];
	
	for(unsigned int i =0;i < Ntypes; i++)
		{
		sums[i].x = 0.0;
		sums[i].y = 0.0;
		}


    // track the number of neighbors added so far


    // each block is going to loop over all Num particles (this assumes memory is padded to a multiple of blockDim.x)
    // in blocks of blockDim.x
    for (int start = 0; start < Num; start += block_size)
        {
        // load data
		float px = 0.0;
		float py = 0.0;
		float pz = 0.0;
		float pw = 0.0;
		if (start + threadIdx.x < Num)
			{
			px = x[start + threadIdx.x];
			py = y[start + threadIdx.x];
			pz = z[start + threadIdx.x];
			pw = w[start + threadIdx.x];
			}
        // make sure everybody is caught up before we stomp on the memory
        __syncthreads();
        sdata[threadIdx.x] = px;
        sdata[threadIdx.x + block_size] = py;
        sdata[threadIdx.x + 2*block_size] = pz;
        sdata[threadIdx.x + 3*block_size] = pw; //< unused, but try to get compiler to fully coalesce reads
        
        // ensure all data is loaded
        __syncthreads();
        
        // now each thread loops over every particle in shmem, but doesn't loop past the end of the particle list (since
        // the block might extend that far)
        int end_offset= block_size;
        end_offset = min(end_offset, Num - start);

        for (int cur_offset = 0; cur_offset < end_offset; cur_offset++)
            {
                // calculate dr
			float p2x = sdata[cur_offset];
			float p2y = sdata[cur_offset + block_size];
			float p2z = sdata[cur_offset + 2*block_size];
			float p2w = sdata[cur_offset + 3*block_size];
			int typi = __float_as_int(p2w);	
			float theta = float(l)*p2x*Lxinv + float(m)*p2y*Lyinv + float(n)*p2z*Lzinv;
			theta *= TWOPI;
			cuFloatComplex expikr;
			expikr.x = cosf(theta);
			expikr.y = sinf(theta);
			sums[typi] = cuCaddf(sums[typi],expikr);
//			if(m==0&&n==0&&typi==0&&l==48)
//				printf("thread %d, %f, %f, %f\n", idx, theta, expikr.x, expikr.y);			
            }
        }
		if(idx < Lmax*Mmax*Nmax)
			{
			for(unsigned int i =0; i<Ntypes;i++)
				{
				ss[i+idx*Ntypes] =sums[i];
	//			printf("thread %d, %d, %f, %f\n", idx, i, sums[i].x,sums[i].y);		
				}
			}
    }



extern "C" cudaError_t gpu_compute_domainsize(float *x,
							float *y,
                            float *z,
							float *w,
                            unsigned int Num,
							int Lmax,
							int Mmax,
							int Nmax,							
							float Lxinv,
							float Lyinv,
							float Lzinv,
							int Kmax,
							float TWOPI,							
							unsigned int Ntypes,
							cuFloatComplex *ss,
							unsigned int block_size)
    {

    // setup the grid to run the kernel
    dim3 grid( (int)ceil((float)(Lmax*Mmax*Nmax) / (float)block_size), 1, 1);
    dim3 threads(block_size, 1, 1);	
	
	unsigned int shared_bytes = sizeof(float)*block_size*4;
	gpu_compute_domainsize_kernel<<<grid,threads, shared_bytes>>>(x,y,z,w,Num,Lmax,Mmax,Nmax,Lxinv,Lyinv,Lzinv,Kmax,TWOPI,Ntypes,ss,block_size); 
 

     cudaDeviceSynchronize();   
    
	checkCUDAErrors("gpu_compute_domainsize_kernel"); 
    
    return cudaSuccess;
    }


