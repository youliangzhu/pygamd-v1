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
#include <stdio.h>

void checkCUDAError(const char *msg)
	{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) 
		{
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(-1);
		}                         
	}

	
__device__ bool exclusion_check(unsigned int* d_n_exclusion, 
									unsigned int* d_exclusion_list, 
									unsigned int idx, 
									unsigned int idx_neighbor,
									unsigned int N)
	{
	unsigned int n_ex = d_n_exclusion[idx];
	for (unsigned int i=0; i< n_ex; i++)
		{
		unsigned int idx_ex = d_exclusion_list[idx+N*i];
		if(idx_ex == idx_neighbor)
			return true;
		}
	return false;
	}
	
__global__ void gpu_compute_rdf_kernel(float4* dpos,
                                     unsigned int N,
									 unsigned int N_total,
                                     float Lx,
									 float Ly,
									 float Lz,
									 float Lxinv,
									 float Lyinv,
									 float Lzinv,
									 float delr,
									 unsigned int *gg,
									 unsigned int maxbin,
									 unsigned int* d_group,
									 unsigned int* d_n_exclusion,
									 unsigned int* d_exclusion_list,
									 unsigned int* d_mol_id_per_particle,
									 bool exclusion_mol,
									 bool exclusion_list,
									 bool bytype,
									 unsigned int block_size)
    {
    // shared data to store all of the particles we compare against
    extern __shared__ float sdata[];
   
    // load in the particle
    int pidx = blockIdx.x * block_size + threadIdx.x;

    float4 pos = make_float4(0.0f,0.0f,0.0f,0.0f);
	unsigned int molid = 0;
	unsigned int tag = pidx;
    if (pidx < N)
		{
		tag = d_group[pidx];
		pos = dpos[pidx];
		if(exclusion_mol)
			molid = d_mol_id_per_particle[tag];
		}
	unsigned int typi = (unsigned int ) pos.w;
	
    for (int start = 0; start < N; start += block_size)
        {
        // load data
			float4 pos2 = make_float4(0.0f,0.0f,0.0f,0.0f);
			if (start + threadIdx.x < N)
				pos2 = dpos[start + threadIdx.x];

        __syncthreads();
        sdata[threadIdx.x] = pos2.x;
        sdata[threadIdx.x + block_size] = pos2.y;
        sdata[threadIdx.x + 2*block_size] = pos2.z;
        sdata[threadIdx.x + 3*block_size] = pos2.w; //< unused, but try to get compiler to fully coalesce reads

        __syncthreads();
        

        int end_offset= block_size;
        end_offset = min(end_offset, N - start);
        
        if (pidx < N)
            {
            for (int cur_offset = 0; cur_offset < end_offset; cur_offset++)
                {
                // calculate dr
				float dx = pos.x - sdata[cur_offset];
				float dy = pos.y - sdata[cur_offset + block_size];
				float dz = pos.z - sdata[cur_offset + 2*block_size];				
				unsigned int typj = (unsigned int ) sdata[cur_offset + 3*block_size];
				unsigned int tag_neighbor = d_group[start + cur_offset];
				if(!bytype||(typi!=typj))
					{
					dx -= Lx*rintf(dx*Lxinv);
					dy -= Ly*rintf(dy*Lyinv);
					dz -= Lz*rintf(dz*Lzinv);				
					float rijSQ = dx*dx + dy*dy + dz*dz;
					float rij = sqrt(rijSQ);
					int bin = int(rij/delr);
					unsigned int idx = bin+pidx*maxbin;
					bool exclusion = false;
					if(exclusion_mol)
						{
						unsigned int molid_neighbor = d_mol_id_per_particle[tag_neighbor];
						if(molid==molid_neighbor)
							exclusion = true;
						}
					if(exclusion_list)
						exclusion = exclusion_check(d_n_exclusion, d_exclusion_list, tag, tag_neighbor, N_total);
	//				printf("thread %d, %d\n", pidx, bin);	
					if ( (bin < maxbin) && ((start + cur_offset) != pidx) && !exclusion)
						gg[idx] += 1;
					}
                }
            }
        }

    }


//! Shared memory used in reducing the sums



__global__ void gpu_compute_partial_sums_share_memory(unsigned int *gg,
													unsigned int *d_scratch,
													unsigned int maxbin)
    {
    extern __shared__ unsigned int compute_sdata[];	
	
    // determine which particle this thread works on
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
	for(unsigned int i=0;i<maxbin;i++)
		compute_sdata[threadIdx.x*maxbin + i] = gg[idx*maxbin + i];
    __syncthreads();
    
    // reduce the sum in parallel
    int offs = blockDim.x >> 1;
    while (offs > 0)
        {
        if (threadIdx.x < offs)
            {
			for(unsigned int i=0;i<maxbin;i++)
				compute_sdata[threadIdx.x*maxbin + i] += compute_sdata[(threadIdx.x+offs)*maxbin + i] ;
            }
        offs >>= 1;
        __syncthreads();
        }
        
    // write out our partial sum
    if (threadIdx.x == 0)
        {
		for(unsigned int i=0;i<maxbin;i++)
			d_scratch[blockIdx.x*maxbin + i] += compute_sdata[i];
//		printf("thread %d, %d\n", idx, d_scratch[blockIdx.x*maxbin + 0]);					
        }
    }

__global__ void gpu_compute_partial_sums(unsigned int *gg,
										unsigned int *d_scratch,
										unsigned int maxbin)
    {
    // determine which particle this thread works on
    int idx = blockIdx.x * blockDim.x + threadIdx.x;

    // reduce the sum in parallel
    int offs = blockDim.x >> 1;
    while (offs > 0)
        {
        if (threadIdx.x < offs)
            {
			for(unsigned int i=0;i<maxbin;i++)
				gg[idx*maxbin + i] += gg[(idx+offs)*maxbin + i] ;
            }
        offs >>= 1;
        __syncthreads();
        }
        
    // write out our partial sum
    if (threadIdx.x == 0)
        {
		for(unsigned int i=0;i<maxbin;i++)
			d_scratch[blockIdx.x*maxbin + i] += gg[idx*maxbin + i];
//		printf("thread %d, %d\n", idx, d_scratch[blockIdx.x*maxbin + 0]);					
        }
    }


extern "C" cudaError_t gpu_compute_rdf(float4* dpos,
									unsigned int N,
									unsigned int N_total,							
									float Lx,
									float Ly,
									float Lz,
									float Lxinv,
									float Lyinv,
									float Lzinv,
									float delr,
									unsigned int *gg,
									unsigned int *d_scratch,
									unsigned int maxbin,
									unsigned int* d_group,
									unsigned int* d_n_exclusion,
									unsigned int* d_exclusion_list,
									unsigned int* d_mol_id_per_particle,
									bool exclusion_mol,
									bool exclusion_list,							
									bool bytype,							
									unsigned int block_size)
    {

    // setup the grid to run the kernel
    dim3 grid( (int)ceil((float)N / (float)block_size), 1, 1);
    dim3 threads(block_size, 1, 1);	
	

		
	unsigned int shared_bytes = sizeof(float)*block_size*4;
	gpu_compute_rdf_kernel<<<grid,threads, shared_bytes>>>(dpos,
															N,
															N_total,
															Lx,
															Ly,
															Lz,
															Lxinv, 
															Lyinv,
															Lzinv,
															delr,
															gg,
															maxbin,
															d_group,
															d_n_exclusion,
															d_exclusion_list,
															d_mol_id_per_particle,
															exclusion_mol,
															exclusion_list,
															bytype,
															block_size); 
    cudaDeviceSynchronize();   
	checkCUDAError("gpu_compute_rdf_kernel");
	
    shared_bytes = sizeof(unsigned int)*block_size*maxbin;
	if(shared_bytes<48000)
		{	
		gpu_compute_partial_sums_share_memory<<<grid, threads, shared_bytes>>>(gg,
																				d_scratch,
																				maxbin);
		 cudaDeviceSynchronize();
		}
	else
		{
		gpu_compute_partial_sums<<<grid, threads>>>(gg,
													d_scratch,
													maxbin);															
		 cudaDeviceSynchronize();   
		}
	checkCUDAError("gpu_compute_partial_sums");    

    
    return cudaSuccess;
    }


