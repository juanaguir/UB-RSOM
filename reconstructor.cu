
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <windows.h>
#include <helper_cuda.h>
#include <helper_string.h>
#include <time.h>

__global__ void pBP(float *recc, int *lineind, float *datacq, int * guid, int lz, int ly, int Sx,  int blockx, int blocky, int nblockx, int nblocky){
 

int jj = blockIdx.x;
int j = threadIdx.x;

		 jj=jj+blocky*nblocky;
		 j=j+blockx*nblockx;
		
			for(int i =0; i < (lz-1); i++ ) {
				for(int ii =  guid[i]-1 ; ii <  guid[i+1] ; ii++ ) {				
					recc[(j*(lz-1)*ly)+(jj*(lz-1))+i] = recc[(j*(lz-1)*ly)+(jj*(lz-1))+i] + datacq [int( lineind[ii] ) - 1 + j + jj * Sx ];
	               
																   }
				  		
			
		
}

 

#define N 512
void main( void ) {
	cudaDeviceReset();  

	// managing and gathering host memory //
FILE *fpunt;	// pointer for write and read from harddrive
	 double *Nsx; // Number of elements from the acquired data
	 double *datacq; //sinogram
	 float *datacqf; //sinogram in float type precision (less memory better for GPU, althoung oduble is not that bad...). Could come like this from MATLAB. I KNOW.
	 double *lineind; // addresses of the voxels to be summed in a line of voxels following the trandsucer axis
	 int *lineindint; // addresses of the voxels to be summed in int type (less memory)
	 double *guid;    // guide to read the voxels addresses.
	 int *guidint;
	 
	 int nblockx; int nblocky; // number of blocks in x direction, y direction.  
	 int blockx; int blocky;  // block being reconstructed.
	 int lbx; int lby;  // Number of lines per blocks (number of parralell threads= lbx*lby)  
	float *recc; // reconstruction
	 double *lxyz; //vector containing the size of the array of addresses and the size of the acquired data
	 int Sx; int Sy; int Sz; // Size of the acquired data 
	 int lx; int ly; int lz; // number of veoxels in the reconstuction grid
	
	 int lindex; // lenght of the array containing the memory addresses of the voxels to be sum,
	 int *A; // variable to guide the parallel for in the cuda function
	 Nsx = (double *)malloc(sizeof(double) * 3);
	lxyz= (double *)malloc(sizeof(double) * 4);
	int zlim; // limit for each iteration of the reconstructor 
	int dum; // "helper" variable
	double dumm; // "helper" variable
	dum=3;
    
	 printf("Loading data and copying it to the GPU memory\n");

	//               Gathering the size of the acquired data    //
	fpunt = fopen("path", "rb");
    fread(Nsx, sizeof(double),dum, fpunt); 
    fclose(fpunt);
	

     Sx=int(Nsx[0]); Sy=int(Nsx[1]); Sz=int(Nsx[2]);
	 

	datacq = (double *)calloc(sizeof(double),Sx*Sy*Sz);
	 fpunt = fopen("path", "rb");
    fread(datacq, sizeof(double),Sx*Sy*Sz, fpunt); 
    fclose(fpunt);
	 datacqf = (float *)calloc(sizeof(float),Sx*Sy*Sz);
	 for(int i=0; i<Sx*Sy*Sz; i++){
	 datacqf[i]=float(datacq[i]);
	 }
	  free(datacq);
	 //Gathering the size of the vector containing the memory addresses and acquired data without zeros 
	
	dum=4;
	fpunt = fopen("path", "rb");
    fread(lxyz, sizeof(double),dum, fpunt); 
    fclose(fpunt);
	 
    ly=int(lxyz[0]); lx=int(lxyz[1]); lz=int(lxyz[2]); lindex=int(lxyz[3]);
   	 
	//Gathering the vector of addresses

	lineind = (double *)malloc(sizeof(double)*lindex); //acquired data
	fpunt = fopen("path", "rb");
    fread(lineind, sizeof(double),lindex, fpunt); 
    fclose(fpunt);
	 lineindint = (int *)calloc(sizeof(int),lindex);
	 for(int i=0; i<lindex; i++){
	 lineindint[i]=int(lineind[i]);
	 }
	 free(lineind);

    //Gathering the guide of the vector of addresses
    guid = (double*)malloc(sizeof(double)*lz); //acquired data
	fpunt = fopen("path", "rb");
     fread(guid, sizeof(double),lz, fpunt); 
    fclose(fpunt);
	 guidint = (int*)malloc(sizeof(int)*lz); //acquired data
	for(int i=0;i<lz;i++){
		guidint[i]=int(guid[i]);
	}
	 free(guid);
     // allocating memory for the reconstruction
	recc = ( float* ) calloc( sizeof( float), ( lz-1 ) * lx * ly );

	
// changing mclasses to svae memmoery 
	

// managing and gathering device memory//
int *devA; 
float *devrecc; 
int *devlineind; 
float *devdatacq; 
int *devguid;  

cudaMalloc( (void**)&devrecc, sizeof(float)* ( lz-1 ) * lx * ly);
cudaMemcpy( devrecc, recc,sizeof(float)* ( lz-1 ) * lx * ly, cudaMemcpyHostToDevice);
cudaMalloc( (void**)&devlineind, sizeof(int)*lindex );
cudaMemcpy( devlineind, lineindint,sizeof(int)*lindex, cudaMemcpyHostToDevice);
cudaMalloc( (void**)&devdatacq, sizeof(float)*Sx*Sy*Sz);
cudaMemcpy( devdatacq, datacqf,sizeof(float)*Sx*Sy*Sz, cudaMemcpyHostToDevice);
cudaMalloc( (void**)&devguid, sizeof(int)*lz);
cudaMemcpy( devguid, guidint,sizeof(int)*lz, cudaMemcpyHostToDevice);


//************************ setting the number of lines to be reconstructed in parallel (lbx*lby) *************************//
lbx=50;  
lby=50;  
//***********************************************************************************************************************//

nblockx= int(floor(double(lx)/double(lbx))); // number of blocks in x direction
nblocky= int(floor(double(ly)/double(lby))); // number of blocks in y direction


dim3 blocks(lbx, 1); // conditioning cuda memory blocks
dim3 grids(lby, 1);  // 

//dim3 blocks(lx, 1); // conditioning cuda meory blocks
//dim3 grids(ly, 1);

int dumblock=0; // dummy variable to count the numer of blocks already reconstructed
printf("Reconstructing\n");
//for(blockx=0; blockx<nblockx; blockx++) {
blockx=0;
blocky=0;
clock_t tic = clock();
for(blockx=0; blockx<nblockx; blockx++) {
	for(blocky=0;blocky<nblocky; blocky++) {
	pBP<<< grids,blocks >>>(devrecc, devlineind, devdatacq, devguid, lz,ly, Sx, blockx,blocky,lbx,lby);
    cudaDeviceSynchronize();
    dumblock+=1;
    printf("block %d out of %d reconstructed \n", dumblock,nblockx*nblocky);
	}
}
clock_t toc = clock();
 printf("Elapsed: %f seconds\n", (double)(toc - tic) / CLOCKS_PER_SEC);

 printf("reconstruction finished, saving data\n");
cudaMemcpy( recc, devrecc, sizeof(float)* ( lz-1 ) * lx * ly, cudaMemcpyDeviceToHost);
cudaFree( devrecc);
cudaFree( devlineind);
cudaFree( devdatacq);
cudaFree( devguid);


double* reccd = ( double* ) calloc( sizeof( double ), ( lz-1 ) * lx * ly );
for(int i=0; i<( lz-1 ) * lx * ly ;i++){
reccd[i]=double(recc[i]);
}

fpunt = fopen("path","w b");

if (fpunt == NULL)
{
printf("The file did not open");
}

fwrite (reccd, sizeof(double),lx*ly*(lz-1),fpunt);
fclose(fpunt);
 printf("Fisnished! press any key");
 getchar();

}