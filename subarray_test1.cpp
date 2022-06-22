#include <mpi.h>
#include <stdio.h>
#include <vector>
#include <iomanip>
#include "block2d.h"
using namespace std;


/*         
				block0                                                 block1
    	 ^ +-----+-----+-----+            +-----+-----+           +-----+-----+-----+                        
 *       | |  3  |  7  |  11 | subarray   |  7  |  11 |  paste    |     |     |     |                       
 *       | +-----+-----+-----+ ------->   +-----+-----+  -------> +-----+-----+-----+                        
 *       | |  2  |  6  |  10 |            |  6  |  10 |  anywhere |     |     |     |                       
 *   sizeY +-----+-----+-----+            +-----+-----+           +-----+-----+-----+                       
 *       | |  1  |  5  |  9  |                                    |     |     |     |      
 *       | +-----+-----+-----+                                    +-----+-----+-----+      
 *       | |  0  |  4  |  8  |                                    |     |     |     |     
 *         +-----+-----+-----+                                    +-----+-----+-----+      
 *        +------ sizeX ----->                        
 *                                                      
 *                                                      
*/


int main(){

	MPI_Init(NULL, NULL);

	int rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);


	Block2d<double> b{3,4};
	
	if (rank==0)	
		b.assign({0,1,2,3,4,5,6,7,8,9,10,11});

	std::cout<<"Rank="<<rank<<"\n"<<endl;
	b.print();

	constexpr int dim = 2; // 2D
	int array_size[dim] = {b.sizeX, b.sizeY};
	int subarray_size[dim] = {2,2};
	int subarray_start[dim] = {1,2}; // address of 6
	
	// set where you want to paste the 
	// subarray in rank 1
	if (rank==1){
		subarray_start[0]=0;
		subarray_start[1]=1;
	}

	MPI_Datatype subtype;

	MPI_Type_create_subarray(dim, array_size, subarray_size, subarray_start,
		MPI_ORDER_C, MPI_DOUBLE, &subtype);
	MPI_Type_commit(&subtype);

	if (rank==0)	
	   		MPI_Send( &b(0,0) , 1 , subtype , 1 , 0 , MPI_COMM_WORLD);
	else if(rank==1)
	        MPI_Recv( &b(0,0), 1 , subtype, 0, 0 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	std::cout<<"After communications:\n"<<endl;
	b.print();
	

	MPI_Finalize();
}

