#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>
#include "site.h"
#include "stack.h"
#include <omp.h>

#define MASTER 0 // The Master Node has the lowest index 
#define MAX_SIZE_TAG 0 // The Tag for Cluster Size 
#define PERC_TAG 1 // The Tag for Percolation

// MPI can't transmit bools, so we'll have to use ints.
#define TRUE 1 // True is represented by 1 as in traditional C standards
#define FALSE 0 // False is reprsented by 0 as in traditional C standards

// Free every row of a 2D matrix
void free_matrix(void** matrix, int size) 
{
	int i; // Declare a matrix index
	#pragma omp parallel for
	for (i=0; i<size; i++) // Iterate though every cell
	{
		free(matrix[i]); // Free the memory allocated to each cell
	}
	free(matrix); // Free the matrix allocation
	matrix = NULL; // Nullify the pointer
}

// Custom modulus function
// This handles negative numbers which C's implementation does not.
int modulus(int a, int b)
{
	if(a%b >= 0) return (a%b); // If the remainder operation result is positive return it verbatim
	else return(b + a); // Else retun the correct modulus result
}

// Return TRUE if a randomly generated number is less than p.
int rng(double p)
{
	int random; // Declare an integer
	random = rand(); // Assign the integer to a random number
	return (random <= (p * RAND_MAX)); // Return the result of whether the algorithm determined true or false
}

// Returns the largest of 2 doubles.
double max(double a, double b) 
{
	return a > b? a: b; // Return a if is larger else return b
}

// Turns a row and col index into 1 index for our 1D lattice
int getIndex(int row, int col, int size)
{
	return(row*size + col); // Return the correct index for a given matrix location
}

// Get a site pointer from the 1d lattice by row and column
SITE* getSite(SITE* matrix, int row, int col,int size)
{
	return &matrix[getIndex(row, col, size)]; // Return the site at a given matrix location
}

// Populate the matrix with SITE structs, calculate bonds for either site or bond percolation.
// This is called by the master node.
void populate_matrix(SITE* matrix, int size, double probability, int site) 
{
	srand(time(NULL)); 
	int numBlocks; 
	int blockSize;

	//Handle each block parallel to each other but sequentially within the block itself
	int block;
	numBlocks = 1;

	// Work out the optimal number of threads to use based on the lattice size
	int divisor;
	for(divisor = 12; divisor > 1; divisor--)
	{
		if(size % divisor == 0)
		{
			numBlocks = divisor;
			break;
		}
	}

	blockSize = size / numBlocks;
	if(blockSize < 1) blockSize = 1;


	#pragma omp parallel for
	for(block = 0; block < numBlocks*numBlocks; block++)
	{
		int row;
		for (row = (block/numBlocks)*blockSize; row < (block/numBlocks + 1)*(blockSize); row++) 
		{
			int col;
			for (col = (block%numBlocks)*blockSize; col < (block%numBlocks + 1)*(blockSize); col++) 
			{
				SITE *currentSite = getSite(matrix, row, col, size);
				currentSite->row = row;
				currentSite->col = col;
				currentSite->populated = rng(probability);
				if(currentSite->populated) 
				{
					if ((col%blockSize > 0) && ((site) || rng(probability))) 
					{
						SITE* west = getSite(matrix, row, modulus(col-1, size), size); 
						currentSite->w = west->populated;
						west->e = currentSite->w;
					}
					if ((row%blockSize > 0) && ((site) || rng(probability))) 
					{
						SITE* north = getSite(matrix, modulus(row-1, size), col, size); 
						currentSite->n = north->populated;
						north->s = currentSite->n;
					}
				}
			}
		}
	}
	
	//printf("Completed Parallel Block Population.\n");

	// Handle the Edge Cases where Boundary Percolation is relevant

	int row;
	int col;
	for(row = 0; row < numBlocks*blockSize; row += blockSize)
	{
		for(col = 0; col < numBlocks*blockSize; col++)
		{
			SITE *currentSite = getSite(matrix, row, col, size);
			currentSite->row = row;
			currentSite->col = col;
			if(currentSite->populated) 
			{
				SITE* current = getSite(matrix, row, col, size);
				if ((row%blockSize == blockSize-1) && ((site) || rng(probability)))
				{
					SITE* south = getSite(matrix, modulus(row+1, size), col, size);
					current->s = south->populated;
					south->n = current->s;
				}
				if ((col%blockSize == blockSize-1) && ((site) || rng(probability)))
				{
					SITE* east = getSite(matrix, row, modulus(col+1, size), size);
					current->e = east->populated;
					east->w = current->e;
				}
			}
		}
	}
	//printf("Completed Row, Column Edge Population.\n");

	for(row = 0; row < numBlocks*blockSize; row++)
	{
		for(col = 0; col < numBlocks*blockSize; col+= blockSize)
		{
			SITE *currentSite = getSite(matrix, row, col, size);
			currentSite->row = row;
			currentSite->col = col;
			if(currentSite->populated) 
			{
				SITE* current = getSite(matrix, row, col, size);
				if ((row%blockSize == blockSize-1) && ((site) || rng(probability)))
				{
					SITE* south = getSite(matrix, modulus(row+1, size), col, size);
					current->s = south->populated;
					south->n = current->s;
				}
				if ((col%blockSize == blockSize-1) && ((site) || rng(probability)))
				{
					SITE* east = getSite(matrix, row, modulus(col+1, size), size);
					current->e = east->populated;
					east->w = current->e;
				}
			}
		}
	}
	//printf("Completed Column, Row Edge Population\n");

}

// Flood the matrix from a starting row and column, mark all seen sites in the visited matrix.
// This is done with a depth first search.
// Return the size of the cluster it formed.
RESULTS flood(SITE* lattice, int** visited, int size, ROUTE routing, int row, int col) {
	int cluster_size = 0;
	int* colsVisited = (int*) calloc(size, sizeof(int));
	int* rowsVisited = (int*) calloc(size, sizeof(int));
	int nRowsVisited = 0;
	int nColsVisited = 0;
	int allRowsVisited;
	int allColsVisited;

	SITE_STACK* stack = (SITE_STACK*) malloc(sizeof(SITE_STACK));
	init_stack(stack);
	visited[row][col] = TRUE;
	SITE* current = getSite(lattice, row, col, size);
	push(stack, current);

	// Perform the depth first search
	while(!is_empty(stack)) {
		SITE* site = pop(stack);
		cluster_size++;
		row = site->row;
		col = site->col;
		if (!colsVisited[col]) {
			colsVisited[col] = TRUE;
			nColsVisited++;
		}
		if (!rowsVisited[row]) {
			rowsVisited[row] = TRUE;
			nRowsVisited++;
		}

		// Push all neighbours onto the stack iff they have not been seen yet.
		if (site->e && !visited[row][modulus(col+1, size)]){
			visited[row][modulus(col+1, size)] = TRUE;
			SITE* east = getSite(lattice, row, modulus(col+1, size), size);
			push(stack, east);
		}
		if (site->s && !visited[modulus(row+1, size)][col]){
			visited[modulus(row+1, size)][col] = TRUE;
			SITE* south = getSite(lattice, modulus(row+1, size), col, size);
			push(stack, south);
		}
		if (site->w && !visited[row][modulus(col-1, size)]){
			visited[row][modulus(col-1, size)] = TRUE;
			SITE* west = getSite(lattice, row, modulus(col-1, size), size);
			push(stack, west);
		}
		if (site->n && !visited[modulus(row-1, size)][col]) {
			visited[modulus(row-1, size)][col] = TRUE;
			SITE* north = getSite(lattice, modulus(row-1, size), col, size);
			push(stack, north);
		}
	}

	// Free all dynamically allocated memory
	free_stack(stack);
	free(colsVisited);
	free(rowsVisited);	 

	allColsVisited = nColsVisited == size;
	allRowsVisited = nRowsVisited == size;

	int percolated = FALSE;
	// Calculate percolation for all percolation types.
	switch(routing) {
		case ROWS:
			percolated = allRowsVisited;
			break;
		case COLS:
			percolated = allColsVisited;
			break;
		case BOTH:
			percolated = (allRowsVisited && allColsVisited);
			break;
	}
	RESULTS result;
	result.percolates = percolated;
	result.size = cluster_size;
	return result;
}

// Performs a flood for every site in a given range of the lattice.
// Sends messages to the master node, called by each slave node.
int get_max_cluster(SITE* lattice, int size, int lowerBound, int upperBound, ROUTE routing) {
	int max = 0;
	int complete = FALSE;

	// The same visited matrix is used by each flood.
	// This is so we can avoid exploring starting points that are a part of a previously visited cluster.
	// Clusters can never overlap so this is fine.
	int** visited = (int**) calloc(size, sizeof(int*));
	for (int i = 0; i<size; i++) {
		visited[i] = (int *) calloc(size, sizeof(int));
	}

	int square_size = size * size;

	for(int count=lowerBound; count<upperBound; count++) {
		int row, col;
		row = count / size;
		col = count % size;
		/*
		* Consider each unvisited site within this process's allocated chunk as a starting point.
		* Perform a depth first search from each starting point, recording the max cluster size.
		* Lattice percolates if any clusters percolate.
		*/
		SITE* thisSite = getSite(lattice, row, col, size);
		if (!visited[row][col] && thisSite->populated) {
			RESULTS result = flood(lattice, visited, size, routing, row, col);
			int cluster_size = result.size;
			if (cluster_size > max) max = cluster_size;
			complete = complete || result.percolates;
		}
	}
	free_matrix((void**) visited, size);

	// Send max size
	MPI_Send(&max, 1, MPI_INT, MASTER, MAX_SIZE_TAG, MPI_COMM_WORLD);

	// Send Percolated
	MPI_Send(&complete, 1, MPI_INT, MASTER, PERC_TAG, MPI_COMM_WORLD);
	return max;
}

// Create an MPI type for our SITE struct.
// Declare the lattice as an array of this type when broadcasting it.
MPI_Datatype create_mpi_site_struct() {
	int nItems = 7;
	int blocklengths[7] = {1,1,1,1,1,1,1};
	MPI_Datatype types[7] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
	MPI_Datatype MPI_SITE;
	MPI_Aint offsets[7];
	for (int i=0; i<7; i++) {
		offsets[i] = sizeof(int) * i;
	}
	MPI_Type_create_struct(nItems, blocklengths, offsets, types, &MPI_SITE);
	MPI_Type_commit(&MPI_SITE);
	return MPI_SITE;
}	

// Print out the lattice (Only populated, not bonds)
void print_lattice(SITE* lattice, int size) {
	for (int i=0; i<size; i++) {
		for (int j=0; j<size; j++) {
			printf("%i ", getSite(lattice, i, j, size)->populated);
		}
		printf("\n");
	}
}

int main(int argc, char* argv[]) {
	//Start timer
	struct timeval start, end;
	gettimeofday(&start, NULL);

	MPI_Init(&argc, &argv);
	int numProcs, numId;
	MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
	MPI_Comm_rank(MPI_COMM_WORLD, &numId);
	MPI_Datatype MPI_SITE = create_mpi_site_struct();
	char type;
	char path;
	int size;
	double p;
	SITE* matrix;

	if (argc >= 4) {
		type = argv[1][0];
		path = argv[2][0];
		size = atoi(argv[3]);
		p = strtod(argv[4], NULL);

		matrix = (SITE*) calloc(size*size, sizeof(SITE));

		// Parse percolation type (rows, columns, all)
		ROUTE routing;
		switch(path) {
			case 'r':
				routing = ROWS;
				break;
			case 'c':
				routing = COLS;
				break;
			case 'd':
				routing = BOTH;
				break;
			default:
				printf("Invalid path argument\n");
				return -1;
		}

		if (numId == MASTER) {
			// The master process will seed the lattice and broadcast it to all other processes.
			populate_matrix(matrix,size, p, type == 's');
			//print_lattice(matrix, size);
		}

		// All other processes receive the lattice as a 1d array.
		MPI_Bcast(matrix, size*size, MPI_SITE, MASTER, MPI_COMM_WORLD);

		if (numId == MASTER) {
			// Master node's branch
			int count = 0;
			int maxSize = 0;	
			int finalPerc = FALSE;

			// Don't wait for processes that have exited early.
			int activeProcs = (size*size - 1) > (numProcs - 1)? (numProcs -1): (size*size - 1);

			printf("Active processes: %i\n", activeProcs);

			// Wait for all active slave processes to report their max size.
			while (count < activeProcs) {
				int size;
				MPI_Status status;
				MPI_Recv(&size, 1, MPI_INT, MPI_ANY_SOURCE, MAX_SIZE_TAG, MPI_COMM_WORLD, &status);
				if (numId == MASTER)
					printf("Process %i sent max_size %i\n", status.MPI_SOURCE, size);
				if (size > maxSize)
					maxSize = size;
				count++;
			}

			count = 0;
			// Wait for all active slave processes to report if they found a cluster that percolates.
			while (count < activeProcs) {
				int perc;
				MPI_Status status;
				MPI_Recv(&perc, 1, MPI_INT, MPI_ANY_SOURCE, PERC_TAG, MPI_COMM_WORLD, &status);
				if (numId == MASTER)
					printf("Process %i received percolates %i\n", numId, perc);

				finalPerc = (finalPerc || perc);
				count++;
			}
			printf("Max cluster: %i\n", maxSize);
			printf("Percolates: %c\n", finalPerc? 'Y': 'N');
		} else {
			// Slave node branch
			double chunk = max(1.0, (double)size*size / (numProcs-1));
			int lower = (int) (chunk * (numId - 1));
			int upper = (int) (lower + chunk);

			// If the calculated chunk is bigger than the total number of sites, do nothing.
			if (upper < size*size && lower < size*size) {
				printf("Process %i exploring site %i - %i\n", numId, lower, upper);
				get_max_cluster(matrix, size, lower, upper, routing);
			}
		}
	} else {
		printf("Invalid arguments\n");
	}

	// end timer
	gettimeofday(&end, NULL);
	double delta = ((end.tv_sec  - start.tv_sec) * 1000000u +
		   	end.tv_usec - start.tv_usec) / 1.e6;

	if (numId == MASTER)
		printf("time=%12.10f\n",delta);	
	MPI_Finalize();

}