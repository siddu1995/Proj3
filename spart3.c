#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define MIN(a,b) ((a)<(b)?(a):(b))

#define BLOCK_LOW(id,p,n) ((id)*(n)/(p))

#define BLOCK_HIGH(id,p,n) ( BLOCK_LOW((id)+1,p,n)-1 )

#define BLOCK_SIZE(id,p,n) (BLOCK_LOW( (id)+1, p, n) - BLOCK_LOW( (id), p, n  ) )

#define BLOCK_OWNER(index,p,n) ((((p)*(index)+1)-1 ) / (n) )

int main (int argc, char *argv[])
{
   int index,count,id,p,node;
   unsigned long long int global_count,i,n,high_value,low_value,size,prime,proc0_size,first;
   unsigned long long int low_value_proc0,high_value_proc0,size_proc0,first_proc0;
   double elapsed_time;
   char *marked,*localmarked;
   MPI_Init(&argc, &argv);
   MPI_Barrier(MPI_COMM_WORLD);
   elapsed_time = -MPI_Wtime();
   MPI_Comm_rank(MPI_COMM_WORLD, &id);
   MPI_Comm_size(MPI_COMM_WORLD, &p);
   if(argc != 2) {
      if (!id) printf("Command line: %s <m>\n", argv[0]);
      MPI_Finalize(); exit(1);
   }
   n = atoll(argv[1]);
  // node = atoi(argv[2]);
   low_value_proc0 = 3;
   high_value_proc0 = 3 + BLOCK_HIGH(0,p,n-2)-(BLOCK_HIGH(0,p,n-2)%2);
   size_proc0  = ((high_value_proc0-low_value_proc0)/2)+1;
   proc0_size = (n-2)/(p*2);
   
   low_value = 3 + BLOCK_LOW(id,p,n-2) + BLOCK_LOW(id,p,n-2) % 2;
   high_value = 3 + BLOCK_HIGH(id,p,n-2) - BLOCK_HIGH(id,p,n-2) % 2;
   size = ( high_value - low_value)/2+1;
   localmarked = (char *) malloc(proc0_size);
   
   if (localmarked == NULL) {
      printf ("Cannot allocate enough memory\n");
      MPI_Finalize();
      exit(1);
   }
   
   for (i = 0; i < size; i++) {localmarked[i] = 0;}

   if ((3 + proc0_size) < (int) sqrt((double) n)) {
      if (!id) printf ("Too many processes\n");
      MPI_Finalize();
      exit(1);
   }
   marked = (char *)malloc(size);
   
   if (marked == NULL) {
      printf ("Cannot allocate enough memory\n");
      MPI_Finalize();
      exit(1);
   }
   for (i = 0; i < size; i++) {marked[i] = 0;}
   index = 0;
   prime = 3;
   do {
      if (prime * prime > low_value)
         first = (prime * prime - low_value)/2;
	else{
		if(!(low_value%prime)){
		first = 0;
		}	
      else {
         if ((low_value % prime)%2==0)
                first = prime - (low_value % prime)/2;
         else
                first = (prime - (low_value % prime))/2;
      }
}

      for (i = first; i < size; i += prime) marked[i] = 1;
      if (!id) {
         while (marked[++index]);
         prime = 2*index+3;
      }
if(id){
first_proc0 = (prime*prime - low_value_proc0)/2;
for(i=first_proc0;i<size_proc0;i+=prime)
localmarked[i]=1;
while(localmarked[++index]);
prime=2*index+3;
}
}while (prime * prime <= n);
   count = 0;
   for (i = 0; i < size; i++)
      if (!marked[i]) count++;
   MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,0, MPI_COMM_WORLD);
   elapsed_time += MPI_Wtime();
   if (!id) {
	global_count++;
      printf ("%llu primes are less than or equal to %llu :: Total elapsed time = %10.6f",global_count, n,elapsed_time);
   }
   MPI_Finalize();
   return 0;
}
