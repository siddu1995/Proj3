#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))
#define BLOCK_LOW(id, p, n) ((id)*(n)/(p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id)+1, p,n)-1)
#define BLOCK_SIZE(id, p, n) (BLOCK_LOW((id)+1,p,n)-BLOCK_LOW((id),p,n))
#define BLOCK_OWNER(index,p,n) ((((p) * index)+1)-1)/(n)

int main(int argc, char *argv[])
{
    int index,count,id,p;
    unsigned long long int global_count,i,k,n,high_value,low_value,size,prime,proc0_size,first;
    unsigned long long int low_proc0,high_proc0,size_proc0,first_proc0;
    double elapsed_time;
    char *marked,*localmarked;
    unsigned long long int size_cache, begin_cache, end_cache, dif_cache, low_cache, high_cache;
    size_cache = 1000000;
    begin_cache=0;

    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    if (argc != 2) {
          if (!id) printf ("Command line: %s <m>\n", argv[0]);
          MPI_Finalize(); exit(1);
    }
    n = atoll(argv[1]);

    low_proc0 = 3;
    high_proc0 = (unsigned long long int)sqrt(n);
    size_proc0 = (high_proc0 - low_proc0) / 2 + 1;
    localMarked = (char*)malloc(size_proc0);

    low_value = 3 + BLOCK_LOW(id,p,n-2) + BLOCK_LOW(id,p,n-2) % 2;
    high_value = 3 + BLOCK_HIGH(id,p,n-2) - BLOCK_HIGH(id,p,n-2) % 2;
    size = (high_value - low_value) / 2 + 1;
    proc0_size = ((n-2)/(2*p));
    marked = (char *)malloc(size);

    if (localMarked == NULL) {
        printf("Cannot allocate memory to local array for seiving primes\n");
        MPI_Finalize();
        exit(1);
    }

    for(k=0;k<size_proc0;k++) localMarked[k] = 0;

    if ((3 + proc0_size) < (int) sqrt((double) n)) {
        if (!id) printf ("Too many processes\n");
        MPI_Finalize();
        exit(1);
    }

    if (marked == NULL) {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }
    for (i = 0; i < size; i++){
        marked[i] = 0;
    }
    index = 0;
    prime = 3;
    do{
        first_proc0 = (prime*prime - low_proc0) / 2;
        for(i=first_proc0;i<size_proc0;i+=prime){
            localMarked[i] = 1;
        }
        while(localMarked[++index]);
        prime = 2 * index + 3;
    }while(prime*prime <= n);

    low_cache = low_value;
    high_cache = high_value;
    begin_cache = 0;
    count = 0;
    if(size%size_cache == 0){
        end_cache = size/size_cache;
    }
    else{
        end_cache = (size/size_cache)+1;
    }
    do{
        low_value = ((low_cache)+(begin_cache * size_cache * 2));
        high_value = MIN(high_cache, (low_value + (2*size_cache -2)));
        dif_cache = (high_value - low_value) / 2 + 1;

        for(k=0;k<dif_cache;k++){
            marked[k] = 0;
        }
        index = 0;
        prime = 3;
        do{
            if(localMarked[index] == 0){
                prime = 2 * index + 3;

                if (prime * prime > low_value){
                    first = (prime * prime - low_value)/2;
                }
                else{
                    if (!(low_value % prime)){
                        first = 0;
                    }
                    else{
                        if ((low_value%prime) % 2 == 0)
        				{
        					first = prime - (low_value%prime) / 2;
        				}
        				else
        				{
        					first = (prime - (low_value%prime)) / 2;
        				}
                    }
                }
                for (i = first; i < dif_cache; i += prime){
                    marked[i] = 1;
                }
            }
            index++;
        }while(index < (((high_proc0 - 3) / 2) + 1));
        for (i = 0; i < dif_cache; i++){
            if (!marked[i]){
                count++;
            }
        }
        begin_cache++;
    }while (begin_cache < end_cache);

    if(p>1){
        MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,0, MPI_COMM_WORLD);
    }
    elapsed_time += MPI_Wtime();
    if (!id) {
        global_count++;
        printf("Total number of primes: %llu, Total time: %10.6f sec\n",global_count,elapsed_time);    }
    MPI_Finalize();
    return 0;
}
