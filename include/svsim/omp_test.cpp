#include <omp.h>
#include <stdio.h>

int main()
{
    // #pragma omp parallel
    //     {
    int num_threads = omp_get_max_threads();
    printf("Number of threads = %d\n", num_threads);

    return 0;
}