//Based on code in slides
#include <stdio.h>
#include "omp.h"

int fib( int n)
{
    int i, j;
    if (n < 2) {
        return n;
    }
    else {
        #pragma omp task shared(i) firstprivate(n)
        {
        i = fib(n-1);
        }

        #pragma omp task shared(j) firstprivate(n)
        {
            j = fib(n-2);
        }

        #pragma omp taskwait
            return i+j;
    } 
}

int main()
{
    int n, f;
    printf("Enter number:");
    scanf("%d", &n);
    #pragma omp parallel shared(n,f)
    {
        #pragma omp single
        f = fib(n);
    }
    printf("\nfib(%d) = %d\n", n, f);
    return 0;
}
