#include <malloc.h>
#include <stdio.h>


/* malloc(3C) wrapper */
	int *MYMALLOC(size)
int *size;
{
        int *newptr;
        newptr = (int *) malloc((unsigned) *size);
        if (newptr == 0) {
                printf("malloc: null pointer\n");
                printf("malloc: called with size = %d\n", *size);
                exit(1);
        }
        return(newptr);
}


	void MYFREE(ptr)
int *ptr;
{
        free((void *) *ptr);
        *ptr = 0;
        return;
}
/* Repeated above VK.OCT09.01 */
        int *mymalloc_(size)
int *size;
{
        int *newptr;
        newptr = (int *) malloc((unsigned) *size);
        if (newptr == 0) {
                printf("malloc: null pointer\n");
                printf("malloc: called with size = %d\n", *size);
                exit(1);
        }
        return(newptr);
}


        void myfree_(ptr)
int *ptr;
{
        free((void *) *ptr);
        *ptr = 0;
        return;
}

/***********************************************************************/
#if defined(CRAYT90) || defined(CRAYT3E)
void MCHECK ()
#else
void mcheck_()
#endif
{
#if defined(CRAYT90) || defined(CRAYT3E)
   if (malloc_check(0) != 0 ) {
      printf("ewdmalloccheck: heap is corrupt\n"); tracebk(); abort();
   }
#endif
}
