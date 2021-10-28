#include<malloc.h>
#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<fcntl.h>
#include"mpc.h"
#define MAXLEN 255
#if defined(CRAYT90) || defined(CRAYT3E)
#define FCDTOCP(x) _fcdtocp(x)
#include <fortran.h>
#else
#define FCDTOCP(x) (x)
#endif

char ch;


#if defined(CRAYT90) || defined(CRAYT3E)
void open_input_file_(keycd)
_fcd keycd;
#else
void open_input_file_(char *keycd, int len)
#endif
{
int i,j,ch,my_pe,ilen,npes;
        char *key;
        key = FCDTOCP(keycd);
#if defined(CRAYT90) || defined(CRAYT3E)
        ilen = _fcdlen(keycd);
#else
        ilen=len;
#endif    
	open_input_file(keycd);
	return;
}

/* FUNCTION TO GET A KEY AND INPUT TERMINATION SIGNAL */
#if defined(CRAYT90) || defined(CRAYT3E)
void RDCHAR(keycd)
_fcd keycd;
#else
void rdchar_(char *keycd, int len)
#endif
{
int i,j,ch,my_pe,ilen,npes;
	char *key;
	key = FCDTOCP(keycd);
#if defined(CRAYT90) || defined(CRAYT3E)
	ilen = _fcdlen(keycd);
#else
	ilen=len;
#endif

#if defined(UNI)
	my_pe=0; npes=1;
#else
MPI_Comm_rank(MPI_COMM_WORLD,&my_pe);
#endif

if(my_pe==0) {
	ch =  getchar();
	/*SKIP BLANKS AND NEWLINE CHARACTERS */
	while(ch==' ' || ch=='\n') ch= getchar();
	/* DISCARD OUTPUT AFTER '!' */
	while(ch=='!' ){ 
		while((ch=getchar())!='\n');
		while(ch==' ' || ch=='\n') ch= getchar();
	}


	/* DISCARD OUTPUT AFTER '#' */
	if(ch=='#') {
		*key=ch;
		while((ch=getchar())!='\n');
		for(i=1;i<ilen;i++) *(key+i)=' ';
		return;
	}

	/* READ IN STRING */
	j = 0;
	while(ch!=' ' && ch!='\n') {
		*(key+j)=ch;
		j++;
		ch= getchar();
	}


	/* BLANK OUT THE REMAINDER OF THE STRING */
	for(i=j;i<ilen;i++)*(key+i)=' ';

/*	printf("READ STRING %s \n",key);   */
}
#if defined(UNI)
#else

	 MPI_Bcast(key,ilen,MPI_CHAR,0,MPI_COMM_WORLD); 
#endif
return;
}		


/* FUNCTION TO READ INTEGER  */
#if defined(CRAYT90) || defined(CRAYT3E)
void RDINT(key)
#else
void rdint_(key)
#endif
int *key;
{
int j,ch,myid;
char temp[MAXLEN];


	/* SKIP BLANKS*/
	while ((ch=getchar())==' ');

	j = 0;
	while(ch!=' ' && ch!='\n') {
		*(temp+j)=ch;
		j ++;
		ch= getchar();
	}

	*(temp+j)='\0';
	*key = atoi(temp);
return;
}		

/* FUNCTION TO READ DOUBLES  */
#if defined(CRAYT90) || defined(CRAYT3E)
void RDREAL(key)
#else
void rdreal_(key)
#endif
double *key;
{
int j,ch,myid;
char temp[MAXLEN];



	/* SKIP BLANKS*/
	while ((ch=getchar())==' ');

	j = 0;
	while(ch!=' ' && ch!='\n') {
		*(temp+j)=ch;
		j++;
		ch= getchar();
	}
	
	*(temp+j)='\0';
	*key = atof(temp);
return;
}		

