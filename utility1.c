#include<stdio.h>
#include<malloc.h>
#include<stdlib.h>
#include<string.h>
#include<unistd.h>
#if defined(UNI)
#else
#include<mpi.h>
#endif
#define MAXLEN 255

static int my_pe,npes;

void setup_parallel_io(void)
{
#if defined(UNI)
	my_pe=0; npes=1;
#else
        MPI_Comm_rank(MPI_COMM_WORLD,&my_pe);
        MPI_Comm_size(MPI_COMM_WORLD,&npes);
#endif
}


FILE *infile;
void open_input_file(fname)
char *fname;
{
    infile = fopen(fname,"r");
    return;
}

void shut_input_file(void)
{
    fclose(infile);
    return;
}


/* FUNCTION TO READ A FILE OF DOUBLES IN PARALLEL */
void frd(fname, data, dim1,dim2)
char* fname; double *data; int dim1,dim2;
{
FILE *fp;
int i,nf,nr,dsize,offset; 

#if defined(UNI)
        my_pe=0; npes=1;
#else
        MPI_Comm_rank(MPI_COMM_WORLD,&my_pe);
        MPI_Comm_size(MPI_COMM_WORLD,&npes);
#endif 

        nf = (dim2 - 1)/npes + 1;
        nr = dim2 - (npes-1)*nf;

        nf = dim1*nf;
        nr = dim1*nr;

/*READ IEEE DATA */

        /* READ DATA ON PN 0 AND SEND TO OTHERS */
        for(i=1;i<npes;i++){

                if(i!=npes-1) dsize = nf; else dsize = nr;

                if(my_pe==i) {
                        offset = i*nf*8;
                        fp = fopen(fname,"r");
                        fseek(fp,offset,SEEK_SET);
                        fread(data,sizeof(double),dsize,fp);
#if defined(I386)
//      			swaparr((char *) data, dsize, sizeof(double));
#endif
                        fclose(fp);
                }
#if defined(UNI)
#else
                MPI_Barrier(MPI_COMM_WORLD);
#endif
        }


/* READ DATA FOR PN 0 */
        if(my_pe==0){
                fp = fopen(fname,"r");
                fread(data,sizeof(double),nf,fp);
#if defined(I386)
//     		swaparr((char *) data, nf, sizeof(double));
#endif
                fclose(fp);
                printf("READ FILE: %s \n",fname);
        }


return;
}


/* FUNCTION TO WRITE A FILE OF DOUBLES IN PARALLEL */
void fwd_(char *fname, double *data, int dim1,int dim2)
{
FILE *fp;
#if defined(UNI)
#else
MPI_Status stat;
#endif
int i,nf,nr,dsize;

#if defined(UNI)
        my_pe=0; npes=1;
#else
        MPI_Comm_rank(MPI_COMM_WORLD,&my_pe);
        MPI_Comm_size(MPI_COMM_WORLD,&npes);
#endif
   nf = (dim2 - 1)/npes + 1;
   nr = dim2 - (npes-1)*nf;
   nf = dim1*nf;
   nr = dim1*nr;

   /* WRITE DATA FOR PN 0 */
  if(my_pe==0)
  {
     fp = fopen(fname,"w");
#if defined(I386)
//     swaparr((char *) data, nf, sizeof(double));
#endif
     fwrite(data,sizeof(double),nf,fp);
#if defined(I386)
//     swaparr((char *) data, nf, sizeof(double));
#endif
  }

  /* TRANSFER DATA TO PN 0  WHICH WRITES TO FILE */
  for(i=1;i<npes;i++)
  {

     if(i!=npes-1) dsize = nf; else dsize = nr;

     if (my_pe==i) {
#if defined(UNI)
#else
         MPI_Send(data,dsize,MPI_DOUBLE,0,100,MPI_COMM_WORLD);
#endif
		   }
     else if(my_pe==0) 
     {
#if defined(UNI)
#else
        MPI_Recv(data,dsize,MPI_DOUBLE,i,100,MPI_COMM_WORLD,&stat);
#endif
#if defined(I386)
//        swaparr((char *) data, dsize, sizeof(double));
#endif
        fwrite(data,sizeof(double),dsize,fp);
#if defined(I386)
//        swaparr((char *) data, dsize, sizeof(double));
#endif
     }
#if defined(UNI)
#else
     MPI_Barrier(MPI_COMM_WORLD);
#endif
   }

   if(my_pe==0){ fclose(fp); printf("WROTE FILE: %s \n",fname);}

return;
}



/* FUNCTION TO READ A FILE OF DOUBLES  (SINGLE PN)*/
void sfrd_(fname,data,dim1,dim2)
char *fname; double *data; int dim1; int  dim2;
{
FILE *fp;
int nf;


    nf = dim1*dim2;
	
    fp = fopen(fname,"r");
    fread(data,sizeof(double),nf,fp);
#if defined(I386)
//    swaparr((char *) data, nf, sizeof(double));
#endif
    fclose(fp);
    printf("READ FILE(Serial): %s \n",fname);

return;
}

/* FUNCTION TO WRITE A FILE OF DOUBLES (SINGLE PN) */
void sfwd(fname,data,dim1,dim2)
char *fname; double *data; int dim1; int  dim2;
{
FILE *fp;
int nf;


    nf = dim1*dim2;
	
    fp = fopen(fname,"w");
#if defined(I386)
//    swaparr((char *) data, nf, sizeof(double));
#endif
    fwrite(data,sizeof(double),nf,fp);
#if defined(I386)
//    swaparr((char *) data, nf, sizeof(double));
#endif
    fclose(fp);
    printf("WROTE FILE(Serial): %s \n",fname);

return;
}



/* FUNCTION TO GET A KEY AND INPUT TERMINATION SIGNAL */
void rdchar(key)
char *key;
{
int i,j,ch;


   if(my_pe==0)
   {
      ch =  getc(infile);

      /*SKIP BLANKS AND NEWLINE CHARACTERS */
      while(ch==' ' || ch=='\n') ch= getc(infile);

      /* DISCARD OUTPUT AFTER '!' */
      while(ch=='!' )
      { 
         /* GOTO END OF LINE */
	    while((ch=getc(infile))!='\n');
	 /*SKIP BLANKS AND NEWLINE CHARACTERS */
   	    while(ch==' ' || ch=='\n') ch= getc(infile);
      }


     /* DISCARD OUTPUT AFTER '#' */
      if(ch=='#') 
      {
         *key=ch;
         /* GOTO END OF LINE */
  	 while((ch=getc(infile))!='\n');
	 for(i=1;i<MAXLEN;i++) *(key+i)='\0';
	 goto  comment; 
      }

      /* READ IN STRING */
      j = 0;
      while(ch!=' ' && ch!='\n') 
      {
         *(key+j)=ch;
	 j++;
	 ch= getc(infile);
      }


      /* BLANK OUT THE REMAINDER OF THE STRING */
      for(i=j;i<MAXLEN;i++)*(key+i)='\0';
   }

comment: 
#if defined(UNI)
#else
   MPI_Bcast(key,MAXLEN,MPI_CHAR,0,MPI_COMM_WORLD);
#endif
return;
}		


/* FUNCTION TO READ INTEGER  */
void rdint(key)
int *key;
{
int j,ch;
char temp[MAXLEN];

   if(my_pe==0)
   {

	/* SKIP BLANKS*/
	while ((ch=getc(infile))==' ');

	j = 0;
	while(ch!=' ' && ch!='\n') {
		*(temp+j)=ch;
		j ++;
		ch= getc(infile);
	}

	*(temp+j)='\0';
	*key = atoi(temp);
   }
#if defined(UNI)
#else
   MPI_Bcast(key,1,MPI_INT,0,MPI_COMM_WORLD);
#endif
return;
}		

/* FUNCTION TO READ DOUBLES  */
void rddbl(key)
double *key;
{
int j,ch;
char temp[MAXLEN];


   if(my_pe==0)
   {

	/* SKIP BLANKS*/
	while ((ch=getc(infile))==' ');

	j = 0;
	while(ch!=' ' && ch!='\n') {
		*(temp+j)=ch;
		j++;
		ch= getc(infile);
	}
	
	*(temp+j)='\0';
	*key = atof(temp);
   }
#if defined(UNI)
#else
   MPI_Bcast(key,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
#endif

return;
}		

#if defined(CRAYT3E)
#elif defined  IEEE
#	include "ieeeint.h"
#if defined(SGI)
#endif
#endif
/*##if defined(I386)
#       include "i386.h"
#endif*/
