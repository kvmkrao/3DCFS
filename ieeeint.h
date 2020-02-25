/* FUNCTION TO READ A FILE OF INTEGERS IN PARALLEL */
void fri(fname, data, dim1,dim2)
char* fname; int *data; int dim1,dim2;
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
                        offset = i*nf*4;
                        fp = fopen(fname,"r");
                        fseek(fp,offset,SEEK_SET);
                        fread(data,sizeof(int),dsize,fp);
//			if (my_pe==1)printf("IEEE:before the I386 loop of fri\n");
#if defined(I386)
//                        swaparr((char *) data, dsize, sizeof(int));
//			if (my_pe==1)printf("IEEE:in the I386 loop of fri\n");
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
                fread(data,sizeof(int),nf,fp);
#if defined(I386)
//                swaparr((char *) data, nf, sizeof(int));
#endif
                fclose(fp);
                printf("READ FILE(PAR) : %s\n",fname);
        }

       
return;
}

/* FUNCTION TO READ A FILE OF INTEGERS (SINGLE PN) */
void sfri(fname,data,dim1,dim2)
char *fname; int *data; int dim1; int  dim2;
{
FILE *fp;
int nf;


    nf = dim1*dim2;
       
    fp = fopen(fname,"r");
    fread(data,sizeof(int),nf,fp);
#if defined(I386)
//    swaparr((char *) data, nf, sizeof(int));
#endif
    fclose(fp);
    printf("READ FILE(SER) : %s \n",fname);

return;
}

