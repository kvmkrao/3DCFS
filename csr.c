#include<stdlib.h>

void CSRSETUP(ien,ia,ja,idiag,eloc,nenf,nef,nnf,nz,flag)
int *ien, *ia, *ja, *idiag,  *eloc;
int *nenf,*nef,*nnf,*nz,*flag;
{
#define ien(i,j) (*(ien+(j-1)*nen+i-1))
#define eloc(i,j,k) (*(eloc+(k-1)*nen*nen+(j-1)*nen+i-1))
#define ia(i)    ia[i-1]
#define ja(i)    ja[i-1]
#define idiag(i) idiag[i-1]
#define hsh(i)   hshptr[i-1]
#define n2e(i)   n2eptr[i-1]
#define n2l(i)   n2lptr[i-1]
int ie,inl,inode,jnl,jnode,iloc,ist,ind,itmp1,itmp2;
int *hshptr,jmin,jmax;
int *n2eptr, *n2lptr;
int n2ls;

int nen,ne,nn;

nen = *nenf; ne = *nef; nn = *nnf; 

/* FORM IA ARRAY AND FIND SIZE OF JA */
    if(*flag==1)
    {

/* FORM N2E ARRAY */
        n2eptr = (int*) malloc((nn+1)*sizeof(int));
//    	iclear(n2eptr,nn+1);
        for(ie=1;ie<=ne;ie++)
        {
	    for(inl=1;inl<=nen;inl++) 
	    { 
	        inode = ien(inl,ie);
                if(inode>0)n2e(inode)+=1;
            }
        }

	itmp1 = n2e(1);
	n2e(1) = 1;
	for(inode=1;inode<=nn;inode++) 
        {
	    itmp2 = n2e(inode+1);
	    n2e(inode+1) = n2e(inode) + itmp1;
	    itmp1 = itmp2; 	
        }


	if((n2e(nn+1)-1) != nen*ne)
        {
          printf("MESH CONTAINS MIXED ELEMENTS! \n");
	  
        }

/* FORM N2L ARRAY */
        n2ls   = n2e(nn+1)-1;
        n2lptr = (int*) malloc(n2ls*sizeof(int));
//        iclear(n2lptr,n2ls);
        for(ie=1;ie<=ne;ie++)
        {
	    for(inl=1;inl<=nen;inl++) 
	    { 
	        inode = ien(inl,ie);
                if(inode>0)
                {
                    iloc = n2e(inode);
                    while(n2l(iloc)!=0) iloc++;
	            n2l(iloc) = ie;
                }
            }
        }

/* FORM IA */
//        iclear(ia,nn+1);
        hshptr   = (int*)malloc(sizeof(int)*nn);
//        iclear(hshptr,nn);
        for(inode=1;inode<=nn;inode++)
        {
            ist = n2e(inode);
            ind = n2e(inode+1) -1;
            for(iloc=ist;iloc<=ind;iloc++)
            {
                ie = n2l(iloc);
                for(inl=1;inl<=nen;inl++)
                {
                    jnode = ien(inl,ie);
                    if(jnode>0)hsh(jnode) = inode;
                }
             } 
 
             for(jnode=1;jnode<=nn;jnode++)
             {    
                 if(hsh(jnode)==inode) ia(inode)+=1;
             }
         }
            
	itmp1 = ia(1);
	ia(1) = 1;
	for(inode=1;inode<=nn;inode++) 
        {
	    itmp2 = ia(inode+1);
	    ia(inode+1) = ia(inode) + itmp1;
	    itmp1 = itmp2; 	
        }
		
        *nz = ia(nn+1) -1;

        free(hshptr);
        free(n2eptr);
        free(n2lptr);
        return;
    }	

/* FORM JA ARRAY */
//    iclear(ja,*nz);
    for(ie=1;ie<=ne;ie++)
    {
        for(inl=1;inl<=nen;inl++)
        {
            inode = ien(inl,ie);
            if(inode>0)
            {
            ist = ia(inode);
            ind = ia(inode+1) -1;		
            for(jnl=1;jnl<=nen;jnl++)
            {
                jnode  = ien(jnl,ie);
                if(jnode>0)
		{
		iloc   = ist;
		while(ja(iloc) !=0 && ja(iloc) !=jnode) iloc++;
		if(iloc>ind)
                {
                printf("ERROR IN JA SETUP, ABORTING! \n");	
		exit(1);
		}
	        if(ja(iloc) ==0) ja(iloc) = jnode;
		}
            }
            }
        }
    }


/* SORT ENTRIES IN JA IN ASCENDING ORDER */
    hshptr   = (int*)malloc(sizeof(int)*nn);
//    iclear(hshptr,nn);
    for(inode=1;inode<=nn;inode++)
    {
        ist = ia(inode);
	ind = ia(inode+1)-1;
        jmin =  nn;
        jmax =   1;

        for(iloc=ist;iloc<=ind;iloc++)
        {
            jnode = ja(iloc);
	    hsh(jnode) = inode;
            if(jnode<jmin) jmin = jnode;
            if(jnode>jmax) jmax = jnode;
        }


        iloc = ist;
        for(jnode=jmin;jnode<=jmax;jnode++)
        {
            if(hsh(jnode)==inode)
            {
                ja(iloc) = jnode;
                iloc ++;
            }
        }
    }
    free(hshptr);

/*FORM IDIAG ARRAY */
    for(inode=1;inode<=nn;inode++)
    {
        iloc = ia(inode);
	ind = ia(inode+1)-1;
        while(ja(iloc)!=inode) iloc++;
        idiag(inode) = iloc;
    }
         
	
/* FORM ELOC ARRAY */
//    iclear(eloc,nen*nen*ne);
    for(ie=1;ie<=ne;ie++)
    {
        for(inl=1;inl<=nen;inl++)
        {
            inode = ien(inl,ie);
            if(inode>0)
            {
                ist = ia(inode);
                ind = ia(inode+1)-1;
                for(jnl=1;jnl<=nen;jnl++)
                {
                    jnode = ien(jnl,ie);
                    if(jnode>0)
                    {
		        iloc  = ist;
			while(ja(iloc)!=jnode) 
                        {
                            iloc++;
                            if(iloc>ind)
                            {
                                printf("ERROR FORMING ELOC, ABORTING !\n");
				exit(1);
                            }
                        }
				
			eloc(jnl,inl,ie) = iloc;
                    }
                 }
              }
         }
                  
      }

return;
}
                
         
