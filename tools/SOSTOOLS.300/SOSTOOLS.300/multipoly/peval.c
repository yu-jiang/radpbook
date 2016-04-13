/*
 * =============================================================
 * peval.c
 * C implementation of peval.m
 *
 * INPUTS 
 *    x: nvars-by-npts matrix each column of which specifies a point
 *       at which to evaluate the polynomial
 *    coef: nterms-by-lp coefficient matrix where lp is the length
 *       of the polynomial
 *    deg: nterms-by-nvars degree matrix
 *
 * OUTPUTS  
 *    val: lp-by-npts matrix where each column specifies the value of
 *       the polynomial evaluated at the corresponding column of x.
 *
 * PJS 8/7/06
 * =============================================================
 */

#include <math.h>
#include "mex.h"

/* Computational subroutine for full (i.e. non-sparse) matrices */
void fullworkFcn(
    double *xp,             /* nvars-by-npts matrix of points, x */
    double *coefp,          /* nterms-by-lp coef matrix */
    double *degp,           /* nterms-by-nvars degree matrix */
    double *val,            /* lp-by-npts matrix of evaluated values */
    int npts,               /* # of points at which to evaluate p */
    int lp,                 /* length of polynomial */
    int nvars,              /* # of variables in the polyomial */
    int nterms             /* # of terms in the polynomial */
    )
{
    
    /* Define variables */
    double monom, p; 
    int i1, i2, i3, i4, exp;
    
    /* Evaluate polynomial at each point */
    for (i1=0; i1<npts; i1++) {        
        /* Evaluate each entry of polynomial vector */
        for (i2=0; i2<lp; i2++){            
            /* Evaluate each term of this entry */
            p = 0;
            for (i3=0; i3<nterms; i3++){
                /* Evaluate each monomial of this term */
                monom  = coefp[i3+nterms*i2];
                if (monom!=0){
                    for (i4=0; i4<nvars; i4++){
                        exp = degp[i3+nterms*i4];
                        if (exp!=0){
                            monom = monom*pow(xp[i4+nvars*i1],exp);
                        }
                    }
                    p = p+monom;                
                }
            }
            val[i2+lp*i1]=p;
            /*mexPrintf("\n i1=%i, i2=%i, val(i2,i1)=%4.3f\n",i1,i2,val[i2+lp*i1]);*/
        }
    }
}


/* Computational subroutine for sparse matrices */
void sparseworkFcn(
    double *xp,             /* nvars-by-npts matrix of points, x */
    double *coefp,          /* nterms-by-lp coef matrix */
    int coefnz,            /* # of nonzero elements of coef */
    mwIndex *coefir,            /* row indices for non-zero elts of coef */
    mwIndex *coefjc,            /* col indices for non-zero elts of coef */
    double *degp,           /* nterms-by-nvars degree matrix */
    int degnz,            /* # of nonzero elements of deg */
    mwIndex *degir,            /* row indices for non-zero elts of deg */
    mwIndex *degjc,            /* col indices for non-zero elts of deg */
    double *val,            /* lp-by-npts matrix of evaluated values */
    int npts,               /* # of points at which to evaluate p */
    int lp,                 /* length of polynomial */
    int nvars,              /* # of variables in the polyomial */
    int nterms             /* # of terms in the polynomial */
    )
{
    
    /* Define variables */
    double *monom;    
    int *coefidx, *degidx;
    int pt, i1, i2, cnt, ncols;
    int term, var, pidx;

        
    /* Create matrix to hold monomial  evaluations */
    monom = (double *) calloc(nterms, sizeof(double));    

    /* Create matrices to hold column indices */    
    coefidx = (int *) calloc(coefnz, sizeof(int));    
    cnt = 0;
    for (i1=0; i1<lp; i1++) {
        ncols = coefjc[i1+1]-coefjc[i1];
        for (i2=0; i2<ncols; i2++){
            coefidx[cnt+i2] = i1; 
        }        
        cnt = cnt+ncols;
    }
    
    degidx = (int *) calloc(degnz, sizeof(int));    
    cnt = 0;
    for (i1=0; i1<nvars; i1++) {
        ncols = degjc[i1+1]-degjc[i1];
        for (i2=0; i2<ncols; i2++){
            degidx[cnt+i2] = i1; 
        }        
        cnt = cnt+ncols;
    }
            
    /* Evaluate polynomial at each point */
    for (pt=0; pt<npts; pt++) {
        /* Compute all monomimials at this point */
        for (i1=0; i1<nterms; i1++){
            monom[i1]=1;
        }
        for (i1=0; i1<degnz; i1++) {
            term = degir[i1];
            var = degidx[i1];
            monom[term] = monom[term]*pow(xp[var+nvars*pt],degp[i1]);
           }
        
        /* Compute each entry of polynomial */
        for (i1=0; i1<coefnz; i1++) {
            term = coefir[i1];
            pidx = coefidx[i1];
            val[pidx+lp*pt] = val[pidx+lp*pt]+coefp[i1]*monom[term];                        
            /*mexPrintf("\n pt = %i, i1=%i, i2=%i, val(i1,i2)=%4.3f\n",pt,pidx,pt,val[pidx+lp*pt]);*/            
        }
    }
    free(coefidx); 
    free(monom);   
    free(degidx);  
}
  

/* The gateway routine: Setup and Error Checking */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    double  *xp, *coefp, *degp; 
    mwIndex coefnz, *coefir, *coefjc;
    mwIndex degnz, *degir, *degjc;
    mwSize szxr, szcoefr, szdegr, szdegc;
    int nvars, nterms;
    mwSize npts, lp;
    double *val;
    
    /* Check for the proper number of arguments. */
    if (nrhs != 3)
        mexErrMsgTxt("Three inputs are required: val = peval(x,coef,deg)");

    /* Check that first input, x, is a real double array */ 
    szxr = mxGetM(prhs[0]);
    npts = mxGetN(prhs[0]);
    if ( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ){
        mexErrMsgTxt("First argument must be a real double array.");    
    }    
    xp = mxGetPr(prhs[0]);    
        
    /* Check that second argument, coef, is a real double array */ 
    szcoefr = mxGetM(prhs[1]);
    lp = mxGetN(prhs[1]);
    if ( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ){
        mexErrMsgTxt("Second argument must be a real double array.");    
    }
    coefp = mxGetPr(prhs[1]);    
    
    /* Check that third argument, deg, is a real double array */ 
    szdegr = mxGetM(prhs[2]);
    szdegc = mxGetN(prhs[2]);
    if ( !mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]) ){
        mexErrMsgTxt("Third argument must be a real double array.");    
    }
    degp = mxGetPr(prhs[2]);    
    
    /* Check that dimensions match */
    if (szxr != szdegc){
        mexErrMsgTxt("Row dimension of x does not match column dimension of deg.");    
    }
    else{
        nvars = szxr;
    }
    if (szcoefr != szdegr){
        mexErrMsgTxt("Row dimension of coef does not match row dimension of deg.");    
    }
    else{
        nterms = szcoefr;
    }
    
    /*  Set the output pointer to the output matrix. */
    plhs[0] = mxCreateDoubleMatrix(lp,npts, mxREAL);

    /*  Create a C pointer to a copy of the output matrix. */
    val = mxGetPr(plhs[0]);
            
    /* Call the C subroutine to evaluate polynomial at sample pts */
    if (mxIsSparse(prhs[1]) && mxIsSparse(prhs[2])){
        //coefnz = mxGetNzmax(prhs[1]);   /* Max # of nonzeros */
        coefir = mxGetIr(prhs[1]);
        coefjc = mxGetJc(prhs[1]);
        coefnz = *(coefjc +lp);       /* Actual # of nonzeros */
        
        //degnz = mxGetNzmax(prhs[2]);   /* Max # of nonzeros */
        degir = mxGetIr(prhs[2]);
        degjc = mxGetJc(prhs[2]);
        degnz = *(degjc + szdegc);       /* Actual # of nonzeros */
        
        sparseworkFcn(xp,coefp,coefnz,coefir,coefjc,degp,degnz,degir,degjc,val,npts,lp,nvars,nterms);
    }
    else if (!mxIsSparse(prhs[1]) && !mxIsSparse(prhs[2])){
        fullworkFcn(xp,coefp,degp,val,npts,lp,nvars,nterms);    
    }
    else {
        mexErrMsgTxt("coef and deg must either both be full or both be sparse matrices");    
    }
    return;
}

