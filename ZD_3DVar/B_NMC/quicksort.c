/*  sort array x in ascending order    */
/*  permute integer arrays mt and nt   */
/*  accordingly                        */

int quicksort_(float *x,int *lo,int *hi,int *mt,int *nt)
{
  int ilo=*lo-1,ihi=*hi-1;
  quicksrt(x,ilo,ihi,mt,nt);
  return 0;
}

int quicksrt(float *x,int lo,int hi,int *mt,int *nt)
{
  float pivot,temp;
  int i=lo, j=hi, itmp;

  pivot=x[(lo+hi)/2];
  do{
    while (i<hi && x[i]<pivot) ++i;
    while (j>lo && x[j]>pivot) --j;
    if(i<=j){
      temp=x[i]; x[i]=x[j]; x[j]=temp;
      itmp=mt[i]; mt[i]=mt[j]; mt[j]=itmp;
      itmp=nt[i]; nt[i]=nt[j]; nt[j]=itmp;
      i++; j--;
    }
  } while(i<=j);

  if(lo<j) quicksrt(x,lo,j,mt,nt);
  if(i<hi) quicksrt(x,i,hi,mt,nt);
  return 0;
}
