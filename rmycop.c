# include <string.h>
void RMYCOPSTRITE(char *out,const char *in,int *n,int *outstr,int *instr)
{
	memcpy(out+*(outstr),in+(*instr),*n);
	/*bcopy(in+(*instr),out+*(outstr),*n);*/
}
