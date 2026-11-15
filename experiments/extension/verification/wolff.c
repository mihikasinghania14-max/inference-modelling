// Wolff single-cluster 2D Ising at T_c on LxL torus. Records |M|, M^2, M^4 every k cluster flips.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdint.h>
static uint64_t s[2];
static inline uint64_t rotl(const uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline uint64_t next(void){uint64_t s0=s[0],s1=s[1],r=s0+s1;s1^=s0;s[0]=rotl(s0,55)^s1^(s1<<14);s[1]=rotl(s1,36);return r;}
static inline double U(void){return (next()>>11)*(1.0/9007199254740992.0);}
int main(int argc,char**argv){
  int L=atoi(argv[1]); long nmeas=atol(argv[2]); s[0]=atol(argv[3])*0x9E3779B97F4A7C15ULL+1; s[1]=0xDEADBEEF12345ULL^s[0];
  int n=L*L; double K=0.5*log(1+sqrt(2.0)); double padd=1-exp(-2*K);
  int *sp=malloc(n*sizeof(int)), *stack=malloc(n*sizeof(int));
  for(int i=0;i<n;i++) sp[i]= (U()<0.5)?1:-1;
  long tot=0; for(int w=0;w<20000;w++){ // warmup + cluster size estimate
    int seed=(int)(U()*n), sg=sp[seed], top=0, cs=0; sp[seed]=-sg; stack[top++]=seed;
    while(top){int x=stack[--top]; cs++; int xi=x%L, yi=x/L; int nb[4]={yi*L+(xi+1)%L, yi*L+(xi+L-1)%L, ((yi+1)%L)*L+xi, ((yi+L-1)%L)*L+xi};
      for(int j=0;j<4;j++){int y=nb[j]; if(sp[y]==sg && U()<padd){sp[y]=-sg; stack[top++]=y;}}}
    if(w>=10000) tot+=cs;}
  double meanC=tot/10000.0; int k=(int)ceil(n/meanC); if(k<1)k=1;
  double a1=0,a2=0,a4=0; // also write block means
  int nb=100; long per=nmeas/nb; double b1[100]={0},b2[100]={0},b4[100]={0};
  for(long t=0;t<nmeas;t++){
    for(int f=0;f<k;f++){int seed=(int)(U()*n), sg=sp[seed], top=0; sp[seed]=-sg; stack[top++]=seed;
      while(top){int x=stack[--top]; int xi=x%L, yi=x/L; int nb4[4]={yi*L+(xi+1)%L, yi*L+(xi+L-1)%L, ((yi+1)%L)*L+xi, ((yi+L-1)%L)*L+xi};
        for(int j=0;j<4;j++){int y=nb4[j]; if(sp[y]==sg && U()<padd){sp[y]=-sg; stack[top++]=y;}}}}
    long M=0; for(int i=0;i<n;i++) M+=sp[i]; double m=fabs((double)M/n);
    int b=t/per; if(b>=nb)b=nb-1; b1[b]+=m; b2[b]+=m*m; b4[b]+=m*m*m*m;
  }
  for(int b=0;b<nb;b++) printf("%.15g %.15g %.15g\n", b1[b]/per, b2[b]/per, b4[b]/per);
  fprintf(stderr,"L=%d k=%d meanC=%.1f\n",L,k,meanC);
  return 0;}
