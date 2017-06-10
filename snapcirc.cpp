/**********************************************************
  名  称:SnapCircuit
  开发者:Dedcirc
  版  本:1.0.0
  描  述:基于复序网矩阵算法的电力系统故障分析程序。支持多线程运算加速
  (唯一指定邮箱:ionsuite@gmail.com)
 *********************************************************/
#include <stdio.h>
#include <math.h>
#include <libgen.h>
#include <unistd.h>
#include <pthread.h>
#define M 100                                /* 最大支路数 */
#define N 100                                /* 最大节点数 */
int F,T,R;
char C;
int w,h,r,
    n,                                       /* 节点数 */
    m,                                       /* 支路数 */
    v,
    LX;                                      /* 故障类型 */
float b1[M][M],b0[M][M];                     /* 电阻导纳的数组 */
float rr,tt,C1[M][M],C0[M][M];
static float Z1[M][M],Z0[M][M];              /* 阻抗矩阵的数组 */
static float Ip1,Ip2,Ip0;                    /* 口电流 */
static float Up1,Up2,Up0;                    /* 口电压 */
static float Zp1,Zp2,Zp0;                    /* 口阻抗 */
static float U1[M],U2[M],U0[M],il[N],I1[M][M],I2[M][M],I0[M][M],UAR[M],UBR[M],UCR[M],UAI[M],UBI[M],UCI[M],IAR[M][M],IBR[M][M],ICR[M][M],IAI[M][M],IBI[M][M],ICI[M][M],UA[M],UB[M],UC[M],IA[M][M],IB[M][M],IC[M][M];
struct zlz                                   /* 正序结构体 */
{
  int h;                                     /* h=1为正序 */
  int p1,p2;                                 /* 支路前后的两个节点 */
  float x;                                   /* 支路的电抗 */
} zlz[M];
struct zlf                                   /* 零序结构体 */
{
  int h;                                     /* h=0为零序 */
  int p1,p2;                                 /* 支路前后的两个节点 */
  float x;                                   /* 支路的电抗 */
} zlf[M];
struct sdl                                   /* 节点输入电流结构体 */
{
  int h,                                     /* h=3为节点注入电流 */
      p1;                                    /* 注入的节点 */
  float i;                                   /* 注入的电流值 */
} sdl[M];
FILE *fp1,*fp2;                              /* 文件指针定义*/
pthread_t thread[2];                         /* 线程申请 */

void Read_data(int argc,char **argv)         /* 读取电网特征文件 */
{
  int i,j;
  chdir(dirname(argv[0]));
  fp1=fopen("./input.txt","r");              /* 打开 input.txt 文件 */
  if(fp1==NULL)
  {
    printf("error:届不到,届不到\n");           /*找不到input.txt或内容是空的*/
  }
  fscanf(fp1,"%d,%d,%d,%d\n",&n,&m,&v,&r);   /* 输入节点数,负支路数,注入电流数 */
  j=1;
  do
  {
    for(i=1; i<=m+v+r; i++)
    {
      fscanf(fp1,"%d",&h);
      if(h==1)
      {
        fscanf(fp1,",%d,%d,%f\n",&zlz[i].p1,&zlz[i].p2,&zlz[i].x);
        j++;
      }
      if(h==0)
      {
        fscanf(fp1,",%d,%d,%f\n",&zlf[i].p1,&zlf[i].p2,&zlf[i].x);
        j++;
      }
      if(h==3)
      {
        fscanf(fp1,",%d,%f\n",&sdl[i].p1,&sdl[i].i);
        j++;
      }
    }
  }
  while(j<=m+v+r);
  fclose(fp1);
  if((fp2=fopen("output.txt","w"))==NULL)
  {
    printf("error:届不到,届不到\n");
  }
  /* 程序原始数据输出 */
  fprintf(fp2,"\n**********************原始数据**********************\n\n");
  fprintf(fp2,"节点数:%2d\t正支路数:%2d\t负支路数:%2d\t注入电流数:%2d\n\n",n,m,v,r);
  for(i=1; i<=m; i++)
  {
    fprintf(fp2,"相关节点:%2d,%2d\tX1=%f \n",zlz[i].p1,zlz[i].p2,zlz[i].x);
  }
  for(i=m+1; i<=m+v; i++)
  {
    fprintf(fp2,"相关节点:%2d,%2d\tX0=%f \n",zlf[i].p1,zlf[i].p2,zlf[i].x);
  }
  for(i=m+v+1; i<=m+v+r; i++)
  {
    fprintf(fp2,"注入电流节点:%2d\til=%f \n",sdl[i].p1,sdl[i].i);
  }
  fprintf(fp2,"\n**********************计算结果**********************\n");
}

void Form()                                  /* 形成节点矩阵并输出 */
{
  int a,b,i,j,k;
  float Z;
  float ZZ[50][50];
  for(i=0;i<=m;i++)
  {
    for(j=0;j<=m;j++)
    {
      b1[i][j]=0;b0[i][j]=0;
      Z1[i][j]=0;Z0[i][j]=0;
    }
  }
  for(i=1;i<=m;i++)
  {
    a=zlz[i].p1;
    b=zlz[i].p2;
    b1[a][b]=1/zlz[i].x;
    b1[b][a]=1/zlz[i].x;
  }
  for(i=m+1;i<=m+v;i++)
  {
    a=zlf[i].p1;
    b=zlf[i].p2;
    b0[a][b]=1/zlf[i].x;
    b0[b][a]=1/zlf[i].x;
  }
  /* 根据节点电压法生成导纳矩阵 */
  /* 正序 */
  for(i=1;i<=n;i++)
  for(j=1;j<=n;j++)
  {
    Z1[i][j]=0;
  }
  for(i=1;i<=n;i++)
  for(j=0;j<=n;j++)
  {
    Z1[i][i]=Z1[i][i]+b1[i][j];
  }
  for(i=1;i<=n;i++)
  for(j=1;j<=n;j++)
  if(i!=j)
  {
    Z1[i][j]=-1*b1[i][j];
  }
  /* 零序 */
  for(i=1;i<=n;i++)
  for(j=1;j<=n;j++)
  {
    Z0[i][j]=0;
  }
  for(i=1;i<=n;i++)
  for(j=0;j<=n;j++)
  {
    Z0[i][i]=Z0[i][i]+b0[i][j];
  }
  for(i=1;i<=n;i++)
  for(j=1;j<=n;j++)
  if(i!=j)
  {
    Z0[i][j]=-1*b0[i][j];
  }
  fprintf(fp2,"\n正序导纳矩阵:");
  for(i=1;i<=n;i++)
  {
    fprintf(fp1,"\n");
    for(j=1;j<=n;j++)
    {
      if(Z1[i][j]>=0)
      fprintf(fp2,"j%f\t",Z1[i][j]);
      else
      fprintf(fp2,"-j%f\t",-Z1[i][j]);
    }
  }
  fprintf(fp2,"\n\n零序导纳矩阵:");
  for(i=1;i<=n;i++)
  {
    fprintf(fp1,"\n");
    for(j=1;j<=n;j++)
    {
      if(Z0[i][j]>=0)
      fprintf(fp2,"j%f\t",Z0[i][j]);
      else
      fprintf(fp2,"-j%f\t",-Z0[i][j]);
    }
  }
  /* 生成正序阻抗矩阵 */
  for(i=1;i<=n;i++)
  for(j=1;j<=n;j++)
  {
    ZZ[i][j]=0;
  }
  for(k=1;k<=n;k++)
  {
    ZZ[k][k]=-1/Z1[k][k];
    for(i=1;i<=n;i++)
    for(j=1;j<=n;j++)
    {
      if(i!=k&&j!=k)
      {
        ZZ[i][k]=-Z1[i][k]/Z1[k][k];
        ZZ[k][j]=-Z1[k][j]/Z1[k][k];
        Z=Z1[i][k]*Z1[k][j]/Z1[k][k];
        ZZ[i][j]=Z1[i][j]-Z;
      }
    }
    for(i=1;i<=n;i++)
    for(j=1;j<=n;j++)
    {
      Z1[i][j]=ZZ[i][j];
    }
  }
  for(i=1;i<=n;i++)
  for(j=1;j<=n;j++)
  {
    Z1[i][j]=-Z1[i][j];
  }
  /* 生成零序阻抗矩阵 */
  for(i=1;i<=n;i++)
  for(j=1;j<=n;j++)
  {
    ZZ[i][j]=0;
  }
  for(k=1;k<=n;k++)
  {
    ZZ[k][k]=-1/Z0[k][k];
    for(i=1;i<=n;i++)
    for(j=1;j<=n;j++)
    {
      if(i!=k&&j!=k)
      {
        ZZ[i][k]=-Z0[i][k]/Z0[k][k];
        ZZ[k][j]=-Z0[k][j]/Z0[k][k];
        Z=Z0[i][k]*Z0[k][j]/Z0[k][k];
        ZZ[i][j]=Z0[i][j]-Z;
      }
    }
    for(i=1;i<=n;i++)
    for(j=1;j<=n;j++)
    {
      Z0[i][j]=ZZ[i][j];
    }
  }
  for(i=1;i<=n;i++)
  for(j=1;j<=n;j++)
  {
    Z0[i][j]=-Z0[i][j];
  }
  fprintf(fp2,"\n\n正序阻抗矩阵:");
  for(i=1;i<=n;i++)
  {
    fprintf(fp1,"\n");
    for(j=1;j<=n;j++)
    {
      if(Z1[i][j]>=0)
      fprintf(fp2,"j%f\t",Z1[i][j]);
      else
      fprintf(fp2,"-j%f\t",-Z1[i][j]);
    }
  }
  fprintf(fp2,"\n\n零序阻抗矩阵:");
  for(i=1;i<=n;i++)
  {
    fprintf(fp1,"\n");
    for(j=1;j<=n;j++)
    {
      if(Z0[i][j]>=0)
      fprintf(fp2,"j%f\t",Z0[i][j]);
      else
      fprintf(fp2,"-j%f\t",-Z0[i][j]);
    }
  }
  fprintf(fp2,"\n");
  for(i=1;i<=n;i++)
  il[i]=0;
  for(i=m+v+1;i<=m+v+r;i++)
  {
    a=sdl[i].p1;
    il[a]=sdl[i].i;
  }
  fprintf(fp2,"\n注入电流值:\n");
  for(i=1;i<=n;i++)
  fprintf(fp2,"il[%d]=%f\n",i,il[i]);
}

void REVISION()                              /* 断相阻抗矩阵修正 */
{
  int i,j;
  printf("请输入发生两相断相故障的2个节点(F&T):");
  scanf("%d,%d",&F,&T);
  for(i=1;i<=m;i++)
  for(j=1;i<=m;i++)
  {
    Z1[i][j]=Z1[i][j]-(Z1[i][F]-Z1[i][T])*(Z1[F][j]-Z1[T][j])/(Z1[F][F]+Z1[T][T]-2*Z1[F][T]-1/b1[F][T]);
  }
  for(i=1;i<=m;i++)
  for(i=1;i<=m;i++)
  {
    Z0[i][j]=Z0[i][j]-(Z0[i][F]-Z0[i][T])*(Z0[F][j]-Z0[T][j])/(Z0[F][F]+Z0[T][T]-2*Z0[F][T]-1/b0[F][T]);
  }
}

void TRY()                                   /* 三相短路故障 */
{
  int j;
  float tt1;
  printf("请输入发生三相故障的节点:");
  scanf("%d",&F);
  tt1=0;
  for(j=1;j<=n;j++)
  {
    tt1+=Z1[F][j]*il[j];
  }
  Zp1=Z1[F][F];
  Ip1=-tt1/Zp1;
  Ip2=0;
  Ip0=0;
  Up1=0;
  Up2=0;
  Up0=0;
}

void SERI()                                  /* 串联型故障 */
{
  int j;
  float zs,tt1;
  if(LX==1)
  {
    printf("请输入发生一相接地短路故障的节点:");
    scanf("%d",&F);
    printf("请输入接地阻抗zs:");
    scanf("%f",&zs);
    tt1=0;
    for(j=1;j<=n;j++)
    {
      tt1+=Z1[F][j]*il[j];
    }
    Zp1=Z1[F][F];
    Zp2=Z1[F][F];
    Zp0=Z0[F][F];
    Ip1=-tt1/(Zp1+Zp2+Zp0+3*zs);
    Ip2=Ip1;
    Ip0=Ip1;
    Up1=tt1+Ip1*Zp1;
    Up2=Ip2*Zp2;
    Up0=Ip0*Zp0;
  }
  if(LX==2)
  {
    Zp1=Z1[F][F]+Z1[T][T]-2*Z1[F][T];
    Zp2=Z1[F][F]+Z1[T][T]-2*Z1[F][T];
    Zp0=Z0[F][F]+Z0[T][T]-2*Z0[F][T];
    tt1=0;
    for(j=1;j<=n;j++)
    {
      tt1+=Z1[F][j]*il[j];
    }
    Ip1=-tt1/(Zp1+Zp2+Zp0);
    Ip2=Ip1;
    Ip0=Ip1;
    Up1=Ip1*(Zp2+Zp0);
    Up2=-Ip2*Zp2;
    Up0=-Ip0*Zp0;
  }
}

void PARA()                                  /* 并联型故障 */
{
  int j;
  float tt1,zp;
  if(LX==3||LX==4)
  {
    printf("请输入发生两相(接地)短路的节点:");
    scanf("%d",&F);
    printf("请输入接地阻抗zp:");
    scanf("%f",&zp);
    Zp1=Z1[F][F];
    Zp2=Z1[F][F];
    Zp0=Z0[F][F];
    if(LX==3)
    {
      tt1=0;
      for(j=1;j<=n;j++)
      {
        tt1+=Z1[F][j]*il[j];
      }
      Ip1=-tt1/(2*Zp1+zp);
      Ip2=-Ip1;
      Ip0=0;
      Up1=tt1+Z1[F][F]*Ip1;
      Up2=Z1[F][F]*Ip2;
      Up0=0;
    }
    if(LX==4)
    {
      tt1=0;
      for(j=1;j<=n;j++)
      {
        tt1+=Z1[F][j]*il[j];
      }
      Ip1=-tt1/(Zp1+Zp1*(Zp0+zp)/(Zp1+Zp0+3*zp));
      Ip2=-Ip1*(Zp0+3*zp)/(Zp1+Zp0+3*zp);
      Ip0=-Ip1*Zp1/(Zp1+Zp0+3*zp);
      Up1=tt1+Z1[F][F]*Ip1;
      Up2=Z1[F][F]*Ip2;
      Up0=Z0[F][F]*Ip0;
    }
  }
  if(LX==5)
  {
    Zp1=Z1[F][F]+Z1[T][T]-2*Z1[F][T];
    Zp2=Z1[F][F]+Z1[T][T]-2*Z1[F][T];
    Zp0=Z0[F][F]+Z0[T][T]-2*Z0[F][T];
    tt1=0;
    for(j=1;j<=n;j++)
    {
      tt1+=Z1[F][j]*il[j];
    }
    Ip1=-tt1/(Zp1+Zp2*Zp0/(Zp2+Zp0));
    Ip2=-Ip1*Zp0/(Zp2+Zp0);
    Ip0=-Ip1*Zp2/(Zp2+Zp0);
    Up1=Ip1*Zp2*Zp0/(Zp2+Zp0);
    Up2=Up1;
    Up0=Up1;
  }
}

void VOLT()                                  /* 计算节点电压 */
{
  int i,j;
  float tt[M];
  for(i=1;i<=n;i++)
  tt[i]=0;
  for(i=1;i<=n;i++)
  for(j=1;j<=n;j++)
  {
    tt[i]=tt[i]+Z1[i][j]*il[j];
  }
  if(LX==0)
  {
    for(i=1;i<=n;i++)
    {
      U1[i]=tt[i]+Z1[i][F]*Ip1;
      U2[i]=0;
      U0[i]=0;
    }
  }
  if(LX==1||LX==2||LX==3||LX==4||LX==5)
  {
    for(i=1;i<=n;i++)
    {
      U1[i]=tt[i]+Z1[i][F]*Ip1;
      U2[i]=Z1[i][F]*Ip2;
      U0[i]=Z0[i][F]*Ip0;
    }
  }
}

void *CURRENT(void *)                        /* 计算各支路电流 */
{
  int a,b,i,j;
  for(i=1;i<=n;i++)
  for(j=1;j<=n;j++)
  {
    I1[i][j]=0;
    I2[i][j]=0;
    I0[i][j]=0;
  }
  for(i=1;i<=m;i++)
  {
    a=zlz[i].p1;
    b=zlz[i].p2;
    I1[a][b]=(U1[a]-U1[b])/zlz[i].x;
    I2[a][b]=(U2[a]-U2[b])/zlz[i].x;
    a=zlf[m+i].p1;b=zlf[m+i].p2;
    I0[a][b]=(U0[a]-U0[b])/zlf[m+i].x;
  }
  pthread_exit(NULL);
}

void BRANCH()                                /* 计算分支系数 */
{
  int i,a,b;
  for(i=1;i<=m;i++)
  {
    a=zlz[i].p1;
    b=zlz[i].p2;
    if(a!=0&&b!=0)
    C1[a][b]=I1[a][b]/Ip1;
    C0[a][b]=I0[a][b]/Ip0;
  }
}

void *PHASE(void *)                          /* 换相 */
{
  int i,a,b;
  for(i=1;i<=n;i++)
  {
    UAR[i]=0;
    UBR[i]=0.866*(U1[i]-U2[i]);
    UCR[i]=0.866*(U2[i]-U1[i]);
    UAI[i]=U0[i]+U1[i]+U2[i];
    UBI[i]=U0[i]-0.5*(U1[i]+U2[i]);
    UCI[i]=U0[i]-0.5*(U1[i]+U2[i]);
  }
  for(i=1;i<=m;i++)
  {
    a=zlz[i].p1;
    b=zlz[i].p2;
    IAR[a][b]=I0[a][b]+I1[a][b]+I2[a][b];
    IBR[a][b]=I0[a][b]-0.5*(I1[a][b]+I2[a][b]);
    ICR[a][b]=I0[a][b]-0.5*(I1[a][b]+I2[a][b]);
    IAI[a][b]=0;
    IBI[a][b]=0.866*(I2[a][b]-I1[a][b]);
    ICI[a][b]=0.866*(I1[a][b]-I2[a][b]);
  }
  pthread_exit(NULL);
}

void MODULE()                                /* 计算相电压/电流的模值 */
{
  int i,a,b;
  for(i=1;i<=n;i++)
  {
    UA[i]=sqrt(UAR[i]*UAR[i]+UAI[i]*UAI[i]);
    UB[i]=sqrt(UBR[i]*UBR[i]+UBI[i]*UBI[i]);
    UC[i]=sqrt(UCR[i]*UCR[i]+UCI[i]*UCI[i]);
  }
  for(i=1;i<=m;i++)
  {
    a=zlz[i].p1;b=zlz[i].p2;
    IA[a][b]=sqrt(IAR[a][b]*IAR[a][b]+IAI[a][b]*IAI[a][b]);
    IB[a][b]=sqrt(IBR[a][b]*IBR[a][b]+IBI[a][b]*IBI[a][b]);
    IC[a][b]=sqrt(ICR[a][b]*ICR[a][b]+ICI[a][b]*ICI[a][b]);
  }
}

void OUTLET()                                /* 输出结果到文件 */
{
  int i,a,b;
  fprintf(fp2,"\n故障类型:%d\n",LX);
  fprintf(fp2,"\n口参数值:\n");
  fprintf(fp2,"Ip1=%f\nIp2=%f\nIp0=%f\n",Ip1,Ip2,Ip0);
  fprintf(fp2,"Up1=%f\nUp2=%f\nUp0=%f\n",Up1,Up2,Up0);
  fprintf(fp2,"\n各节点序电压值:");
  for(i=1;i<=n;i++)
  {
    if(U1[i]>=0)
    fprintf(fp2,"\nU1[%d]=j%f\n",i,U1[i]);
    else
    fprintf(fp2,"\nU1[%d]=-j%f\n",i,-U1[i]);
    if(U2[i]>=0)
    fprintf(fp2,"U2[%d]=j%f\n",i,U2[i]);
    else
    fprintf(fp2,"U2[%d]=-j%f\n",i,-U2[i]);
    if(U0[i]>=0)
    fprintf(fp2,"U0[%d]=j%f\n",i,U0[i]);
    else
    fprintf(fp2,"U1[%d]=-j%f\n",i,-U0[i]);
  }
  fprintf(fp2,"\n各支路序电流值:");
  for(i=1;i<=m;i++)
  {
    a=zlz[i].p1;b=zlz[i].p2;
    fprintf(fp2,"\nI1[%d][%d]=%f\nI2[%d][%d]=%f\nI0[%d][%d]=%f\n",a,b,I1[a][b],a,b,I2[a][b],a,b,I0[a][b]);
  }
  fprintf(fp2,"\n分布系数:\n");
  for(i=1;i<=m;i++)
  {
    a=zlz[i].p1;b=zlz[i].p2;
    if(a!=0&&b!=0)
    {
      fprintf(fp2,"C1[%d][%d]=%f\nC0[%d][%d]=%f\n",a,b,C1[a][b],a,b,C0[a][b]);
    }
  }
  fprintf(fp2,"\n各节点各相电压值:");
  for(i=1;i<=n;i++)
  {
    if(UAI[i]>0)
    fprintf(fp2,"\nUA[%d]=%f+j%f\t",i,UAR[i],UAI[i]);
    if(UAI[i]<0)
    fprintf(fp2,"\nUA[%d]=%f-j%f\t",i,UAR[i],-UAI[i]);
    if(UAI[i]==0)
    fprintf(fp2,"\nUA[%d]=%f\t",i,UAR[i]);
    if(UBI[i]>0)
    fprintf(fp2,"UB[%d]=%f+j%f\t",i,UBR[i],UBI[i]);
    if(UBI[i]<0)
    fprintf(fp2,"UB[%d]=%f-j%f\t",i,UBR[i],-UBI[i]);
    if(UBI[i]==0)
    fprintf(fp2,"UB[%d]=%f\t",i,UBR[i]);
    if(UCI[i]>0)
    fprintf(fp2,"UC[%d]=%f+j%f\t",i,UCR[i],UCI[i]);
    if(UCI[i]<0)
    fprintf(fp2,"UC[%d]=%f-j%f\t",i,UCR[i],-UCI[i]);
    if(UCI[i]==0)
    fprintf(fp2,"UC[%d]=%f\t",i,UCR[i]);
  }
  fprintf(fp2,"\n\n各支路各相电流值:");
  for(i=1;i<=m;i++)
  {
    a=zlz[i].p1;b=zlz[i].p2;
    if(IAI[a][b]>0)
    fprintf(fp2,"\nIA[%d][%d]=%f+j%f\t",a,b,IAR[a][b],IAI[a][b]);
    if(IAI[a][b]<0)
    fprintf(fp2,"\nIA[%d][%d]=%f-j%f\t",a,b,IAR[a][b],-IAI[a][b]);
    if(IAI[a][b]==0)
    fprintf(fp2,"\nIA[%d][%d]=%f\t",a,b,IAR[a][b]);
    if(IBI[a][b]>0)
    fprintf(fp2,"IB[%d][%d]=%f+j%f\t",a,b,IBR[a][b],IBI[a][b]);
    if(IBI[a][b]<0)
    fprintf(fp2,"IB[%d][%d]=%f-j%f\t",a,b,IBR[a][b],-IBI[a][b]);
    if(IBI[a][b]==0)
    fprintf(fp2,"IB[%d][%d]=%f\t",a,b,IBR[a][b]);
    if(ICI[a][b]>0)
    fprintf(fp2,"IC[%d][%d]=%f+j%f\t",a,b,ICR[a][b],ICI[a][b]);
    if(ICI[a][b]<0)
    fprintf(fp2,"IC[%d][%d]=%f-j%f\t",a,b,ICR[a][b],-ICI[a][b]);
    if(ICI[a][b]==0)
    fprintf(fp2,"IC[%d][%d]=%f\t",a,b,ICR[a][b]);
  }
  fprintf(fp2,"\n\n各节点各相电压模值:");
  for(i=1;i<=n;i++)
  {
    fprintf(fp2,"\n|UA|[%d]=%f\t|UB|[%d]=%f\t|UC|[%d]=%f",i,UA[i],i,UB[i],i,UC[i]);
  }
  fprintf(fp2,"\n\n输出各支路各相电流模值:");
  for(i=1;i<=m;i++)
  {
    a=zlz[i].p1;b=zlz[i].p2;
    fprintf(fp2,"\n|IB|[%d][%d]=%f\t|IC|[%d][%d]=%f",a,b,IA[a][b],a,b,IB[a][b],a,b,IC[a][b]);
  }
  fprintf(fp2,"\n\n");
}

int main(int argc, char *argv[])
{
  printf("\n*************************************************************\n                                          _\n      ____                               /_/\n     /  __%c      _          _____ _____   _  _____  _____\n     %c  %c   / %c | |   /%c   | -__/ | ___%c | | %c   _/ | ___%c\n    __%c  %c  | |%c| |  / -%c  | |    | |___ | | | /    | |___\n    %c____/  |_| %c / / /%c %c %c/     |____/ |_%c |/     |____/\n  \n A BETTER ELECTICAL POWER SYSTEM ERROR CULCULATOR\n POWERED BY ionsuite IN 2017\n\n*************************************************************\n\n",92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92,92);
  Read_data(argc,argv);
  Form();
  for(w=1;w<=100;w++)
  {
    printf("0.三相短路\n1.单相接地\n2.两相断相\n3.两相短路\n4.两相接地短路\n5.单相断相\n请选择故障类型:");
    scanf("%d",&LX);
    if(LX==0)
    TRY();
    if(LX==1||LX==2)
    {
      if(LX==2)
      REVISION();
      SERI();
    }
    if(LX==3||LX==4||LX==5)
    {
      if(LX==5)
      REVISION();
      PARA();
    }
    pthread_create(&thread[0],NULL,CURRENT,(void *)(NULL));
    VOLT();
    pthread_join(thread[0],NULL);
    pthread_create(&thread[1],NULL,PHASE,(void *)(NULL));
    BRANCH();
    pthread_join(thread[1],NULL);
    MODULE();
    OUTLET();
    printf("是否存在其它的故障(Y/N):");
    getchar();
    scanf("%c",&C);
    if(C=='Y'||C=='y')
    continue;
    else
    break;
  }
  printf("\n*************************************************************\n\n原始数据由input.txt输入,计算结果在output.txt中\n");
  return 0;
}
