#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int Lx,Ly,N_sites;
int J=1;
int state_arr[5000];

void init_state()
{
    for(int i=0;i<N_sites;i++)
    {
        state_arr[i]=0;
    }
}
void init_randomize_state()
{
    srand(20);
    for(int i=0;i<N_sites/2;i++)
    {
        int pos=(int)((float)rand()/RAND_MAX*N_sites);
        //printf("%d ",pos);
        state_arr[pos]=1;
    }

}

void nbr2D(int nbrarr[5000][4])
{
    for(int i=0;i<N_sites;i++)
    {
        int k1=i/Ly;
        int k2=i%Ly;

        if(k1>=1 && k1<=Lx-2)
        {
            nbrarr[i][2]=i-Ly;
            nbrarr[i][3]=i+Ly;
        }
        else
        {
            int a=(k1==0);
            nbrarr[i][2]=Ly*(Lx-2+a)+k2;
            nbrarr[i][3]=Ly*a+k2;
        }
        if(k2>=1 && k2<=Ly-2)
        {
            nbrarr[i][0]=i-1;
            nbrarr[i][1]=i+1;
        }
        else
        {
            int b=(k2==0);
            nbrarr[i][0]=(i-1+Ly*b);
            nbrarr[i][1]=(i+1-Ly*(1-b));
        }
    }

}

float ord_param()
{
    float m=0;
    for(int j=0;j<Ly;j++)
    {
        float rho_y=0;
        for(int i=0;i<Lx;i++)
        {
            rho_y+=state_arr[i*Ly+j];
        }
        m+=abs(rho_y/Lx - 0.5);
    }
    return 2*m/Ly;
}

int main()
{
    int L;float T;
    scanf("%d",&L);
    Lx=L;Ly=2*L;N_sites=Lx*Ly;

    int nbrarr[N_sites][4];
    nbr2D(nbrarr);
    init_state();
    init_randomize_state();



    return 0;
}


