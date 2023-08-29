#include<stdio.h>
#include<math.h>
#include<stdlib.h>

int Lx,Ly,N_sites;
int J=1;
int state_arr[5000];
int nbrarr[5000][4];
float T;
int L;

void init_state()
{
    for(int i=0;i<N_sites;i++)
    {
        state_arr[i]=0;
    }
}

void init_ordered_state()
{
    for(int i=0;i<Lx;i++)
    {
        for(int j=0;j<Ly;j++)
        {
            state_arr[i*Ly+j]=1;
        }
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
        m+=fabs(rho_y/Lx - 0.5);
        //printf("%f\n",abs(-2));
    }
    return 2*m/Ly;
}

float H()
{
    float E=0;
    for(int i=0;i<N_sites;i++)
    {
        for(int j=0;j<4;j++)
            E+=state_arr[i]*state_arr[nbrarr[i][j]];
    }
    return -J*E;
}


float MC_update(float E)
{
    float beta=1/T;
    int i=(int)((float)rand()/RAND_MAX*N_sites);
    if(state_arr[i])
    {
        int j=(int)((float)rand()/RAND_MAX*4);
        if(!state_arr[nbrarr[i][j]])
        {
            float delE = 0;
            for(int k=0;k<4;k++)
                delE+= state_arr[nbrarr[i][k]] - state_arr[nbrarr[nbrarr[i][j]][k]];
            delE=2*J*(delE+1);

            if(delE<=0 || ((float)rand()/RAND_MAX)<exp(-beta*delE))
            {
                E+=delE;
                state_arr[i]=0;
                state_arr[nbrarr[i][j]]=1;
            }
        }
    }
    return E;
}



void get_data()
{
    int Trlax=10000000;
    float E=H();
    float m_mean=0;
    float m2_mean=0;
    float m4_mean=0;
    int N_steps=10000000;
    srand(99820);
    for(int i=0;i<Trlax;i++) //reach steady state
        for(int j=0;j<N_sites;j++)
            E=MC_update(E);
    printf("steady state reached, ord param =%f\n",ord_param());
    for(int i=0;i<N_steps;i++)
    {
        for(int j=0;j<N_sites;j++)
            E=MC_update(E);
        //printf("MC step %d\n",i);
        float m=ord_param();
        m_mean+=m;
        m2_mean+=m*m;
        m4_mean+=m*m*m*m;
    }
    float data_arr[]={m_mean/N_steps,m2_mean/N_steps,m4_mean/N_steps};
    printf("m=%f\nm2=%f\nm4=%f",data_arr[0],data_arr[1],data_arr[2]);
    //return data_arr;
}

int main()
{
    scanf("%d %f",&L,&T);
    Lx=L;Ly=2*L;N_sites=Lx*Ly;
    nbr2D(nbrarr);

    init_state();
    init_randomize_state();

    printf("ord param of initial state = %f\n",ord_param());
    get_data();
    //float* data_arr=get_data();

    //printf("m_mean = %f\nm2_mean = %f\nm4_mean = %f",data_arr[0],data_arr[1],data_arr[2]);
    return 0;
}


