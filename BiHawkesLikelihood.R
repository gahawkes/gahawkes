library(alabama)
library(multicore)
BiHawkesLikelihood<-function(theta,points)
{

    mu_1<-theta[1];mu_2<-theta[2];
    alpha_11<-theta[3];alpha_12<-theta[5];alpha_21<-theta[4];alpha_22<-theta[6];
    beta_1<-theta[7];beta_2<-theta[8]    
    #points:dataset,type list, with list 1 being the total time(Nx1 matrix), and list2 being a list of the time of subprocesses
    t<-points[[1]]; 
    #sub_process<-list();
    t1<-points[[2]][[1]];
    N<-dim(t1)[1];
    t2<-points[[2]][[2]];
    #print(t2)
    M<-dim(t2)[1];
    datanum_total<-dim(t)[1];

    #calculate R cubic
    R_11<-matrix(NaN,N,1);
    R_12<-matrix(NaN,N,1);
    R_21<-matrix(NaN,M,1);
    R_22<-matrix(NaN,M,1);
    
    #Set the first floor of R as zero
    R_11[1,1]<-0;R_12[1,1]<-0;R_21[1,1]<-0;R_22[1,1]<-0;
    #fill the R matrix
    for(i in 2:N)
    {
        R_11[i,1]<-(1+R_11[i-1,1])*(exp(-beta_1*(t1[i,1]-t1[i-1,1])));   
    }

    for(i in 2:N)
    {
        temp_sum<-0;
        for(j in 1:M)
        {
            if(t2[j,1]>=t1[i,1]){break;}
            if(t2[j,1]>=t1[i-1,1]&&t2[j,1]<t1[i,1]){temp_sum<-temp_sum+exp(-beta_1*(t1[i,1]-t2[j,1]));}
        }
        #print(R_12[i-1,1]*exp(-beta_1*(t1[i,1]-t1[i-1,1])))
        R_12[i,1]<-R_12[i-1,1]*exp(-beta_1*(t1[i,1]-t1[i-1,1]))+temp_sum;
    }

    for(j in 2:M)
    {
        R_22[j,1]<-(1+R_22[j-1,1])*(exp(-beta_2*(t2[j,1]-t2[j-1,1])));
    }

    for(j in 2:M)
    {
        temp_sum<-0;
        for(i in 1:N)
        {
            if(t1[i,1]>=t2[j,1]){break;}
            if(t1[i,1]>=t2[j-1,1]&&t1[i,1]<t2[j,1]){temp_sum<-temp_sum+exp(-beta_2*(t2[j,1]-t1[i,1]));}
        }
        R_21[j,1]<-R_21[j-1,1]*exp(-beta_2*(t2[j,1]-t2[j-1,1]))+temp_sum;
    }
    
    #R filled
    #Begin calculate likelihood for each coordinate
    
    T<-t[datanum_total,1];

    #Calculate L1
    part1<-0;
    part2<-0;
    part3<-0;
    #part1
    for(i in 2:N)
    {
        part1<-part1+log(mu_1+alpha_11*R_11[i,1]+alpha_12*R_12[i,1]);
    }
    #print(part1)
    #part2
    for(i in 1:N)
    {
        part2<-part2+(1-exp(-beta_1*(T-t1[i,1])));
    }
    part2<-part2*alpha_11/beta_1;
    #part3
    for(j in 1:M)
    {
        part3<-part3+(1-exp(-beta_1*(T-t2[j,1])));
    }
    part3<-part3*alpha_12/beta_1;
    L1<-part1-mu_1*T-part2-part3;
    
    #Calculate L2
    part1<-0;
    part2<-0;
    part3<-0;
    #part1
    for(j in 2:M)
    {
        part1<-part1+log(mu_2+alpha_21*R_21[j,1]+alpha_22*R_22[j,1]);
    }
    #print(part1)
    #part2
    for(j in 1:M)
    {
        part2<-part2+(1-exp(-beta_2*(T-t2[j,1])));
    }
    part2<-part2*alpha_22/beta_2;
    #part3
    for(i in 1:N)
    {
        part3<-part3+(1-exp(-beta_2*(T-t1[i,1])));
    }
    part3<-part3*alpha_21/beta_2;
    L2<-part1-mu_2*T-part2-part3;
    L<-L1+L2

    return (L);
}
#MLEconstrain<-function(theta,points,process_dim)
#{
#    M<-process_dim
#    gamma<-theta[(M+1):(M^2+M)]/theta[(M^2+M+1):(2*M^2+M)];
#    return (gamma);
#}
MLEhin<-function(theta,points)
{
    mu_1<-theta[1];mu_2<-theta[2];
    alpha_11<-theta[3];alpha_12<-theta[5];alpha_21<-theta[4];alpha_22<-theta[6];
    beta_1<-theta[7];beta_2<-theta[8] 
    h<-c(mu_1,mu_2,alpha_11,alpha_12,alpha_21,alpha_22,beta_1,beta_2,1-alpha_11/beta_1,1-alpha_12/beta_1,1-alpha_21/beta_2,1-alpha_22/beta_2);
    h;
}


#Optimization
#Benchmark
loop_num<-1;
mu<-array(c(0.5,0.5),dim=c(2,1));
alpha<-array(c(0.2,0.5,0.2,0.5),dim=c(2,2));
beta<-array(c(1,1,1,1),dim=c(2,2));
beta_r<-array(c(1,1),dim=c(2,1));
#theta<-c(mu,alpha)
theta<-c(mu,alpha,beta_r)

T<-100;
for(i in 1:loop_num)
{
   #simulation
   #cat('Start simulation ',i,'\n')
   #start<-proc.time()
    #set.seed(1)
   x<-MultiHawkesSim(mu,alpha,beta,T);
   true<-BiHawkesLikelihood(theta,x)
    #print(result)
   par<-constrOptim.nl(par=theta,fn=BiHawkesLikelihood,hin=MLEhin,control.outer=list(trace=FALSE),control.optim=list(trace=FALSE,fnscale=-0.01),points=x);
    print(par)
    write(c(true,par$value,par$par),file='pars_100_bh_parallel_test.csv',ncolumns=10,sep=",",append=TRUE);
    print(proc.time()-start)
    cat('End simulation ',i,'\n')
}

# TestSimulation<-function(x)
# {
#     #set.seed(1)
#     T<-50;
#     mu<-array(c(0.5,0.5),dim=c(2,1));
#     alpha<-array(c(0.2,0.5,0.2,0.5),dim=c(2,2));
#     beta<-array(c(1,1,1,1),dim=c(2,2));
#     beta_r<-array(c(1,1),dim=c(2,1));
#     y<-MultiHawkesSim(mu,alpha,beta,T);
#     par<-constrOptim.nl(par=theta,fn=BiHawkesLikelihood,hin=MLEhin,control.outer=list(trace=FALSE),control.optim=list(trace=FALSE,fnscale=-0.01),points=y);
#     return (par);
# }
# start<-proc.time()
# for(i in 1:5)
# {
#     TestSimulation(1);
# }
# end<-proc.time()-start;
# print(end)
