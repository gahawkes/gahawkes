library(alabama)
MultiHawkesLikelihood<-function(theta,points,process_dim)
{

    #mu:background intensity M x 1 matrix
    #alpha,beta: M x M matrix, with rows i being the parameters of process i
    #Get parameters
    M<-process_dim
    mu<-matrix(theta[1:M],M,1);
    #print(mu)
    alpha<-matrix(theta[(M+1):(M^2+M)],M,M);
    #beta<-array(c(1,1,1,1),dim=c(2,2));
    beta<-matrix(theta[(M^2+M+1):(2*M^2+M)],M,M);

    #Get datas
    #points:dataset,type list, with list 1 being the total time(Nx1 matrix), and list2 being a list of the time of subprocesses
    t<-points[[1]]; 
    sub_process<-list();
    #print('test')
    for(m in 1:M)
    {
        sub_process[[m]]<-points[[2]][[m]];
    }
    #datanum_sub: M x 1 matrix to store the length of total process and sub processes
    #datanum_total: 1 x 1 matrix to store the length of total process
    datanum_sub<-matrix(0,M,1);
    for(m in 1:M)
    {
        datanum_sub[m,1]<-dim(sub_process[[m]])[1];
    }
    datanum_total<-dim(t)[1];
    
    #calculate R cubic
    dimR<-max(datanum_sub);
    R<-array(NaN,c(M,M,dimR));
    #Set the first floor of R as zero
    R[,,1]<-0;
    #fill the R matrix
    for(m in 1:M)
    {
        for(n in 1:M)
        {
            if(m!=n)
            {
                for(i in 2:datanum_sub[m,1])
                {
                    temp_sum<-0;
                    for(k in 1:datanum_sub[n,1])
                    {
                        if(sub_process[[n]][k,1]>=sub_process[[m]][i,1]){break;}
                        if((sub_process[[n]][k,1]>=sub_process[[m]][i-1,1])&&(sub_process[[n]][k,1]<sub_process[[m]][i,1]))
                        {
                            temp_sum<-temp_sum+exp(-beta[m,n]*(sub_process[[m]][i,1]-sub_process[[n]][k,1]));
                        }
                    }
                    R[m,n,i]<-exp(-beta[m,n]*(sub_process[[m]][i,1]-sub_process[[m]][i-1,1]))*R[m,n,i-1]+temp_sum;
                }
            }
            else
            {
                for(i in 2:datanum_sub[m,1])
                {
                    R[m,n,i]<-exp(-beta[m,n]*(sub_process[[m]][i,1]-sub_process[[m]][i-1,1]))*(R[m,n,i-1]+1);
                }
            }
        }
    }
    #R filled
    #Begin calculate likelihood for each coordinate
    #sub_le is a M x 1 matrix to store the likelihood of each coordinate
    
    sub_le<-matrix(0,M,1);
    T<-t[datanum_total,1];

    for(m in 1:M)
    {
        part1<-0;
        #for(i in 1:datanum_sub[m])#datanum_total
        #{
        #    for(n in 1:M)
        #    {
        #        part1<-part1+alpha[m,n]/beta[m,n]*(1-exp(-beta[m,n]*(T-sub_process[[m]][i,1])));#t[i,1]
        #    }
        #}
        for(n in 1:M)#datanum_total
        {
            for(i in 1:datanum_sub[n])
            {
                part1<-part1+alpha[m,n]/beta[m,n]*(1-exp(-beta[m,n]*(T-sub_process[[n]][i,1])));#t[i,1]
            }
        }

        part2<-0;
        for(i in 2:datanum_sub[m])
        {
            temp_sum<-0;
            for(n in 1:M)
            {
                temp_sum<-temp_sum+alpha[m,n]*R[m,n,i];
            }
            part2<-part2+log(mu[m,1]+temp_sum);            
        }
        #print(part2)
        sub_le[m,1]<--mu[m,1]*T-part1+part2;
    }

    return (sum(sub_le));

}
MLEconstrain<-function(theta,points,process_dim)
{
    M<-process_dim
    gamma<-theta[(M+1):(M^2+M)]/theta[(M^2+M+1):(2*M^2+M)];
    return (gamma);
}
MLEhin<-function(theta,points,process_dim)
{

    M<-process_dim
    mu<-theta[1:M];
    alpha<-theta[(M+1):(M^2+M)];
    #beta<-array(c(1,1,1,1),dim=c(2,2));
    beta<-theta[(M^2+M+1):(2*M^2+M)]
    gamma<-1-theta[(M+1):(M^2+M)]/theta[(M^2+M+1):(2*M^2+M)];
    #mu_lb<-mu-0.0001
    #mu_ub<-10-mu
    #alpha_lb<-alpha-0.0001
    #alpha_ub<-10-alpha;
    #beta_ub<-10-beta
    #beta_lb<-beta-0.0001
    h<-c(gamma,mu,alpha,beta)
    h;
}


#Optimization
#Benchmark
loop_num<-20
result<-matrix(0,1,10)

mu<-array(c(0.5,0.5),dim=c(2,1));
alpha<-array(c(0.2,0.5,0.2,0.5),dim=c(2,2));
beta<-array(c(1,1,1,1),dim=c(2,2));
#theta<-c(mu,alpha)
theta<-c(mu,alpha,beta)
T<-10;
M<-2;
set.seed(1)
#for(i in 1:loop_num)
#{
  #simulation
  #cat('Start simulation ',i,'\n')
  start<-proc.time()
  #x<-MultiHawkesSim(mu,alpha,beta,T);
    result1<-MultiHawkesLikelihood(theta,x,2)
    #print(result1)
  #theta<-runif(10)
  par<-constrOptim.nl(par=theta,fn=MultiHawkesLikelihood,hin=MLEhin,control.optim=list(trace=0,fnscale=-0.01),points=x,process_dim=M);
  #print(par$par)
  #write(par$par,file='pars_100_changed.csv',ncolumns=10,sep=",",append=TRUE);
  #print(proc.time()-start)
  #cat('End simulation ',i,'\n')
#}
result<-result[-1,]