LambdaValue<-function(time,history,sub_process,mu,alpha,beta)
{
    #value is 1xM matrix to store ambda value of each process
    M<-dim(mu)[1]
    value<-mu;
    history_before_time<-history[history<=time,];
    sub_process_before_time<-list();
    #trim the sub_process to the length before time
    #print(sub_process[[2]])
    for(m in 1:M)
    {
        #coerce as matrix
        temp<-sub_process[[m]][sub_process[[m]]<=time,];
        #print(temp)
        #print(length(temp))
        temp_dim<-length(temp);
        if(temp_dim==0)
        {
            sub_process_before_time[[m]]<-array(0,dim=c(1,1))
        }
        else
        {
            sub_process_before_time[[m]]<-array(temp,dim=c(temp_dim,1))
        }
        #sub_process[[m]]<-array(sub_process[[m]],dim=c(temp_dim,1))
        #print(sub_process[[m]])
        #print(sub_process_before_time[[m]])

    }

    if(time<history[1,1])
    {
        #if time is smaller than the first event, the lambda will be mu
        value<-mu;
        #print(value)
        return (value);
    }
    else
    {
        for(m in 1:M)
        {

            #for lambda m
            for(j in 1:M)
            {
                #for calculating stimulation of sub_process j 
                for(i in 1:dim(sub_process_before_time[[j]])[1])
                {

                    #print(sub_process_before_time[[j]])    
                    value[m,1]<-value[m,1]+alpha[m,j]*exp(-beta[m,j]*(time-sub_process_before_time[[j]][i,1]));
                }
                    
            }
        }
    return (value);
    }
}

MultiHawkesSim<-function(mu,alpha,beta,T)
{

    #mu:background intensity M x 1 matrix
    #alpha,beta: M x M matrix, with rows i being the parameters of process i
    #T: scalar, end time of simulation
    
    #M is the total number of processes in the model
    M<-dim(mu)[1];
    total_index<-1;sub_index<-array(1,dim=c(1,M));
    total_intensity<-sum(mu);

    #First Event
    s<--log(runif(1))/total_intensity;
    #CC print(s)
    #t:Store event time,col1=real time,col2=Index of process, from 1 to M,col3=Index of event within a process
    t<-array(0,dim=c(1,1));
    #sub_process stores the sub history of each process
    sub_process<-list();
    n0=0;
    for(m in 1:M)
    {
        sub_process[[m]]<-array(0,dim=c(1,1));
    }
    #print(s)
    if(s>T)
    {
        return (list(t,sub_process));
        #print("quit!")
    }

    #Attribution test
    randomD<-runif(1);
    if(mu[1,1]/total_intensity>=randomD)
    {
        #print(randomD)
        t[1,1]<-s;
        #print(t)
        sub_process[[1]][1,1]<-s;
        n0<-1;
    }
    else
    {
        for(j in 2:M)
        {
            if ((sum(mu[c(1:j-1),1])/total_intensity<randomD)&&(sum(mu[c(1:j),1])/total_intensity>=randomD))
            {
                t[1,1]<-s;
                #print(t)
                sub_process[[j]][1,1]<-s;
                n0<-j;
            }
        }
    }
    #CC print(t)

    #General Routine
    out_range<-FALSE;
    while(out_range==FALSE)
    {
        total_index<-total_index+1;
        sub_index[1,n0]<-sub_index[1,n0]+1;
        #print(sub_process)
        total_intensity<-sum(LambdaValue(t[total_index-1,1],t,sub_process,mu,alpha,beta))+sum(alpha[,n0]);

        s<-s-log(runif(1))/total_intensity;
        #print(s)
        #CC
        if(s>T){break;}

        while(out_range==FALSE)
        { 
            #Attribution - Rejection test
            value=LambdaValue(s,t,sub_process,mu,alpha,beta);
            #print(value)
            intensity_s=sum(value);
            randomD=runif(1)#
            if(intensity_s/total_intensity>=randomD)
            {
                if(value[1,1]/total_intensity>=randomD)
                {             
                    sub_process[[1]]<-rbind(sub_process[[1]],c(s));
                    #print(sub_process[[1]])
                    t<-rbind(t,c(s));
                    n0<-1
                }
                else
                {
                    for(j in 2:M)
                    {
                        if ((sum(value[c(1:j-1),1])/total_intensity<randomD)&&(sum(value[c(1:j),1])/total_intensity>=randomD))
                        {
                            t<-rbind(t,c(s));
                            sub_process[[j]]<-rbind(sub_process[[j]],c(s));
                            n0<-j;
                        }
                    }
                }
                break; 
            }
            else
            {
                total_intensity<-sum(value);
                s<-s-log(runif(1))/total_intensity;
                if(s>T)
                {
                    out_range<-TRUE;
                    #print(out_range)
                    break;
                }
            }
        }
    }
    
    #Remove the entry t=0 in sub_process
    for(i in 1:M)
    {
        if (sub_process[[i]][1,1]==0)
            sub_process[[i]]=as.matrix(sub_process[[i]][-1,])
            #print(i)
    }
    return (list(t,sub_process));
}

MultiHawkesPlot<-function(history,mu,alpha,beta)
{
    process=history[[1]];
    sub_process=history[[2]];
    number=length(sub_process);
    time=seq(0,max(process),by=0.1);
    y<-array();
    for(i in 1:length(time))
    {
        y<-rbind(y,c(LambdaValue(time[i],process,sub_process,mu,alpha,beta)));
        print(time[i])
        print(y[i+1,])
    }
    #print(y)
    y=y[-1,]
    matplot(time,cbind(y[,1],y[,2]),type='l',xlab="Time",ylab="Intensity",lwd=2,cex.axis=1.5,cex.lab=1.5)
    legend("topleft",c("Dim 1","Dim 2"),lty=1:2,col=c(1,2))
}
T<-100
mu<-array(c(0.2,0.2),dim=c(2,1));
alpha<-array(c(0.2,0.7,0.4,0.4),dim=c(2,2));
beta<-array(c(1,1,1,1),dim=c(2,2));
x<-MultiHawkesSim(mu,alpha,beta,T);
#MultiHawkesPlot(x,mu,alpha,beta);