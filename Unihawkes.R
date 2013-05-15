LambdaValue<-function(time,history,mu,alpha,beta)
{
    value<-mu;
    history_before_time<-matrix();
    history_before_time<-history[history[,1]<=time];
    for(i in 1:length(history_before_time))
    {
        if(time<history[1,1])
        {
            value<-mu;
        }
        else
        {
            value<-value+alpha*exp(-beta*(time-history_before_time[i]));
        }
    }
    return (value);
}

UniHawkesSim<-function(mu,alpha,beta,T)
{
    lambda<-mu;
    n<-1;s<-(-1/lambda*log(runif(1)));
    #CC print(s)
    t<-matrix(0,1,1);
    if (s<=T)
    {
        t[1,1]<-s;
    }
    else
    {
        return (t); 
    }
    #print(t)
    #CC
    
    out_range<-FALSE;
    while(out_range==FALSE)
    {
        n<-(n+1);
        lambda<-LambdaValue(t[n-1,1],t,mu,alpha,beta)+alpha;
        s<-s-1/lambda*log(runif(1));#
        #CC print(s)
        #print(s)
        if (s>T){break;}
        
        while(out_range==FALSE)
        { 
            #Rejection test
            lambda_s<-LambdaValue(s,t,mu,alpha,beta);
            #print(lambda_s)
            if(lambda_s/lambda>=runif(1))#r
            {
                t<-rbind(t,c(s));
                break; 
            }
            else
            {
                lambda<-lambda_s;
                s<-s-1/lambda*log(runif(1));
                if(s>=T)
                {
                    out_range<-TRUE;
                    break;
                }
            }
        }
    }
    return (t);
}

UniHawkesPlot<-function(history,mu,alpha,beta)
{
    time<-seq(0,max(history[,1]),by=0.1);
    y<-c();
    lambda_max<-c();
    for(i in 1:length(time))
    {
        y[i]<-LambdaValue(time[i],history,mu,alpha,beta);
    }
    plot(time,y,xlab="Time",ylab="Intensity",type='l',lwd=2,cex.axis=1.5,cex.lab=1.5);
}
set.seed(2)
x1=UniHawkesSim(1.2,0.6,0.8,50);
UniHawkesPlot(x1,1.2,0.6,0.8);