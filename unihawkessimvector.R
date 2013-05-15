LambdaValue<-function(time,history,mu,alpha,beta)
{
    value<-mu;
    history_before_time<-c();
    history_before_time<-history[history<=time];
    for(i in 1:length(history_before_time))
    {
        if(time<history[1])
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

UniHawkesSimVector<-function(mu,alpha,beta,T)
{
    lambda<-mu;
    n<-1;s<-(-1/lambda*log(runif(1)));
    t<-c();
    if (s<=T)
    {
        t<-c(t,s);
    }
    else
    {
        return (t); 
    }
    out_range<-FALSE;
    while(out_range==FALSE)
    {
        n<-(n+1);
        lambda<-LambdaValue(t[n-1],t,mu,alpha,beta)+alpha;
        s<-s-1/lambda*log(runif(1));
        if (s>T)
        {
            break;
        }
        while(out_range==FALSE)
        { 
            #Rejection test
            lambda_s<-LambdaValue(s,t,mu,alpha,beta);
            if(lambda_s/lambda>=runif(1))
            {
                t<-c(t,s);
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
    time<-seq(0,max(history),by=0.1);
    y<-c();
    lambda_max<-c();
    for(i in 1:length(time))
    {
        y[i]<-LambdaValue(time[i],history,mu,alpha,beta);
    }
    plot(time,y,type='l');
}
set.seed(2)
print('vector')
start<-proc.time()
x<-UniHawkesSimVector(1.2,0.6,0.8,1000);
dur1<-proc.time()-start
print(dur1)
#UniHawkesPlot(x,1.2,0.6,0.8);