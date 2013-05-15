PriceSimulation<-function(x)
{
    #Input:
    #   x - The result produced from MultiHawkesSim
    #Output:
    #   Price - The time and price in a nX2 matrix
    #   The initial price are set to be 100
    scale=0.01
    rise<-x[[2]][[1]]
    fall<-x[[2]][[2]]
    rise_dim<-length(rise)
    fall_dim<-length(fall)
    price<-array(0,dim=c(rise_dim+fall_dim,2))
    #print(price)
    i<-1;j<-1;
    for(n in 1:(rise_dim+fall_dim))
    {
        if(rise[i]>fall[j])
        {
            price[n,1]<-fall[j];
            if(n==1){price[n,2]<--scale;}
            else{price[n,2]<-price[n-1,2]-scale;}
            j<-min(j+1,fall_dim);
        }
        else if(rise[i]==fall[j])
        {
            price[n,1]<-rise[i];
            if(n==1){price[n,2]<-0;}
            else{price[n,2]<-price[n-1,2];}
            j<-min(j+1,fall_dim);
            i<-min(i+1,rise_dim);
        }
        else
        {
            price[n,1]<-rise[i];
            if(n==1){price[n,2]<-scale;}
            else{price[n,2]<-price[n-1,2]+scale;}
            i<-min(i+1,rise_dim);
        }
    }
    price[,2]=price[,2]+100;
    return (price);
}
