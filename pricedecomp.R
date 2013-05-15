PriceDecomposition<-function(price)
{
    n=dim(price)[1]
    #x will be transmitted to mle function to estimate parameters
    x=list();
    #Store the total event time into x[[1]]
    x[[1]]=as.matrix(price[-1,1])
    #print(x[[1]])
    
    #sub_process:used to store the time of N+ and N-
    sub_process=list();
    sub_process[[1]]=array(0,dim=c(n,1));#positive coordinate
    sub_process[[2]]=array(0,dim=c(n,1));#negative coordinate
    n_positive=0;
    n_negative=0;
    
    #Decompose the price into N+ and N-
    for(i in 2:n)
    {
        if (price[i,2]>price[i-1,2])
        {
            n_positive=n_positive+1;
            sub_process[[1]][n_positive,1]=price[i,1]
        }
        else if(price[i,2]==price[i-1,2])
        {
            n_positive=n_positive+1
            n_negative=n_negative+1
            sub_process[[1]][n_positive,1]=price[i,1]
            sub_process[[2]][n_negative,1]=price[i,1]
        }
        else
        {
            n_negative=n_negative+1
            sub_process[[2]][n_negative,1]=price[i,1]
        }
    }
    
    #print(sub_process[[1]])
    #Delete the unfilled entries in sub_process
    sub_process[[1]]=as.matrix(sub_process[[1]][-(n_positive+1):-n,1]);
    sub_process[[2]]=as.matrix(sub_process[[2]][-(n_negative+1):-n,1]);
    #print(sub_process[[1]])
    #delete the first row of x[[1]]
    x[[1]]=as.matrix(x[[1]][-1,]);
    x[[2]]=sub_process;
    
    return (x)
}