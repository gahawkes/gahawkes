goog_simresult=array(0,dim=c(500,1))
for(i in 1:500)
{
    for(j in 1:1000)
    {
        if(is.na(goog_samplepaths[j,i]))
        {
            goog_simresult[i,1]=goog_samplepaths[j-1,i]
            break
        }
    }
}