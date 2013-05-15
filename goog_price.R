library(alabama)

#Read in data and extract the time range
data<-read.table('goog_2010_dec.csv',header=TRUE,sep=',')
test<-as.POSIXct(as.character(data$DATE),format="%Y%m%d")
data$TIME<-as.POSIXct(paste(test,data$TIME))
data$DATE<-NULL
idx<-tapply(1:NROW(data),data$TIME,"[",1)
data_red<-data[idx,c(2:3)]  
day_01_idx<-(data_red$TIME>=as.POSIXct("2010-12-01 10:45:00",format="%Y-%m-%d %H:%M:%S"))&(data_red$TIME<=as.POSIXct("2010-12-01 10:55:00",format="%Y-%m-%d %H:%M:%S"))
day_01<-data_red[t(day_01_idx),]
#plot(day_01[,2]~day_01[,1],type='s',xlab='time',ylab='Price')

#transform time type to double
n=dim(day_01)[1]
goog_data=array(0,dim=c(n,2))
goog_data[,1]=as.double(day_01[,1])
goog_data[,2]=as.double(day_01[,2])
goog_data[,1]=goog_data[,1]-goog_data[1,1]

# #Divide the price into training set and testing set
# #Traing set: 300 seonds
# #Testing set: 300 seconds
goog_trainingset=goog_data[goog_data[,1]<=300,]
goog_testingset=goog_data[goog_data[,1]>300]

#divide the price into positive movement and negativement
goog_trainingset_decomposed=PriceDecomposition(goog_trainingset)

#Learn the parameters of the Hawkes process
#Set initial estimation of parameters
goog_mu_0<-array(c(0.2,0.2),dim=c(2,1));
goog_alpha_0<-array(c(0.1,0.1,0.1,0.1),dim=c(2,2));
goog_beta_0<-array(c(0.5,0.5),dim=c(2,1));
goog_theta<-c(goog_mu_0,goog_alpha_0,goog_beta_0);
goog_par<-constrOptim.nl(par=goog_theta,fn=BiHawkesLikelihood,hin=MLEhin,control.outer=list(trace=FALSE),
                    control.optim=list(trace=FALSE,fnscale=-0.01),points=goog_trainingset_decomposed);
#print(par)
goog_mu_estimate=array(goog_par$par[1:2],dim=c(2,1))
goog_alpha_estimate=array(goog_par$par[3:6],dim=c(2,2))
goog_beta_estimate=array(c(goog_par$par[7],goog_par$par[8],goog_par$par[7],goog_par$par[8]),dim=c(2,2))


#Use the estimated parameters to do monte carlo simulation
print('Begin Simulation')
T=300;
loop=500;
goog_samplepaths=matrix(nrow=1000,ncol=loop)
for(i in 1:loop)
{
    goog_sim1<-MultiHawkesSim(goog_mu_estimate,goog_alpha_estimate,goog_beta_estimate,T);
    goog_price1=PriceSimulation(goog_sim1);
    n=dim(goog_price1)[1]
    for(j in 1:n)
    {
        goog_samplepaths[j,i]=goog_price1[j,2]
    }
    
    if(i==1)
    {
        plot(goog_price1[,2]~goog_price1[,1],type='s',xlab='seconds',ylab='Price',xlim=c(0,305),ylim=c(99,101));
        par(new=TRUE)
    }
    plot(goog_price1[,2]~goog_price1[,1],type='s',axes=FALSE,xlab='',ylab='',xlim=c(0,305),ylim=c(99,101));
    par(new=TRUE)
}

