library(alabama)
#Code for monte carlo simulation using sythetic data

#Simulate an upward trend
#Set parameters
mu<-array(c(0.12,0.12),dim=c(2,1));
alpha<-array(c(0.02,0.07,0.07,0),dim=c(2,2));
beta<-array(c(0.11,0.11,0.11,0.11),dim=c(2,2));
T<-2000;
scale<-0.01;

#Simulate the 500 seonds sythetic price movement y
set.seed(2)
x<-MultiHawkesSim(mu,alpha,beta,T);
price=PriceSimulation(x);

#Delete repeated time
n=dim(price)[1]
price=price[-(n-2):-n,]

plot(price[,2]~price[,1],type='s',xlab='seconds',ylab='Price');

#Divide the sythetic price into training set and testing set
#Traing set: 300 seonds
#Testing set: 300 seconds
trainingset=price[price[,1]<=300,]
testingset=price[price[,1]>300]

#divide the price into positive movement and negativement
trainingset_decomposed=PriceDecomposition(trainingset)

#Learn the parameters of the Hawkes process
#Set initial estimation of parameters
mu_0<-array(c(0.5,0.5),dim=c(2,1));
alpha_0<-array(c(0.2,0.5,0.2,0.5),dim=c(2,2));
beta_0<-array(c(1,1),dim=c(2,1));
theta<-c(mu_0,alpha_0,beta_0);
par<-constrOptim.nl(par=theta,fn=BiHawkesLikelihood,hin=MLEhin,control.outer=list(trace=FALSE),
                    control.optim=list(trace=FALSE,fnscale=-0.01),points=trainingset_decomposed);
#print(par)
mu_estimate=array(par$par[1:2],dim=c(2,1))
alpha_estimate=array(par$par[3:6],dim=c(2,2))
beta_estimate=array(c(par$par[7],par$par[8],par$par[7],par$par[8]),dim=c(2,2))
#print(beta_estimate)

#Use the estimated parameters to do monte carlo simulation
T=300;
loop=500;
samplepaths=matrix(nrow=1000,ncol=loop)
for(i in 1:loop)
{
    sim1<-MultiHawkesSim(mu_estimate,alpha_estimate,beta_estimate,T);
    price1=PriceSimulation(sim1);
    n=dim(price1)[1]
    for(j in 1:n)
    {
        samplepaths[j,i]=price1[j,2]
    }
 
    if(i==1)
    {
        plot(price1[,2]~price1[,1],type='s',xlab='seconds',ylab='Price',xlim=c(0,305),ylim=c(99,101));
        par(new=TRUE)
    }
    plot(price1[,2]~price1[,1],type='s',axes=FALSE,xlab='',ylab='',xlim=c(0,305),ylim=c(99,101));
    par(new=TRUE)
}
