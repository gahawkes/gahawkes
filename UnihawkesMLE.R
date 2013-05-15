library(neldermead);
library(optimx)
library(Rsolnp)
library(alabama)
UniHawkesLikelihood<-function(theta,points)
{
    alpha<-theta[2];
    beta<-theta[3];
    lambda0<-theta[1];
    data_num=dim(points)[1]
    tn<-points[data_num,1]
    R<-matrix(0,data_num,1);
    for(i in 1:dim(R)[1])
    {
        if(i==1){R[i,1]<-0;}
        else
        {
            R[i,1]<-exp(-beta*(points[i,1]-points[i-1,1]))*(1+R[i-1,1])
        }
    }
    part1<-tn-lambda0*tn;
    part2_exp<-0;
    for(i in 1:data_num)
    {
        part2_exp<-part2_exp+exp(-beta*(tn-points[i,1]));
    }
    part2<-(alpha/beta)*(data_num-part2_exp);
    part3<-0;
    for(i in 1:data_num)
    {
        part3<-part3+log(lambda0+alpha*R[i,1]);
    }
    logle<-part1-part2+part3;
    return (logle);
}
MLEconstrain<-function(theta,points)
{
    alpha<-theta[2];
    beta<-theta[3]
    return (alpha/beta);
}
MLEhin<-function(theta,points)
{
    alpha<-theta[2];
    beta<-theta[3]
    h<-rep(NA,1);
    h[1]<-1-alpha/beta;
    h[2]<-alpha;
    h[3]<-beta;
    h;
}
#Load data
data<-read.table('goog_2010_dec.csv',header=TRUE,sep=',')
test<-as.POSIXct(as.character(data$DATE),format="%Y%m%d")
data$TIME<-as.POSIXct(paste(test,data$TIME))
data$DATE<-NULL
idx<-tapply(1:NROW(data),data$TIME,"[",1)
data_red<-data[idx,c(2:3)]  
day_01_idx<-(data_red$TIME>=as.POSIXct("2010-12-01 09:30:00",format="%Y-%m-%d %H:%M:%S"))&(data_red$TIME<=as.POSIXct("2010-12-01 10:00:00",format="%Y-%m-%d %H:%M:%S"))
day_01<-data_red[t(day_01_idx),]

#Separate positive and negative movement
pos_mov<-array(0,dim=c(1,1));
neg_mov<-array(0,dim=c(1,1));
for(i in 2:length(day_01[,1]))
{
    if(day_01[i,2]>day_01[i-1,2])
    {
        #print(c(as.double(day_01[i,1])))
        pos_mov<-rbind(pos_mov,as.double(day_01[i,1]));
        #print(pos_mov)
    }
    else if(day_01[i,2]<day_01[i-1,2])
    {
        neg_mov<-rbind(neg_mov,as.double(day_01[i,1]));
    }
    else
    {
        pos_mov<-rbind(pos_mov,as.double(day_01[i,1]));
        neg_mov<-rbind(neg_mov,as.double(day_01[i,1]));
    }
}
pos_mov<-as.matrix(pos_mov[-1,]);
pos_mov<-pos_mov-pos_mov[1,1];
pos_mov<-as.matrix(pos_mov[-1,]);
neg_mov<-as.matrix(neg_mov[-1,]);
neg_mov<-neg_mov-neg_mov[1,1];
neg_mov<-as.matrix(neg_mov[-1,]);
x<-UniHawkesSim(1,0.5,0.8,5000);
print('complete simulation')
print(x)
#data_points<-as.matrix(x);
#MLE 
print('start simulation')
#theta<-c(10,4,5);
theta<-c(1.2,0.6,0.8);
result<-UniHawkesLikelihood(theta,data_points);
#print(result)
theta<-c(1.2,0.6,0.8);
#result<-UniHawkesLikelihood(theta,data_points);
#print(result)

#par<-optimx(theta,UniHawkesLikelihood,points=data_points,method= "BFGS", control=list(trace=2,REPORT=1,maximize=TRUE))
#par<-nlm(UniHawkesLikelihood,c(0.4,0.6,0.5),points=data_points)
#par<-gosolnp(fun=UniHawkesLikelihood,LB=c(0,0,0),UB=c(2,2,2),ineqfun=MLEconstrain,ineqUB=c(1),n.sim=100,control=list(trace=1),points=x)

par<-auglag(par=theta,fn=UniHawkesLikelihood,hin=MLEhin,control.optim=list(trace=2,fnscale=-1),points=data_points)