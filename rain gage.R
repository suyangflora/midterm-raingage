
library(stringr)
theFiles<-dir("C:/Users/Yangsu/Desktop/rain data/",pattern="\\.txt")
theFiles
for (a in theFiles){
  nameToUse<-str_sub(string=a,start=1,end=7)
  temp<-read.csv(file=file.path("C:/Users/Yangsu/Desktop/rain data",a),skip=2,stringsAsFactors = F)
  assign(x=nameToUse,value=temp)
}

rain<-rbind(L0001.t,	L0002.t,	L0003.t,	L0004.t,	L0005.t,	L0006.t,	L0007.t,	L0008.t,	L0009.t,	
          L0010.t,	L0011.t,	L0012.t,  L0101.t,	L0102.t,	L0103.t,	L0104.t,	L0105.t,	L0106.t,	
          L0107.t,	L0108.t,	L0109.t,	L0110.t,	L0111.t,	L0112.t,  L0201.t,	L0202.t,	L0203.t,	
          L0204.t,	L0205.t,	L0206.t,	L0207.t,	L0208.t,	L0209.t,	L0210.t,	L0211.t,	L0212.t,
          L0301.t,	L0302.t,	L0303.t,	L0304.t,	L0305.t,	L0306.t,	L0307.t,	L0308.t,	L0309.t,	
          L0310.t,	L0311.t,	L0312.t,  L0401.t,	L0402.t,	L0403.t,	L0404.t,	L0405.t,	L0406.t,
          L0407.t,	L0408.t,	L0409.t,	L0410.t,	L0411.t,	L0412.t)

dim(rain)
colnames(rain) <- 0:24
head(rain)

## to calculate each rain gage data during a rain storm, we add up numbers between two "----"
## "----" means no raining
## "T   " means a little rain but uncountable, so if T is between two numbers, then it is included in one storm
## I set "T   " to be the minor number so that it helps me with later calculation
## for "M   " I do not know what it means but it is always between two "----" so I set it to be 0

rain[rain=="----"] <- 0
rain[rain=="M   "] <- 0
rain[rain=="M"] <- 0
rain[rain=="T   "] <- 10^(-8)
head(rain)

## delete the first column
r01<-rain[,(2:25)]
head(r01)

## change class of character to class of numeric
bos <- as.data.frame(sapply(r01, as.numeric))
bosrain<-bos[complete.cases(bos), ]
View(bosrain)

## change data frame into vector
brain <- as.vector(t(bosrain))

## use function to get rain data for each storm
## we also build up a vector to put the results in

sum <- 0
j=1
vector<-0
for(i in 1:length(brain)) 
{
  if(brain[i] != 0)
  {
    sum=sum+brain[i]
  }
  if(brain[i]==0 && sum!=0)
  {
    vector[j]=sum
    j=j+1
    sum=0
  }
  if(brain[i]!=0 & i==length(brain))
  {
    vector[j]=sum 
  }
}

## in order to delete those T without surrounding by numbers, we choose to keep only two digits parts
vector1<-round(vector, 2)
v2<-vector1[vector1 != 0.00]

## so that v2 is boston logan airport rain gage data which looks the same as illinois rain data
## next we start with rain gage distribution analysis
class(v2)
logan <- data.frame(v2)
colnames(logan) <- "x"
library(ggplot2)
qplot(x, data=logan, geom = "histogram",binwidth=.15)

## it looks like gamma distribution

mean(logan$x)
var(logan$x)


alpha <- mean(logan$x)^2/var(logan$x)  # alpha = 0.36
lambda <- mean(logan$x)/var(logan$x)   # lambda = 1.28


gam<-(lambda^(alpha)/gamma(alpha))*(logan$x^(alpha-1))*exp(-lambda*logan$x)
gam1<-data.frame(gam)


## gamma distribution density plot
library(ggplot2)
ggplot(gam1, aes(x=gam1$gam)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot

## data distribution density plot
ggplot(logan, aes(x=logan$x)) + 
  geom_histogram(aes(y=..density..),      # Histogram with density instead of count on y-axis
                 binwidth=.5,
                 colour="black", fill="white") +
  geom_density(alpha=.2, fill="#FF6666")  # Overlay with transparent density plot


# for Variance & confidence interval

## using MEM
lam<-mean(logan$x)/(sd(logan$x)^2)
alp<-(mean(logan$x))^2/(sd(logan$x)^2)

B<-1000
Tboot1<-rep(0,B)
Tboot2<-rep(0,B)
for(i in 1:B){
  x <- sample(logan$x,1000,replace=TRUE)
  Tboot1[i] <- mean(x)/(sd(x)^2)
  Tboot2[i] <- (mean(x))^2/(sd(x)^2)
}

Percentile1 <- c(quantile(Tboot1,.025),quantile(Tboot1,.975))
pivotal1 <- c((2*lam - quantile(Tboot1, .975)),(2*lam - quantile(Tboot1, .025))) 

cat("Method       95% Interval\n")
cat("Pivotal1     (", pivotal1[1], ",     ", pivotal1[2], ") \n")
cat("Percentile1  (", Percentile1[1], ",    ", Percentile1[2], ") \n")


Percentile2 <- c(quantile(Tboot2,.025),quantile(Tboot2,.975))
pivotal2 <- c((2*alp - quantile(Tboot2, .975)),(2*alp - quantile(Tboot2, .025))) 

cat("Method       95% Interval\n")
cat("Pivotal2     (", pivotal2[1], ",     ", pivotal2[2], ") \n")
cat("Percentile2  (", Percentile2[1], ",    ", Percentile2[2], ") \n")



##  for MLE method

  mle.x <- logan$x
  n <- length(logan$x)
  # first we need to have alpha and lambda from MEM
  mem.alp <- mean(mle.x)^2/var(mle.x)
  mem.lam <- (mean(mle.x))/var(mle.x)
  mem.alp
  mem.lam
  
  # second we use MLE to get parameter value
  minus.likelihood <- function(theta) {-(n*theta[1]*log(theta[2])-n*lgamma(theta[1])+(theta[1]-1)*sum(log(mle.x))-theta[2]*sum(mle.x))}
  
  max.likelihood <- nlminb(start=c(mem.alp, mem.lam), obj = minus.likelihood) 
  
  max.likelihood$par


 