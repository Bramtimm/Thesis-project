y <- function(x) 2/sqrt(2*pi*(x-1)^3)*exp(-(2-1*(x-1))^2/(2*(x-1)))
lik.fun <- function(x) (1/2)*log(2^2/(2*pi))-(3/2)*log(x-1)-2^2*((x-1-1)^2/(2*(x-1))) 

x<-seq(1,100,.01)
a<-seq(0,100,0.01)
plot(x,y(x),type='l',ylim=c(0,2),xlim=c(0,10))
lines(x,exp(lik.fun(x)),col="red")
lines(density(y(x)))
y
density(y(x))
a<-y(x)
plot(x,a)
y
