BK7n <- function(a){
	a1 <- 2.271176
	a2 <- -9.700709/10^9
	a3 <- 0.0110971/10^6
	a4 <- 4.622809/10^11
	a5 <- 1.616105/10^11
	a6 <- -8.285043/10^13
	n4 <- (a1+a2*a^2+a3/a^2+a4/a^4+a5*1/a^6+a6/a^8)^2
	n2 <- sqrt(n4)
	n <- sqrt(n2)
	#print(n)
	return(n)
}

R <- 0.15
D <- 0.1
n1 <- 1.0
n2 <- 1.5

Gaussian <- function(r,h,n1,n2,b){
	a1 <- 0					#Paraxial Approximation
	if(!b){
		a1 <- alph1(h,r)		#None Approximated Angle
	}
	a2 <-alph2(r,a1,n1,n2)
	f <- r*sin(a2)/cos(a2)+r
	return(f)
}

#Incomming light angle in relation of the norm.
alph2 <- function(r,a1,n1,n2){
	a2 <- asin(sin(a1)*(n1/n2))
	return(a2)
}
alph1 <- function(h,r){
	a1 <- asin(h/r)
	return(a1)
}

Gauss_Approx <- function(x) Gaussian(R,x,n1,n2,TRUE);
fa <- Vectorize(Gauss_Approx);
curve(fa,from=0,to=0.1, xlab="h", ylab="f");

Gauss <- function(x) Gaussian(R,x,n1,n2,FALSE);
f <- Vectorize(Gauss);
curve(f,from=0,to=0.1, xlab="h", ylab="f", add=TRUE);


print("Done!")

"
Gaussian(R,h,n1,n2,TRUE)
Gaussian(R,h,n1,n2,FALSE)

BK7n(400/10^9)
BK7n(500/10^9)
BK7n(600/10^9)
BK7n(700/10^9)
"
