
R <- 0.15					#Radious
D <- 0.1					#Diameter of lens
n1 <- 1.0					#Refraction index
n2 <- 1.5					#Refraction index
R <- D*0.51
#Gaussian function (Radious,Hight,incoming refraction index, material refraction index, Use approximation BOOLEAN)
Gaussian <- function(r,h,n1,n2,b){
	a1 <- 0					#Paraxial Approximation
	if(!b){
		a1 <- alph1(h,r)		#None Approximated Angle
	}
	a2 <-alph2(r,a1,n1,n2)
	f <- r*sin(a2)/cos(a2)+r
	return(f)
}

#Refraction angle to norm of surface
alph2 <- function(r,a1,n1,n2){
	a2 <- asin(sin(a1)*(n1/n2))
	return(a2)
}
#Light angle without Paraxial Approxation to norm of surface
alph1 <- function(h,r){
	a1 <- asin(h/r)
	return(a1)
}

#Refreaction index calculation of Glass material BK7
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
	return(n)
}
par(mfrow = c(2,2));

#Paraxoide Approximation applied
Gauss_Approx <- function(x) Gaussian(R,x,n1,n2,TRUE);
fa <- Vectorize(Gauss_Approx);
plot.function(fa, from=0, to=0.1, xlab="Hight", ylab="Focus Point",
		ylim=c(fa(0)-0.01,0.25));
#No Paraxoide Approximation
Gauss <- function(x) Gaussian(R,x,n1,n2,FALSE);
f <- Vectorize(Gauss);
plot.function(f,from=0,to=0.1, add=TRUE, col="red",
		);
n2v <- Vectorize(BK7n);
plot.function(n2v, from=(400/10^9), to=(700/10^9), ylab="Reflection Index", xlab="Wavelength");

#BK7 replace material of lens	
h <- 0.025;
f_chrom <- function(x) Gaussian(R,h,n1,BK7n(x),FALSE);
v_chrom <- Vectorize(f_chrom);

plot.function(v_chrom, from=(400/10^9), to=(700/10^9), xlab="", ylab="");

