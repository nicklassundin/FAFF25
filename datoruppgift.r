
R <- 0.15					#Radious
D <- 0.1					#Diameter of lens
n1 <- 1.0					#Refraction index
n2 <- 1.5					#Refraction index
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
	if(is.nan(a1)){
		return
	}
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
plot.function(fa, from=0, to=D/2, xlab="Hight", ylab="Focus Point",
		ylim=c(fa(0)-0.01,0.25));
#No Paraxoide Approximation
Gauss <- function(x) Gaussian(R,x,n1,n2,FALSE);
f <- Vectorize(Gauss);
plot.function(f,from=0,to=D/2, add=TRUE, col="red",
		);
n2v <- Vectorize(BK7n);
plot.function(n2v, from=(400/10^9), to=(700/10^9), ylab="Reflection Index", xlab="Wavelength");

#BK7 replace material of lens	
h <- 0.025;
f_chrom <- function(a){
	bk7n2 <- BK7n(a);
	f <- Gaussian(R, h, n1, bk7n2, FALSE);
	return(f)
}
v_chrom <- Vectorize(f_chrom);
plot.function(v_chrom, from=(400/10^9), to=(700/10^9), xlab="", ylab="");

#Assignment no. 2
L <- 0.2			#Length
D <- 0.008 			#Diameter
t <- 200/10^6			#pulse duration
tau <- 230/10^6			#Lifespan
N0 <- 1.4*10^20			#Number of Ions cm^-3
sigma <- 2.8/10^23 		#
c <- 299792458			#Speed of Light m/s

V <- L*pi*(D/2)^2;		#cavity Volyme
B <- sigma*c/V;			#Probability of stimulated emission ion and photon

N_inf <- 0.01*N0;
P <- N_inf/tau;			#Pump strength

#Assignment 2:b definitions
R1 <- 1;
R2 <- 0.05;
T <- 200/10^6;

tau_c <- function(r1,r2) {	#Lifespan in cavity for photons
	tau_r <- -2*L/(c*(log(r1)+log(r2)));
	return(tau_r)
}

#Differential eqvations:
N_prim <- function(Phi,N){		#Number of Ions
	y <- P-B*N*Phi-N/tau;
	return(y)
}

Phi_prim <- function(Phi, N) {	#
	y <- B*V*N*(Phi+1)-Phi/tau_c(R1,R2);
	if(is.infinite(y)){
		y <- 0;
	}
	return(y)
}

#Differential Solver
Num_Solv_Diff <- function(Tn, s){
	m <- matrix(3, 1:(2*s));
	i <- 1;
	h <- Tn/s;
	N <- N0+0;
	Phi <- 0;
	m[0,0] <- 0;
	m[0,1] <- N;
	m[0,2] <- Phi;
	while(i<s){
		N <- N + N_prim(Phi, N);
		Phi <- Phi + Phi_prim(Phi, N);
		m[i,0] <- i*h;
		m[i,1] <- N;
		m[i,2] <- Phi;
		i <- i+1;
	}
	n_v <- as.vector(m[0,], mode="any");
	phi_v <- as.vector(m[1,], mode="any");	
	return(c(n_v, phi_v));
}

m <- Num_Solv_Diff(T,30)
print(m)
plot.function(m[0]);
