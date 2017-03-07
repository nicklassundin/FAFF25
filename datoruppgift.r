
R <- 0.15					#Radious
D <- 0.1					#Diameter of lens
n1 <- 1.0					#Refraction index
n2 <- 1.5					#Refraction index
#Gaussian function (Radious,Hight,incoming refraction index, material refraction index, Use approximation BOOLEAN)
Gaussian <- function(r,h,n1,n2,b){
	a1 = 0					#Paraxial Approximation
	if(!b){
		a1 = alph1(h,r)		#None Approximated Angle
	}
	a2 = alph2(r,a1,n1,n2)
	f = r*sin(a2)/cos(a2)+r
	return(f)
}

#Refraction angle to norm of surface
alph2 <- function(r,a1,n1,n2){
	a2 = asin(sin(a1)*(n1/n2))
	return(a2)
}
#Light angle without Paraxial Approxation to norm of surface
alph1 <- function(h,r){
	a1 = asin(h/r)
	return(a1)
}

#Refreaction index calculation of Glass material BK7
BK7n <- function(x){
	a1 = 2.271176
	a2 = -9.700709*(10^-3)*(10^-6)^-2
	a3 = 0.0110971*(10^-6)*(10^-6)^2
	a4 = 4.622809*(10^-5)*(10^-6)^4
	a5 = 1.616105*(10^-5)*(10^-6)^6
	a6 = -8.285043*(10^-7)*(10^-6)^8
	n2 = a1+a2*x^2+a3*x^(-2)+a4*x^(-4)+a5*x^(-6)+a6*x^(-8)
	n = a1;
	n = abs(sqrt(as.complex(n2)));
	return(n);
}
par(mfrow = c(2,3));

#Paraxoide Approximation applied
Gauss_Approx <- function(x) Gaussian(R,x,n1,n2,TRUE);
fa <- Vectorize(Gauss_Approx);
plot.function(fa, from=0, to=D/2, xlab="Hight", ylab="Focus Point",
		ylim=c(fa(0)-0.01,0.20));

#No Paraxoide Approximation
Gauss <- function(x) Gaussian(R,x,n1,n2,FALSE);
f <- Vectorize(Gauss);
plot.function(f,from=0,to=D/2, add=TRUE, col="red");

#BK7n Reflection index
n2v <- Vectorize(BK7n);
plot.function(n2v, from=(400/(10^9)), to=(700/(10^9)), ylab="Reflection Index", xlab="Wavelength");

#BK7 replace material of lens	
h <- 0.025;
f_chrom <- function(x){
	f = Gaussian(R, h, n1, BK7n(x), FALSE);
	return(f)
}
v_chrom <- Vectorize(f_chrom);
plot.function(v_chrom, from=(400/(10^9)), to=(700/(10^9)), xlab="", ylab="");

#Assignment no. 2
L <- 0.2			#Length
D <- 0.008 			#Diameter
tb <- 200/10^6			#pulse duration
tau <- 230/10^6			#Lifespan
N0 <- 1.4*10^20			#Number of Ions cm^-3
sigma <<- 2.8/10^23 		#
c <- 299792458			#Speed of Light m/s

V <- L*pi*(D/2)^2;		#cavity Volyme
B <- sigma*c/V;			#Probability of stimulated emission ion and photon

N_inf <<- 0.01*N0;
P <- N_inf/tau;			#Pump strength

#Assignment 2:b definitions
R1 <- 1;
R2 <- 0.05;
tb <- 200/10^6;

tau_c <- function(r1,r2) {	#Lifespan in cavity for photons
	tau_r = -2*L/(c*(log(r1)+log(r2)));
	return(tau_r)
}


#Differential eqvations:
N_prim <- function(N, Phi){		#Number of Ions
	y = P-B*N*Phi-N/tau;
	return(y)
}

Phi_prim <- function(Phi, N) {	#
	y = B*V*N*(Phi+1)-Phi/tau_c(R1,R2);
	return(y)
}


Solv <- function(f0, f_prim, g_prim, t){
	h = (t[1]-t[2]);
	f = rep(0, length(t));
	g = rep(0, length(t));
	f[1] <- (f0[1]);
	g[1] <- (f0[2]);
	for(i in as.single(1:length(t))){
		f[i+1] = f[i] - f_prim(f[i],g[i])*h;
		g[i+1] = g[i] - g_prim(g[i],f[i])*h;
	}
	return(f);
}

N <- function(x){
	return(Solv(c(N0,0), N_prim, Phi_prim, x))
}
Phi <- function(x){
	return(Solv(c(0,N0), Phi_prim, N_prim, x))
}
x0 = N(seq(0, 0.0002, length=20));
x1 = Phi(seq(0, 0.0002, length=20));
plot(x0);
plot(x1);
dev.off()

#Write to file #1
pdf("para_approx.pdf", width=7, height=5)
plot.function(f,from=0,to=D/2, col="red");
plot.function(fa, from=0, to=D/2, xlab="Hight",add=TRUE, ylab="Focus Point");
dev.off()

print(2)
#Write to file #2
pdf("BK7_index.pdf", width=7, height=5)
plot.function(n2v, from=(400/(10^9)), to=(700/(10^9)), ylab="Reflection Index", xlab="Wavelength");
dev.off()

#Write to file #3
pdf("BK7_abo.pdf", width=7, height=5)
plot.function(v_chrom, from=(400/(10^9)), to=(700/(10^9)), xlab="", ylab="");
dev.off()

#Write to file #4
pdf("N.pdf", width=7, height=5)
plot(x0)
dev.off()

