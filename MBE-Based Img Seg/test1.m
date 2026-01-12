
clear all;clc


Img = 2*diag(ones(10,1))+diag(ones(9,1),-1)+diag(ones(9,1),+1);
u0 = 4*diag(ones(10,1))+diag(ones(9,1),-1)+diag(ones(9,1),+1);

N = 64;
Img = rand(N);
u0 = rand(N);




timestep = .001;
iterNum = 1;

lambda1 = 1;
lambda2 = 1;
% nu = 0.01*255*255;
nu = 1;
% mu = 1;
mu = 1;

epsilon = 1.0;
sigma = 3.0;   

K=fspecial('gaussian',round(2*sigma)*2+1,sigma);   
I = Img;
KI=conv2(Img,K,'same'); 

KONE=conv2(ones(size(Img)),K,'same'); 

%%%
delta = 0;


u = u0;
for n=1:iterNum
    u = Reg_RSF_MBE_Fourier(u,I,K,KI,KONE,nu,delta,timestep,1,mu,lambda1,lambda2,epsilon,n);
end

phi_0 = u0;
pa.nu = nu;
pa.mu = mu;
pa.sigma = sigma;
pa.I = I;
pa.Img = Img;
pa.K = K;
pa.KI = KI;
pa.KONE = KONE;
pa.delta = delta;
pa.epsilon = epsilon;
pa.timestep = timestep;
pa.iterNum = iterNum;
pa.lambda1 = lambda1;
pa.lambda2 = lambda2;
pa.dt = timestep;

for n=1:pa.iterNum
    phi_0 = SAV_fun(phi_0,Img,pa);
end


x = linspace(0,1,N^2);
plot(x,phi_0(:),'r-',x,u(:),'b--')
norm(u - phi_0,inf)
