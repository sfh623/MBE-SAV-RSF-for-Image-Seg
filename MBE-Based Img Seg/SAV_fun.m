function phi_0 = SAV_fun(phi_0,Img,pa)
%FFT prepare
[M, N] = size(phi_0);
ii = 1:1:M;
ii = ii';
jj = 1:1:N;
Z1 = 2*pi*(ii-1)/M;
Z2 = 2*pi*(jj-1)/N;
Z1_c = cos(Z1);
Z2_c = cos(Z2);
Z1_2 = repmat(Z1_c, 1, N);
Z2_2 = repmat(Z2_c, M, 1);
Lap2_FFT = (Z1_2 + Z2_2 - 2).^2;
% Lap_FFT = (Z1_2 + Z2_2 - 2)./2;
A = 1+pa.dt*pa.delta*pa.mu*Lap2_FFT ;

eng_E1 = fun_E1(phi_0,Img,pa);
r = sqrt(eng_E1);
b = fun_U(phi_0,Img,pa)/r;

Gb= -b;
vb=b(:);
vphi_0=phi_0(:);
c = phi_0 + pa.dt*r*Gb - (pa.dt/2)*(vb'*vphi_0)*Gb;
fft_c=fftn(c);
res1=real(ifftn(fft_c./A));

fft_Gb=fftn(Gb);
res2=real(ifftn(fft_Gb./A));

%res1 = AinverseOpt1(c,pa);
% res2 = AinverseOpt1(Gb,pa);
gamma = -sum(sum((b.*res2)));

d = (sum(sum((b.*res1))))/(1+(pa.dt*gamma/2));
e = c + (pa.dt/2)*d*Gb;
fft_e = fftn(e);
Phi = fft_e./A;
phi = real(ifftn(Phi));
%phi = AinverseOpt1(e,pa);

phi_0=phi;
end


