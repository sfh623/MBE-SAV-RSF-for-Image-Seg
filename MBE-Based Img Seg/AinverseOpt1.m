function [ res ] = AinverseOpt1(f,pa)
%AINVERSEOPT solves u + 2/3*dt (u''''-beta/eps^2 u'')  = f
%global N lambda 
[M,N]=size(f);
phi_0(1:M,1:M)=f;

phi_0(1:M,M:2*M-1)=fliplr(phi_0(1:M,1:M));
phi_0(M:2*M-1,1:M)=flipud(phi_0(1:M,1:M));
phi_0(M:2*M-1,M+1:2*M-1)=flipud(phi_0(1:M,M+1:2*M-1));
f = phi_0;

[M,N]=size(f);
rhs_f = fft2 (f);
index = [0:N/2  -N/2:-1]';
%index = [0:N/2-1 0 -N/2+1:-1]';
four = index.^4;
sec = index.^2;    % Don't forget negative sign
six = index.^6;
% temp = 1+ tau_n*lambda + eps*tau_n*((repmat(six,1,N)+repmat(six',N,1)+3*(repmat(sec,1,N).*repmat(four',N,1)+repmat(four,1,N).*repmat(sec',N,1))));%+...
            %beta/eps^2*(repmat(sec,1,N)+repmat(sec',N,1)));

% temp = 1 + pa.dt*pa.delta*(repmat(four,1,N)+2*repmat(sec,1,N).*repmat(sec',N,1)+repmat(four',N,1));
AA=repmat(four,1,N);
BB=repmat(sec,1,N);
CC=repmat(four',N,1);
DD=repmat(sec',N,1);
temp = 1 + pa.dt*pa.delta*(repmat(four,1,N)+2*repmat(sec,1,N).*repmat(sec',N,1)+repmat(four',N,1));

res1 = ifft2(rhs_f ./temp);
res=real(res1(1:100,1:100));
figure(23), imshow(res,[])

end