% Fourier for the RSF with MBE.
function phi = Reg_RSF_MBE(phi,Img,Ksigma,KI,KONE,nu,delta,dt,iterN,mu,lambda1,lambda2,epsilon)
% (u0,Img,Ksigma,KI,KONE,nu,timestep,mu,lambda1,lambda2,epsilon,numIter)
% %% FFT prepare
% [M, N] = size(phi);
% ii = 1:1:M;
% ii = ii';
% jj = 1:1:N;
% Z1 = 2*pi*(ii-1)/M;
% Z2 = 2*pi*(jj-1)/N;
% Z1_c = cos(Z1);
% Z2_c = cos(Z2);
% Z1_2 = repmat(Z1_c, 1, N);
% Z2_2 = repmat(Z2_c, M, 1);
% Lap2_FFT = 4 * (Z1_2 + Z2_2 - 2).^2;
% Lap_FFT = 2 * (Z1_2 + Z2_2 - 2);


%% Iteration
for i=1:iterN
    
    
    % r
    phi_lap = phi([end 1:end-1],:) + phi([2:end 1],:) + phi(:,[end 1:end-1]) ...
        + phi(:,[2:end 1]) - 4*phi;
    phi_lap2 = phi_lap([end 1:end-1],:) + phi_lap([2:end 1],:) ...
        + phi_lap(:,[end 1:end-1]) + phi_lap(:,[2:end 1]) - 4*phi_lap;
    
    phi_x = (phi([end 1:end-1],:) - phi([2:end 1],:))/2;
    phi_y = (phi(:,[end 1:end-1]) - phi(:,[2:end 1]))/2;
    phi_g2 = phi_x.^2 + phi_y.^2;
    d_x = phi_g2.*phi_x;
    d_y = phi_g2.*phi_y;
    dd_x = (d_x([end 1:end-1],:) - d_x([2:end 1],:))/2;
    dd_y = (d_y(:,[end 1:end-1]) - d_y(:,[2:end 1]))/2;
    A_phi = dd_x + dd_y;
    
    T_mbe = -delta * phi_lap2 + A_phi - phi_lap;
    
    %% RSF
    K=curvature_central(phi);
    DrcU=(epsilon/pi)./(epsilon^2.+phi.^2);               % eq.(9)
    
    [f1, f2] = localBinaryFit(Img, phi, KI, KONE, Ksigma, epsilon);
    
    
    %%% compute lambda1*e1-lambda2*e2
    s1=lambda1.*f1.^2-lambda2.*f2.^2;                   % compute lambda1*e1-lambda2*e2 in the 1st term in eq. (15) in IEEE TIP 08
    s2=lambda1.*f1-lambda2.*f2;
    dataForce=(lambda1-lambda2)*KONE.*Img.*Img+conv2(s1,Ksigma,'same')-2.*Img.*conv2(s2,Ksigma,'same');
    % eq.(15)
    A=-DrcU.*dataForce;                                 % 1st term in eq. (15)
    L=DrcU.*K;                                      % 2nd term in eq. (15)
    L=nu.*L;

    phi = phi + dt*(mu*T_mbe + A + L);
    
    % L
%     L_coef = 1/dt + mu*3/4*delta*Lap2_FFT + mu*Lap_FFT;
%     r_ = phi/dt - mu*delta/4*phi_lap2 + mu*A_phi + A + L;
%     r_ = phi/dt - mu*delta/4*phi_lap2 + mu*A_phi;
    
    % fftn and ifftn
%     r_fftn = fftn(r_);
%     Phi = r_fftn ./ L_coef;
%     phi = real(ifftn(Phi));
    
end

function [f1, f2]= localBinaryFit(Img, u, KI, KONE, Ksigma, epsilon)
% compute f1 and f2
Hu=0.5*(1+(2/pi)*atan(u./epsilon));                     % eq.(8)

I=Img.*Hu;
c1=conv2(Hu,Ksigma,'same');
c2=conv2(I,Ksigma,'same');                              % the numerator of eq.(14) for i = 1
f1=c2./(c1);                                            % compute f1 according to eq.(14) for i = 1
f2=(KI-c2)./(KONE-c1);                                  % compute f2 according to the formula in Section IV-A,
                                                        % which is an equivalent expression of eq.(14) for i = 2.

function k = curvature_central(u)
% compute curvature
[ux,uy] = gradient(u);
normDu = sqrt(ux.^2+uy.^2+1e-10);                       % the norm of the gradient plus a small possitive number
                                                        % to avoid division by zero in the following computation.
Nx = ux./normDu;
Ny = uy./normDu;
[nxx,junk] = gradient(Nx);
[junk,nyy] = gradient(Ny);
k = nxx+nyy;                                            % compute divergence

