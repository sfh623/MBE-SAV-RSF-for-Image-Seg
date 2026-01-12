% Split Algorithm
% Fourier for the MBE.
function phi = Reg_MBE_Fourier(phi, delta, dt, iterN)

%% FFT prepare
[M, N] = size(phi);
ii = 1:1:M;
ii = ii';
jj = 1:1:N;
Z1 = 2*pi*(ii-1)/M;
Z2 = 2*pi*(jj-1)/N;
Z1_c = cos(Z1);
Z2_c = cos(Z2);
Z1_2 = repmat(Z1_c, 1, N);
Z2_2 = repmat(Z2_c, M, 1);
Lap2_FFT = 4 * (Z1_2 + Z2_2 - 2).^2;
Lap_FFT = 2 * (Z1_2 + Z2_2 - 2);


%% Iteration
for i=1:iterN
    % L
    L_coef = 1/dt + 3/4*delta*Lap2_FFT + 2*Lap_FFT;
    
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
    
    r_ = phi/dt - delta/4*phi_lap2 + A_phi;
    
    % fftn and ifftn
    r_fftn = fftn(r_);
    Phi = r_fftn ./ L_coef;
    phi = real(ifftn(Phi));
    
end
end