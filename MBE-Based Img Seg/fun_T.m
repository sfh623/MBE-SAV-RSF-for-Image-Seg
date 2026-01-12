function T = fun_T(phi)
%  phi_lap = phi([end 1:end-1],:) + phi([2:end 1],:) + phi(:,[end 1:end-1]) ...
%          + phi(:,[2:end 1]) - 4*phi;
%     phi_lap2 = phi_lap([end 1:end-1],:) + phi_lap([2:end 1],:) ...
%         + phi_lap(:,[end 1:end-1]) + phi_lap(:,[2:end 1]) - 4*phi_lap;
%  [gx, gy] = gradient(phi);
  gx = (phi([2:end 1],:)-phi([end 1:end-1],:))/2;
  gy = (phi(:,[2:end 1])-phi(:,[end 1:end-1]))/2;

    phi_g2 =(gx.^2 + gy.^2);
    d_x = (phi_g2-1).*gx;
    d_y = (phi_g2-1).*gy;

    dd_x = (-d_x([end 1:end-1],:) + d_x([2:end 1],:))/2;
    dd_y = (-d_y(:,[end 1:end-1]) + d_y(:,[2:end 1]))/2;
%     A_phi = dd_x + dd_y;
A_phi = dd_x + dd_y;
    T=A_phi;
%     T = A_phi - phi_lap;