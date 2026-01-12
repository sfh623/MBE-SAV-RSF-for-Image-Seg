 function A= fun_A(phi,pa)
phi_lap = phi([end 1:end-1],:) + phi([2:end 1],:) + phi(:,[end 1:end-1]) ...
        + phi(:,[2:end 1]) - 4*phi;
   phi_lap2 = phi_lap([end 1:end-1],:) + phi_lap([2:end 1],:) ...
        + phi_lap(:,[end 1:end-1]) + phi_lap(:,[2:end 1]) - 4*phi_lap;

    A =eye(size(pa.dt.*phi_lap2))+pa.dt.*phi_lap2.*delta;

