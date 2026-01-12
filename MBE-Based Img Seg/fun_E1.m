
function eng_E1 =fun_E1(phi,Img,pa)

%     Hphi = 0.5*(1 + 2/pi*atan(phi/pa.epsilon));
%    
%     M1 = Hphi;
%     M2 = 1 - Hphi;
if phi>=0,M1=1;
else M1=0;
end

if phi>=0,M2=0;
else M2=1;
end
    
    [f1, f2]= localBinaryFit(Img, phi,pa);
    Seg1 = pa.lambda1*(Img.^2.*M1.*pa.KONE - 2*Img.*M1.*conv2(f1,pa.K,'same') + M1.*conv2(f1.^2,pa.K,'same'));
    Seg2 = pa.lambda2*(Img.^2.*M2.*pa.KONE - 2*Img.*M2.*conv2(f2,pa.K,'same') + M2.*conv2(f2.^2,pa.K,'same'));
    Seg = Seg1 + Seg2;
  
    
    % Length Term
    [gx, gy] = gradient(phi);
    gnorm = sqrt(gx.^2 + gy.^2);
    

    
    
    DrcPhi=(pa.epsilon/pi)./(pa.epsilon^2.+phi.^2);
    Length = pa.nu*DrcPhi.*gnorm;
 
  
    % MBE Term
    MBE =pa.mu*0.25*(gnorm.^2 - 1).^2;

    eng_E1 = sum(Seg(:))+sum(Length(:))+ sum(MBE(:));
    

    
                    