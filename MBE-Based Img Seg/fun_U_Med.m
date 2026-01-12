function phi = fun_U_Med(phi,Img,pa)
   
 T = fun_T(phi);
    
% RSF
    M=curvature_central(phi);
    DrcU=(pa.epsilon/pi)./(pa.epsilon^2.+phi.^2);               % eq.(9)
    
    [f1, f2] = localBinaryFit(Img, phi,pa );
    
   
    s1=pa.lambda1.*f1.^2-pa.lambda2.*f2.^2;                   % compute lambda1*e1-lambda2*e2 in the 1st term in eq. (15) in IEEE TIP 08
    s2=pa.lambda1.*f1-pa.lambda2.*f2;
%     dataForce=(pa.lambda1-pa.lambda2)*pa.KONE.*Img.*Img+conv2(s1,pa.K,'same')-2.*Img.*conv2(s2,pa.K,'same');
   dataForce=(pa.lambda1-pa.lambda2)*pa.KONE.*Img.*Img...
        +imfilter(s1,pa.K,'replicate')...
        -2.*Img.*imfilter(s2,pa.K,'replicate');%for medical
    A=DrcU.*dataForce;                                
    L=DrcU.*M;                                      
    L1=pa.nu.*L;
    phi =(-pa.mu*T - L1+A);


