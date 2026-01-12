function [f1, f2]= localBinaryFit2(Img, u, pa)
if u>=0,Hu=1;
else Hu=0;
end
% eq.(8)
I=Img.*Hu;
c1=conv2(Hu,pa.K,'same');
c2=conv2(I,pa.K,'same');                              % the numerator of eq.(14) for i = 1
f1=c2./(c1);                                            % compute f1 according to eq.(14) for i = 1
f2=(pa.KI-c2)./(pa.KONE-c1);  
