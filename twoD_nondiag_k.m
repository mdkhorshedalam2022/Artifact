function [kx,kz,kxz,kzx]=twoD_nondiag_k(k1,k3,alpha)
% Calculate 2D nondiagonalized hydraulic conductivity tensor (k) elements
theta=deg2rad(alpha); %converting angle into radians
kx=((k1+k3)/2)+(((k1-k3)/2)*cos(2*theta));
kz=((k1+k3)/2)-(((k1-k3)/2)*cos(2*theta));
kxz=(((k1-k3)/2)*sin(2*theta)); 
kzx=kxz; 
end