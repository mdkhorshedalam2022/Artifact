function [coeff_1,coeff_2,coeff_3,coeff_4,coeff_5,coeff_6,coeff_7,coeff_8,coeff_9,coeff_10]=coreEqns(kxx,kxx_1,kzz,kzz_N,kxz,kxz_1,kzx,kzx_N,dx,dz,dt,mv,h)
%cofficients of core equations
coeff_1=kxx/dx^2; %a(ii,ii-1)
coeff_2=-(((kxx_1-kxx)/dx^2)+2*(kxx/dx^2)+((kxz_1-kxz)/(dx*dz))+((kzx_N-kzx)/(dz*dx))+((kzz_N-kzz)/dz^2)+2*(kzz/dz^2)+mv/dt); %a(ii,ii)
coeff_3=(((kxx_1-kxx)/dx^2)+(kxx/dx^2)+((kzx_N-kzx)/(dz*dx)));%a(ii,ii+1)
coeff_4=kzz/dz^2; %a(ii,ii-N)
coeff_5=(((kxz_1-kxz)/(dx*dz))+((kzz_N-kzz)/dz^2)+(kzz/dz^2)); %a(ii,ii+N)
coeff_6=(1/4)*((kxz/(dx*dz))+ (kzx/(dz*dx))); %a(ii,ii+N+1)
coeff_7=-(1/4)*((kxz/(dx*dz))+ (kzx/(dz*dx))); %a(ii,ii-N+1)
coeff_8=-(1/4)*((kxz/(dx*dz))+ (kzx/(dz*dx))); %a(ii,ii+N-1)
coeff_9=(1/4)*((kxz/(dx*dz))+ (kzx/(dz*dx))); %a(ii,ii-N-1)
coeff_10=-(mv/dt)*h; % b(ii)
end