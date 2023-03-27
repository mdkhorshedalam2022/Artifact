function [hplot,htplot,hpplot,Kxx,Kzz,Kxz,Kzx]=mapping(h,he,kxx,kzz,kxz,kzx,tt)
%mapping total head, pressure head and hydarulic conductivities 1D to 2D
%mapping total head 1D to 3D for specific time 
global M N htplot hpplot
    for ii=1:M
        for jj=1:N
            hplot(ii,jj)=h((ii-1)*N+jj);%mapping total head in 2D
            Kxx(ii,jj)=kxx((ii-1)*N+jj);%mapping kxx in 2D
            Kzz(ii,jj)=kzz((ii-1)*N+jj);%mapping kzz in 2D
            Kxz(ii,jj)=kxz((ii-1)*N+jj);%mapping kxz in 2D
            Kzx(ii,jj)=kzx((ii-1)*N+jj);%mapping kzx in 2D
            htplot(ii,jj,tt)=h((ii-1)*N+jj);%mapping total head in 3D
            hpplot(ii,jj,tt)=h((ii-1)*N+jj)-he((ii-1)*N+jj);%mapping pressure head in 3D
        end
    end
end