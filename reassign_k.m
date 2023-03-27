function [kxxold, kzzold, kxzold, kzxold]=reassign_k(kxx,kzz,kxz,kzx)
%storing previous values of the elements of k as old k 
%before calucating new k
        kxxold=kxx; 
        kzzold=kzz;
        kxzold=kxz;
        kzxold=kzx;
end
        