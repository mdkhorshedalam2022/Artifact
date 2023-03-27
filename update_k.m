function [mv,kxx,kzz,kxz,kzx]=update_k(h,he,h_old,kxx0,kzz0,kxz0,kzx0)
%Updating hydraulic conductivity
global a1 a2
    if h-he>=0 %for saturated soil
        mv=0.00001;
        kxx= kxx0;
        kzz= kzz0;
        kxz= kxz0;
        kzx= kzx0;
    else
        mv=0.001;%for unsaturated soil
        D=(1+a1*(abs((h+h_old)/2-he))^a2);
        kxx = kxx0/D; kzz = kzz0/D; 
        kxz = kxz0/D; kzx = kzx0/D;           
    end

end