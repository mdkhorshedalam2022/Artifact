function [kxx,kzz,kxz,kzx]=calculate_k(kxx0,kzz0,kxz0,kzx0,h_old,h,he)
 %calculate components of hydraulic conductivity tensor 
 %from intial values each component based on saturation of soil
global a1 a2
    if h-he>=0 %for saturated soil
        kxx= kxx0;
        kzz= kzz0;
        kxz= kxz0;
        kzx= kzx0;
    else %for unsaturated soil
        %water rentention formula
        % kxx = kxx0./(1+a1*(abs((h_old+h_old)/2-he)).^a2);
        % kzz = kzz0./(1+a1*(abs((h_old+h_old)/2-he)).^a2);
        % kxz = kxz0./(1+a1*(abs((h_old+h_old)/2-he)).^a2);
        % kzx = kzx0./(1+a1*(abs((h_old+h_old)/2-he)).^a2);  
        kxx = kxx0./(1+a1*(abs((h_old+h)/2-he)).^a2);
        kzz = kzz0./(1+a1*(abs((h_old+h)/2-he)).^a2);
        kxz = kxz0./(1+a1*(abs((h_old+h)/2-he)).^a2);
        kzx = kzx0./(1+a1*(abs((h_old+h)/2-he)).^a2); 
    end
end