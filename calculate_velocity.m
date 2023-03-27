function [vx,vy]= calculate_velocity(hplot,Kx,Ky,Kxy,Kyx,dx,dy)
%Calculate velocity components
    global M N
    
   for ii=1
       for jj=1:N
           if jj==1 % bottom left node (x&y both forward difference)
               vx(ii,jj)=-Kx(ii,jj)*(hplot(ii,jj+1)-hplot(ii,jj))/dx-Kxy(ii,jj)*(hplot(ii+1,jj)-hplot(ii,jj))/dy; 
               vy(ii,jj)=-Kyx(ii,jj)*(hplot(ii,jj+1)-hplot(ii,jj))/dx-Ky(ii,jj)*(hplot(ii+1,jj)-hplot(ii,jj))/dy;
           elseif jj==N % bottom right node (x backward &y forward difference)
               vx(ii,jj)=-Kx(ii,jj)*(hplot(ii,jj)-hplot(ii,jj-1))/dx-Kxy(ii,jj)*(hplot(ii+1,jj)-hplot(ii,jj))/dy; 
               vy(ii,jj)=-Kyx(ii,jj)*(hplot(ii,jj)-hplot(ii,jj-1))/dx-Ky(ii,jj)*(hplot(ii+1,jj)-hplot(ii,jj))/dy;
           else % bottom (x central and y forward) 
               vx(ii,jj)=-Kx(ii,jj)*(hplot(ii,jj+1)-hplot(ii,jj-1))/(2*dx)-Kxy(ii,jj)*(hplot(ii+1,jj)-hplot(ii,jj))/dy; 
               vy(ii,jj)=-Kyx(ii,jj)*(hplot(ii,jj+1)-hplot(ii,jj-1))/(2*dx)-Ky(ii,jj)*(hplot(ii+1,jj)-hplot(ii,jj))/dy;    
           end
       end
   end
   
   for ii=M
       for jj=1:N
           if jj==1 % top left node (x forward &y backward difference)
               vx(ii,jj)=-Kx(ii,jj)*(hplot(ii,jj+1)-hplot(ii,jj))/dx-Kxy(ii,jj)*(hplot(ii,jj)-hplot(ii-1,jj))/dy; 
               vy(ii,jj)=-Kyx(ii,jj)*(hplot(ii,jj+1)-hplot(ii,jj))/dx-Ky(ii,jj)*(hplot(ii,jj)-hplot(ii-1,jj))/dy;
           elseif jj==N % top right node (x and y backward difference)
               vx(ii,jj)=-Kx(ii,jj)*(hplot(ii,jj)-hplot(ii,jj-1))/dx-Kxy(ii,jj)*(hplot(ii,jj)-hplot(ii-1,jj))/dy; 
               vy(ii,jj)=-Kyx(ii,jj)*(hplot(ii,jj)-hplot(ii,jj-1))/dx-Ky(ii,jj)*(hplot(ii,jj)-hplot(ii-1,jj))/dy;
           else % top (x central and y backward) 
               vx(ii,jj)=-Kx(ii,jj)*(hplot(ii,jj+1)-hplot(ii,jj-1))/(2*dx)-Kxy(ii,jj)*(hplot(ii,jj)-hplot(ii-1,jj))/dy; 
               vy(ii,jj)=-Kyx(ii,jj)*(hplot(ii,jj+1)-hplot(ii,jj-1))/(2*dx)-Ky(ii,jj)*(hplot(ii,jj)-hplot(ii-1,jj))/dy;    
           end
       end
   end
   
   for ii=2:M-1
       for jj=1:N
           if jj==1 % left side (x forward &y central)
               vx(ii,jj)=-Kx(ii,jj)*(hplot(ii,jj+1)-hplot(ii,jj))/dx-Kxy(ii,jj)*(hplot(ii+1,jj)-hplot(ii-1,jj))/(2*dy); 
               vy(ii,jj)=-Kyx(ii,jj)*(hplot(ii,jj+1)-hplot(ii,jj))/dx-Ky(ii,jj)*(hplot(ii+1,jj)-hplot(ii-1,jj))/(2*dy);
           elseif jj==N % right side (x backward &y central)
               vx(ii,jj)=-Kx(ii,jj)*(hplot(ii,jj)-hplot(ii,jj-1))/dx-Kxy(ii,jj)*(hplot(ii+1,jj)-hplot(ii-1,jj))/(2*dy); 
               vy(ii,jj)=-Kyx(ii,jj)*(hplot(ii,jj)-hplot(ii,jj-1))/dx-Ky(ii,jj)*(hplot(ii+1,jj)-hplot(ii-1,jj))/(2*dy);
           else % core(central difference)
               vx(ii,jj)=-Kx(ii,jj)*(hplot(ii,jj+1)-hplot(ii,jj-1))/2/dx-Kxy(ii,jj)*(hplot(ii+1,jj)-hplot(ii-1,jj))/2/dy; 
               vy(ii,jj)=-Kyx(ii,jj)*(hplot(ii,jj+1)-hplot(ii,jj-1))/2/dx-Ky(ii,jj)*(hplot(ii+1,jj)-hplot(ii-1,jj))/2/dy; 
           end
       end
   end
       
       
   
   
