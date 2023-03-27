function plotQuiver(ht,tt,vx,vz)
%Plot flow Velocity with Equipotential Lines
global L M N T dt
xx=linspace(0,L,N);
zz=linspace(0,T,M);
[x_dom,z_dom]=meshgrid(xx,zz);
quiver(x_dom,z_dom,vx,vz,'b','Linewidth',1,'AutoScaleFactor',2);
hold on
[C,ht]=contour(x_dom,z_dom,ht);
hold off
clabel(C,ht,'Linestyle','--');
colorbar();
xlabel('Length [m]');
ylabel('Height of seepage [m]');
time = strcat('time = ', num2str(tt*dt),' sec');
text(1,1,time);
MV(tt)=getframe(); %movie 
title ("Flow velocity with Equipotential lines")
end