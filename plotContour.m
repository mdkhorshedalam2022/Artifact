function plotContour(ht,tt)
%Plot contour lines of pressure head
global L M N T dt
xx=linspace(0,L,N);
yy=linspace(0,T,M);
[x_dom,y_dom]=meshgrid(xx,yy);
[C,ht]=contourf(x_dom,y_dom,ht);
clabel(C,ht,'LineStyle','--');
colorbar();
xlabel('Length [m]');
ylabel('Height of pressure head [m]');
time = strcat('time = ', num2str(dt*tt),' sec');
text(1,1,time);
MV(tt)=getframe(); %movie 
title ("Contour lines of Pressure head")
end

