%  Author: Karl Obermeyer
%  Last modified: 10/31/04
%  This Matlab program plots various trajectories on a T^2
%  given parameter values
clear all;format long;clc;

%parameter values
Omega_1=.5726;
Omega_2=.81958;
k=1;

STEPS=600;               %steps to interate initial conditions
RES_1=70;                %Create IC grid w/specified res
RES_2=70;
ORBITS=RES_1*RES_2;      %Number of plot calls
X_1=linspace(0,1,RES_1); X_2=linspace(0,1,RES_2);
x=zeros(STEPS,2,ORBITS); %initialize place to store orbits

figure; axis([0 1 0 1]); hold on; set(gcf,'Position',[200 200 700 600],'Color',[1 1 1]);
%title(sprintf('Omega_1=%1.4f; Omega_2=%1.4f; k=%1.3f; STEPS=%g', Omega_1, Omega_2, k, STEPS) ); 
%Get initial conditions and plot orbits on phase space
i=0;
for m=1:RES_1
    for n=1:RES_2   
        i=i+1;
        x_0=[X_1(m) X_2(n)];
        for j=2:STEPS 
            x(1,:,i)=x_0;   %reset initial condition
            x(j,1,i)=mod(x(j-1,1,i)+Omega_1+(k/(2*pi))*sin(2*pi*x(j-1,2,i)),1);
            x(j,2,i)=mod(x(j-1,2,i)+Omega_2+(k/(2*pi))*sin(2*pi*x(j-1,1,i)),1); 
        end
        plot(x(:,1,i),x(:,2,i),'.','MarkerSize',1);     
     end
end
 axis off;hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
