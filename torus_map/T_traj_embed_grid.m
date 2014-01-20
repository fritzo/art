%  Author: Karl Obermeyer
%  Last modified: 10/31/04
%  This Matlab program plots various trajectories on a T^2
%  given parameter values

clear all;format long;clc;

%Omega_1=.5725;
Omega_1=.5726;
%Omega_2=.81975;
Omega_2=.81958;

%function T_trajectories=T_trajectories()

k=1;
%STEPS=400;            %steps to interate initial condition
STEPS=600;    
%RES_1=35;             %Create IC grid w/specified res
RES_1=70;
%RES_2=35;
RES_2=70;
ORBITS=RES_1*RES_2;   %Number of plot calls
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
        %plot(x(5*(1:STEPS/5),1),x(5*(1:STEPS/5),2),'.');       
        plot(x(:,1,i),x(:,2,i),'.','MarkerSize',1); 
        %xlabel('x_1');ylabel('x_2');
     end
end
 axis off;hold off;
 

 %%figure; hold on; set(gcf,'Position',[200 200 700 600],'Color',[1 1 1]); %initialize next figure;
 %title(sprintf('Omega_1=%1.4f; Omega_2=%1.4f; k=%1.3f; STEPS=%g', Omega_1, Omega_2, k, STEPS) ); 
 
 %Now embed T^2 in R^3 and plot
 %%y=zeros(STEPS,3,ORBITS); %initialize place to store embedded orbits
 %%R=3;r=1; %major, minor radius
 %%x(:,1,:)=2*pi*x(:,1,:); %convert to radians
 %%x(:,2,:)=2*pi*x(:,2,:);
 
 %%y(:,1,:) = R*cos(x(:,1,:))+r*sin(x(:,2,:)).*cos(x(:,1,:));
 %%y(:,2,:) = R*sin(x(:,1,:))+r*sin(x(:,2,:)).*sin(x(:,1,:));
 %%y(:,3,:) = r*cos(x(:,2,:));
 %%for i=1:ORBITS
 %%    plot3(y(:,1,i),y(:,2,i),y(:,3,i),'.','MarkerSize',4); 
 %%end
 %%axis([-4 4 -4 4 -1.5 1.5]);view(3);axis off;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
