%Inversion modelling Magnetic Crustal block fix
% using Simulated Annealing  algorithm
%Mohammad Rheza Zamani
%Reference : Grandis,H.(2009): Pengantar pemodelan inversi geofisika, Jakarta:HAGI
clear all;
clc;
%Input Paramater
mu0 = 1.26*10^-3; %H/km (Permeabilitas magnetik)
delta_Mz = 1*10^-3; %Kontras densitas
z1 = 2.5; %km
z2 = 3; %km
x0 = 20; %km
m = 20; %km
%titik pengukuran
x = -100:1:100;
%data sintetik
[d_Bzobs] = fwd_mgt(x,x0,z1,z2,mu0,delta_Mz,m);

%Paramater inversi
nitr = 500;
T = 5;
dec = 0.05;
%Batas atas dan batas bawah
x0_min = 1;
x0_max = 30;
z1_min = 1;
z1_max = 5;
z2_min = 1;
z2_max = 5;
delta_Mz_min = 1*10^-3;
delta_Mz_max = 5*10^-3;
m_min = 1;
m_max = 30;

model1(1,1) = 20;
model1(1,2) = 3;
model1(1,3) = 2;
model1(1,4) = 2;
model1(1,5) = 20;
%model1(1,1) = x0_min + rand*(x0_max-x0_min);
%model1(1,2) = z1_min + rand*(z1_max-z1_min);
%model1(1,3) = z1_min + rand*(z2_max - z2_min);
%model1(1,4) = delta_Mz_min + rand*(delta_Mz_max - delta_Mz_min);
%model1(1,5) = m_min + rand*(m_max-m_min);
dBz1(1,:) = fwd_mgt(x,model1(1),model1(2),model1(3),mu0,model1(4),model1(5));
misfit1 = misfit_magnetic(d_Bzobs,dBz1);
%Membuat Video
v = VideoWriter('SA Inversi Magnetic.avi');
open(v);
for itr = 1 : nitr
    model2(1,1) = x0_min + rand*(x0_max-x0_min);
    model2(1,2) = z1_min + rand*(z1_max-z1_min);
    model2(1,3) = z2_min + rand*(z2_max - z2_min);
    model2(1,4) = delta_Mz_min + rand*(delta_Mz_max - delta_Mz_min);
    model2(1,5) = m_min + rand*(m_max-m_min);
    dBz2(1,:) = fwd_mgt(x,model2(1),model2(2),model2(3),mu0,model2(4),model2(5));
    misfit2 = misfit_magnetic(d_Bzobs(1,:),dBz2(1,:));
    delta_E = misfit2 - misfit1;
    if delta_E<0
        model1 = model2;
        misfit1 = misfit2;
        dBz1 = dBz2;

           if model1(1)<x0_min 
                model1(1) = x0_min;
            end
            if model1(2)<z1_min
                model1(2) = z1_min;
            end
            if model1(3)<z2_min
                model1(3) = z2_min;
            end
            if model1(4)<delta_Mz_min
                model1(4) = delta_Mz_min;
            end
            if model1(5)<m_min
                model1(5) = m_min;
            end
            if model1(1)>x0_max 
                model1(1) = x0_max;
            end
            if model1(2)>z1_max
                model1(2) = z1_max;
            end
            if model1(3)>z2_max
                model1(3) = z2_max;
            end
            if model1(4)>delta_Mz_max
                model1(4) = delta_Mz_max;
            end
            if model1(5)>m_max
                model1(5) = m_max;
            end
    else
        P = exp(-delta_E/T); 
        if P>= rand
           model1 = model2;
           misfit1 = misfit2;
           dBz1 = dBz2;
           if model1(1)<x0_min 
                model1(1) = x0_min;
            end
            if model1(2)<z1_min
                model1(2) = z1_min;
            end
            if model1(3)<z2_min
                model1(3) = z2_min;
            end
            if model1(4)<delta_Mz_min
                model1(4) = delta_Mz_min;
            end
            if model1(5)<m_min
                model1(5) = m_min;
            end
            if model1(1)>x0_max 
                model1(1) = x0_max;
            end
            if model1(2)>z1_max
                model1(2) = z1_max;
            end
            if model1(3)>z2_max
                model1(3) = z2_max;
            end
            if model1(4)>delta_Mz_max
                model1(4) = delta_Mz_max;
            end
            if model1(5)>m_max
                model1(5) = m_max;
            end
        else
           model1 = model1;
           misfit1 = misfit1;
           dBz1 = dBz1;
           if model1(1)<x0_min 
                model1(1) = x0_min;
            end
            if model1(2)<z1_min
                model1(2) = z1_min;
            end
            if model1(3)<z2_min
                model1(3) = z2_min;
            end
            if model1(4)<delta_Mz_min
                model1(4) = delta_Mz_min;
            end
            if model1(5)<m_min
                model1(5) = m_min;
            end
            if model1(1)>x0_max 
                model1(1) = x0_max;
            end
            if model1(2)>z1_max
                model1(2) = z1_max;
            end
            if model1(3)>z2_max
                model1(3) = z2_max;
            end
            if model1(4)>delta_Mz_max
                model1(4) = delta_Mz_max;
            end
            if model1(5)>m_max
                model1(5) = m_max;
            end
    end
    end
   Egen(itr) = misfit1;
   T = T*(1-dec);
   Temperature(itr) = T;
%persiapan plot
x1 = x0 - (m);
x2 = x0 - (m);
x3 = x0 + (m);
x4 = x0 + (m);
y1 = z1;
y2 = z2;
y3 = z2;
y4 = z1;

x11 = model1(1) - model1(5);
x22 = model1(1) - model1(5);
x33 = model1(1) + model1(5);
x44 = model1(1) + model1(5);
y11 = model1(2);
y22 = model1(3);
y33 = model1(3);
y44 = model1(2);

figure(1)
subplot(2,1,1)
plot(x,d_Bzobs*10^9,'b.',x,dBz1*10^9,'r-','MarkerSize',20,'MarkerFaceColor','r','LineWidth',2.5)
xlabel('Distance (km)','FontWeight','bold')
ylabel(['\DeltaB_{z} (nT)'],'FontWeight','bold')
title(['\bf \fontsize{10}\fontname{Times}Geomagnetic model || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)],'FontWeight','bold');
subtitle(['x0 = ',num2str(model1(1)),' || z1 = ',num2str(model1(2)),' || z2 = ',num2str(model1(3)),' || \DeltaM_{z} = ',num2str(model1(4)),' || m = ',num2str(model1(5))],'FontWeight','bold')
legend({'Observed Data','Calculated Data'},'Location','Southeast')
grid on
subplot(2,1,2)
fill([x1 x2 x3 x4],[y1 y2 y3 y4],'k')
hold on
fill([x11 x22 x33 x44],[y11 y22 y33 y44],'r')
hold off
axis([min(x) max(x) 0 5])
xlabel('Distance (km)','FontWeight','bold')
ylabel('Depth (km)','FontWeight','bold')
title('Subsurface Model','FontWeight','bold')
set(gca,'Ydir','reverse')
legend({'Synthetic Model','Inversion Model'},'Location','Southeast')
grid on 
set(gcf, 'Position', get(0, 'Screensize'));
frame = getframe(gcf);
writeVideo(v,frame);
end
close(v);

figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RSME','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Misfit ');
grid on

%Plot Temperature
figure(3)
plot(1:nitr,Temperature,'b','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('Temperature','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Penurunan Temperature ');
grid on
