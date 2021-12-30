%Inversion modelling Magnetic Crustal block fix
% using  Very Fast Simulated Annealing  algorithm (VFSA)
%Mohammad Rheza Zamani
%Reference : Grandis,H.(2009): Pengantar pemodelan inversi geofisika, Jakarta:HAGI
%Reference : W. Srigutomo, M. Heriyanto, and M. Hilmi Aufa. Gravity Inversion of Talwani Model using Very Fast Simulated Annealing. Journal of Mathematical and Fundamental Sciences, Vol. 51, No. 2, 2019, 177-190. doi: 10.5614/j.math.fund.sci.2019.51.2.7.
clear all;
clc;
%Input Paramater
mu0 = 1.26*10^-3; %H/km (Permeabilitas magnetik)
delta_Mz = 1; %Kontras densitas
z1 = 2.5; %km
z2 = 3; %km
x0 = 2; %km
m = 20; %km
%titik pengukuran
x = -100:1:100;
%data sintetik
[d_Bzobs] = fwd_mgt(x,x0,z1,z2,mu0,delta_Mz,m);

%Paramater inversi
nitr = 200;
T = 1;
dec =1;
%Batas atas dan batas bawah
x0_min = 1;
x0_max = 5;
z1_min = 1;
z1_max = 10;
z2_min = 1;
z2_max = 10;
delta_Mz_min = 1;
delta_Mz_max = 10;
m_min = 1;
m_max = 50;

%Tebakan awal
%model(1,1) = 20;
%model(1,2) =7;
%model(1,3) = 7;
%model(1,4) = 5;
%model(1,5) = 20;
model(1,1) = x0_min + rand*(x0_max-x0_min);
model(1,2) = z1_min + rand*(z1_max-z1_min);
model(1,3) = z2_min + rand*(z2_max - z2_min);
model(1,4) = delta_Mz_min + rand*(delta_Mz_max - delta_Mz_min);
model(1,5) = m_min + rand*(m_max-m_min);
dBz(1,:) = fwd_mgt(x,model(1),model(2),model(3),mu0,model(4),model(5));
misfit = misfit_magnetic(d_Bzobs,dBz);
%Membuat Video
v = VideoWriter('VFSA Inversi Magnetic.avi');
open(v);
for itr = 1 : nitr
    model1(1,1) = x0_min + rand*(x0_max-x0_min);
    model1(1,2) = z1_min + rand*(z1_max-z1_min);
    model1(1,3) = z2_min + rand*(z2_max - z2_min);
    model1(1,4) = delta_Mz_min + rand*(delta_Mz_max - delta_Mz_min);
    model1(1,5) = m_min + rand*(m_max-m_min);
    ui = rand;
    yi = sign(ui-0.5)*T*((((1 + (1/T)))^abs(2*ui-1))-1);
    model2(1,1) = model1(1,1) + yi*(x0_max-x0_min);
    model2(1,2) = model1(1,2) + yi*(z1_max-z1_min);
    model2(1,3) = model1(1,3) + yi*(z2_max - z2_min);
    model2(1,4) = model1(1,4) + yi*(delta_Mz_max - delta_Mz_min);
    model2(1,5) = model1(1,5) + yi*(m_max-m_min);
    dBz2(1,:) = fwd_mgt(x,model2(1),model2(2),model2(3),mu0,model2(4),model2(5));
    misfit2 = misfit_magnetic(d_Bzobs(1,:),dBz2(1,:));
    delta_E = misfit2 - misfit;
    if delta_E<= 0
        model = model2;
        misfit = misfit2;
        dBz = dBz2;

           if model(1)<x0_min 
                model(1) = x0_min;
            end
            if model(2)<z1_min
                model(2) = z1_min;
            end
            if model(3)<z2_min
                model(3) = z2_min;
            end
            if model(4)<delta_Mz_min
                model(4) = delta_Mz_min;
            end
            if model(5)<m_min
                model(5) = m_min;
            end
            if model(1)>x0_max 
                model(1) = x0_max;
            end
            if model(2)>z1_max
                model(2) = z1_max;
            end
            if model(3)>z2_max
                model(3) = z2_max;
            end
            if model(4)>delta_Mz_max
                model(4) = delta_Mz_max;
            end
            if model(5)>m_max
                model(5) = m_max;
            end
    else
        P = exp(-delta_E/T); 
        R1 = rand;
        if P> R1
           model = model2;
           misfit = misfit2;
           dBz = dBz2;
           if model(1)<x0_min 
                model(1) = x0_min;
            end
            if model(2)<z1_min
                model(2) = z1_min;
            end
            if model(3)<z2_min
                model(3) = z2_min;
            end
            if model(4)<delta_Mz_min
                model(4) = delta_Mz_min;
            end
            if model(5)<m_min
                model(5) = m_min;
            end
            if model(1)>x0_max 
                model(1) = x0_max;
            end
            if model(2)>z1_max
                model(2) = z1_max;
            end
            if model(3)>z2_max
                model(3) = z2_max;
            end
            if model(4)>delta_Mz_max
                model(4) = delta_Mz_max;
            end
            if model(5)>m_max
                model(5) = m_max;
            end
    end
    end
   Egen(itr) = misfit;
   T = T*exp(-dec*(itr)^(1/length(model)));
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

x11 = model(1) - model(5);
x22 = model(1) - model(5);
x33 = model(1) + model(5);
x44 = model(1) + model(5);
y11 = model(3);
y22 = model(4);
y33 = model(4);
y44 = model(3);

figure(1)   
subplot(2,1,1)
plot(x,d_Bzobs,'b.',x,dBz,'r-','MarkerSize',20,'MarkerFaceColor','r','LineWidth',2.5)
xlabel('Distance (km)','FontWeight','bold')
ylabel(['\DeltaB_{z} (nT)'],'FontWeight','bold')
title(['\bf \fontsize{10}\fontname{Times}Geomagnetic model || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)],'FontWeight','bold');
subtitle(['x0 = ',num2str(model(1)),' || z1 = ',num2str(model(2)),' || z2 = ',num2str(model(3)),' || \DeltaM_{z} = ',num2str(model(4)),' || m = ',num2str(model(5))],'FontWeight','bold')
legend({'Observed Data','Calculated Data'},'Location','Southeast')
set(gcf, 'Position', get(0, 'Screensize'));
grid on
subplot(2,1,2)
fill([x1 x2 x3 x4],[y1 y2 y3 y4],'k')
hold on
fill([x11 x22 x33 x44],[y11 y22 y33 y44],'r')
hold off
axis([min(x) max(x) 0 10])
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
