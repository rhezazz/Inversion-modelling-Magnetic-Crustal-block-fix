%Inversion modelling Magnetic Crustal block fix
% using Genetic Algorithm
%Reference : Grandis,H.(2009): Pengantar pemodelan inversi geofisika, Jakarta:HAGI
clear all;
clc;
%Input Paramater
mu0 = 1.26*10^-3; %H/km (Permeabilitas magnetik)
delta_Mz = 1; %Kontras densitas
z1 = 2.5; %km
z2 = 3; %km
x0 = 20; %km
m = 20; %km
%titik pengukuran
x = -100:1:100;
%data sintetik
[d_Bzobs] = fwd_mgt(x,x0,z1,z2,mu0,delta_Mz,m);

%Parameter inversi
prob = 0.9;
nitr = 100;
npop = 200;
%Batas atas dan batas bawah
x0_min = 1;
x0_max = 30;
z1_min = 1;
z1_max = 5;
z2_min = 1;
z2_max = 5;
delta_Mz_min = 1;
delta_Mz_max = 5;
m_min = 1;
m_max = 30;
%Membuat model acak
for i = 1 : npop
    model(i,1) = x0_min + rand*(x0_max-x0_min);
    model(i,2) = z1_min + rand*(z1_max-z1_min);
    model(i,3) = z1_min + rand*(z2_max - z2_min);
    model(i,4) = delta_Mz_min + rand*(delta_Mz_max - delta_Mz_min);
    model(i,5) = m_min + rand*(m_max-m_min);
    dBz(i,:) = fwd_mgt(x,model(1),model(2),model(3),mu0,model(4),model(5));
    misfit(i) = misfit_magnetic(d_Bzobs,dBz(i,:));
end

%Membuat Video
v = VideoWriter('GA Inversi Magnetic.avi');
open(v);

%Inversion
for itr = 1 : nitr
    %fitness
    fitness = 1./misfit;
    %fitness ternomalisasi
    fitness_norm = fitness./sum(fitness);
    sc = 0;
    for i = 1 : npop
        sc = sc + fitness_norm(i);
        cumm(i) = sc;
    end
    indx = 1;
    %Roulette wheel pemilihan induk
    for i=1:npop/2
      R1 = rand;
      for j=1:npop
        if R1 < cumm(j)
            ipar1 = j;
            break
        end
    end
    
    R2 = rand;
    for j=1:npop
        if R2 < cumm(j)
            ipar2 = j;
            break
        end
    end
    
    %Crossover dan offspring
    R3 = rand;
    if R3 < prob
        i1 = rand;
        i2 = rand;
        i3 = rand;
        i4 = rand;
        i5 =  rand;
        model_new(indx,1) = i1*model(ipar1,1)+(1-i1)*model(ipar2,1);
        model_new(indx+1,1) = i1*model(ipar2,1)+(1-i1)*model(ipar1,1);
        model_new(indx,2) =  i2*model(ipar1,2)+(1-i2)*model(ipar2,2);
        model_new(indx+1,2) = i2*model(ipar2,2)+(1-i2)*model(ipar1,2);
        model_new(indx,3) = i3*model(ipar1,3)+(1-i3)*model(ipar2,3);
        model_new(indx+1,3) = i3*model(ipar2,3)+(1-i3)*model(ipar1,3);
        model_new(indx,4) = i4*model(ipar1,4)+(1-i4)*model(ipar2,4);
        model_new(indx+1,4) = i4*model(ipar2,4)+(1-i4)*model(ipar1,4);
        model_new(indx,5) = i5*model(ipar1,5)+(1-i5)*model(ipar2,5);
        model_new(indx+1,5) = i5*model(ipar2,5)+(1-i5)*model(ipar1,5);
    else
        model_new(indx,:) = model(ipar1,:);
        model_new(indx+1,:) = model(ipar2,:);
    end
    %Menambah index untuk setiap iterasi
     indx = indx + 2; 
    end
    model = model_new;
    for i=1:npop
        dBz_new(i,:) = fwd_mgt(x,model(i,1),model(i,2),model(i,3),mu0,model(i,4),model(i,5));
        misfit(i) = misfit_magnetic(d_Bzobs(1,:),dBz_new(i,:));
    end
    Egen(itr) = misfit(i);
    fitgen(itr) = fitness(i);


%persiapan plot
x1 = x0 - (m);
x2 = x0 - (m);
x3 = x0 + (m);
x4 = x0 + (m);
y1 = z1;
y2 = z2;
y3 = z2;
y4 = z1;
pgon = polyshape([x1 x2 x3 x4],[y1 y2 y3 y4]);
x11 = model(i,1) - model(i,5);
x22 = model(i,1) - model(i,5);
x33 = model(i,1) + model(i,5);
x44 = model(i,1) + model(i,5);
y11 = model(i,3);
y22 = model(i,4);
y33 = model(i,4);
y44 = model(i,3);
pgon_mod = polyshape([x11 x22 x33 x44],[y11 y22 y33 y44]);
figure(1)
subplot(2,1,1)
plot(x,d_Bzobs,'b.',x,dBz_new,'r-','MarkerSize',20,'MarkerFaceColor','r','LineWidth',2.5)
xlabel('Distance (km)','FontWeight','bold')
ylabel(['\DeltaB_{z} (nT)'],'FontWeight','bold')
title(['\bf \fontsize{10}\fontname{Times}Geomagnetic model || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)],'FontWeight','bold');
subtitle(['x0 = ',num2str(model(i,1)),' || z1 = ',num2str(model(i,2)),' || z2 = ',num2str(model(i,3)),' || \DeltaM_{z} = ',num2str(model(i,4)),' || m = ',num2str(model(i,5))],'FontWeight','bold')
legend({'Observed Data','Calculated Data'},'Location','Southeast')
grid on
subplot(2,1,2)
plot(pgon,'Facecolor','k')
hold on
plot(pgon_mod,'Facecolor','r')
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

%Plot grafik misfit
figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RSME','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Misfit ');
grid on

%Plot grafik fitness
figure(3)
plot(1:nitr,fitgen,'b','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('Fitness','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Fitness ');
grid on
