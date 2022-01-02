%Inversion modelling Magnetic Crustal block fix
% using modified Symbiotic Organism Search
%Mohammad Rheza Zamani
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
[dBdz_obs] = fwd_mgt(x,x0,z1,z2,mu0,delta_Mz,m);

%Paramater inversi
nitr = 500;
npop = 100;
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

for i = 1 : npop
    model(i,1) = x0_min + rand*(x0_max-x0_min);
    model(i,2) = z1_min + rand*(z1_max-z1_min);
    model(i,3) = z2_min + rand*(z2_max - z2_min);
    model(i,4) = delta_Mz_min + rand*(delta_Mz_max - delta_Mz_min);
    model(i,5) = m_min + rand*(m_max-m_min);
    dBdz(i,:) = fwd_mgt(x,model(1),model(2),model(3),mu0,model(4),model(5));
    E(i) = misfit_magnetic(dBdz_obs,dBdz);
end

%Proses Inversi
for itr = 1 : nitr
    for i = 1 : npop
        idx = find(E ==min(E));
        model_best = model(idx(1),:);
        %Mutualisme
        j = randi(npop,1);
        k = randi(npop,1);
        if j==i || k==i
            j = randi(npop,1);
            k = randi(npop,1);
        end
        model_mut = [model(i,:);model(j,:)];
        mv_m =(model(i,:)+model(j,:))/2;
        bf = 1;
        for l = 1 : 2
            mod_mut(l,:) = model_mut(l,:) + rand*(model(k)-mv_m*bf);
            if mod_mut(l,1)<x0_min 
                mod_mut(l,1) = x0_min;
            end
            if mod_mut(l,2)<z1_min
                mod_mut(l,2) = z1_min;
            end
            if mod_mut(l,3)<z2_min
                mod_mut(l,3) = z2_min;
            end
            if mod_mut(l,4)<delta_Mz_min
                mod_mut(l,4) = delta_Mz_min;
            end
            if mod_mut(l,5)<m_min
                mod_mut(l,5) = m_min;
            end
            if mod_mut(l,1)>x0_max 
                mod_mut(l,1) = x0_max;
            end
            if mod_mut(l,2)>z1_max
                mod_mut(l,2) = z1_max;
            end
            if mod_mut(l,3)>z2_max
                mod_mut(l,3) = z2_max;
            end
            if mod_mut(l,4)>delta_Mz_max
                mod_mut(l,4) = delta_Mz_max;
            end
            if mod_mut(l,5)>m_max
                mod_mut(l,5) = m_max;
            end
        end
        %Hitung model untuk prosedur mutualisme
        for l = 1 : 2
            [dBdz_mut] = fwd_mgt(x,mod_mut(l,1),mod_mut(l,2),mod_mut(l,3),mu0,mod_mut(l,4),mod_mut(l,5));
            err_mut = misfit_magnetic(dBdz_obs,dBdz_mut);
            Em(l) = err_mut;
            %Update model jika  nilai misfit lebih baik proses mutualisme
            if l == 1
                if Em(l)<E(i)
                    model(i,:) = mod_mut(l,:);
                    E(i) = Em(l);
                    dBdz(i,:) = dBdz_mut;
                end
            else
                if Em(l)<E(j)
                    model(j,:) = mod_mut(l,:);
                    E(i) = Em(l);
                    dBdz(i,:) = dBdz_mut;
                end
            end
        end
        %Komensalisme
        j = randi(npop,1);
        if j == i
            j = randi(npop,1);
        end
        mod_com = model(i) +(0.4+0.9*rand)*(model_best-model(j));
            if mod_com(1)<x0_min 
                mod_com(1) = x0_min;
            end
            if mod_com(2)<z1_min
                mod_com(2) = z1_min;
            end
            if mod_com(3)<z2_min
                mod_com(3) = z2_min;
            end
            if mod_com(4)<delta_Mz_min
                mod_com(4) = delta_Mz_min;
            end
            if mod_com(5)<m_min
                mod_com(5) = m_min;
            end
            if mod_com(1)>x0_max 
                mod_com(1) = x0_max;
            end
            if mod_com(2)>z1_max
                mod_com(2) = z1_max;
            end
            if mod_com(3)>z2_max
                mod_com(3) = z2_max;
            end
            if mod_com(4)>delta_Mz_max
                mod_com(4) = delta_Mz_max;
            end
            if mod_com(5)>m_max
                mod_com(5) = m_max;
            end
        %Perhitungan misfit untuk prosedur komensalisme
        [dBdz_com] = fwd_mgt(x,mod_com(1),mod_com(2),mod_com(3),mu0,mod_com(4),mod_com(5));
         Ec = misfit_magnetic(dBdz_obs,dBdz_com);
         %Update model jika  nilai misfit lebih baik proses komensalisme
         if Ec < E(i)
             model(i,:) = mod_com(1,:);
             E(i) = Ec;
             dBdz(i,:) = dBdz_com(1,:);
         end
         %Parasitisme
         j = randi(npop,1);
         if j == i 
             j = randi(npop,1);
        end
         mod_par = model(i,:);
         p1 = randi(5,1);
         if p1 == 1
            mod_par(1) = x0_min + rand*(x0_max-x0_min);
         elseif p1 == 2
            mod_par(2) = z1_min + rand*(z1_max - z1_min);
         elseif p1 == 3
            mod_par(3) = z2_min + rand*(z2_max - z2_min);
         elseif p1 == 4
            mod_par(4) = delta_Mz_min + rand*(delta_Mz_max - delta_Mz_min);
         else
            mod_par(5) = m_min + rand*(m_max-m_min);
         end
         %Perhitungan misfit untuk tahap parasitisme
         [dBdz_par] = fwd_mgt(x,mod_par(1),mod_par(2),mod_par(3),mu0,mod_par(4),mod_par(5));
         Ep = misfit_magnetic(dBdz_obs,dBdz_par);
         %Update model jika  nilai misfit lebih baik proses parasitisme
         if Ep < E(j)
             model(j,:) = mod_par(1,:);
             E(j) = Ep;
             dBdz(j,:) = dBdz_par(1,:);
         end
    end
    %Update model terbaik untuk setiap iterasi
    Emin = 100;
    for ipop = 1 : npop
        Emin = E(ipop);
        model_baru = model(ipop,:);
        dBdz_model = dBdz(ipop,:);
    end
    %Nilai misfit terbaik
    Egen(itr)=Emin;

    %persiapan plot
x1 = x0 - (m);
x2 = x0 - (m);
x3 = x0 + (m);
x4 = x0 + (m);
y1 = z1;
y2 = z2;
y3 = z2;
y4 = z1;

x11 = model_baru(1) - model_baru(5);
x22 = model_baru(1) - model_baru(5);
x33 = model_baru(1) + model_baru(5);
x44 = model_baru(1) + model_baru(5);
y11 = model_baru(3);
y22 = model_baru(4);
y33 = model_baru(4);
y44 = model_baru(3);

figure(1)
subplot(2,1,1)
plot(x,dBdz_obs,'b.',x,dBdz_model,'r-','MarkerSize',20,'MarkerFaceColor','r','LineWidth',2.5)
xlabel('Distance (km)','FontWeight','bold')
ylabel(['\DeltaB_{z} (nT)'],'FontWeight','bold')
title(['\bf \fontsize{10}\fontname{Times}Geomagnetic model || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)],'FontWeight','bold');
subtitle(['x0 = ',num2str(model_baru(1)),' || z1 = ',num2str(model_baru(2)),' || z2 = ',num2str(model_baru(3)),' || \DeltaM_{z} = ',num2str(model_baru(4)),' || m = ',num2str(model_baru(5))],'FontWeight','bold')
legend({'Observed Data','Calculated Data'},'Location','Southeast')
grid on
subplot(2,1,2)
fill([x1 x2 x3 x4],[y1 y2 y3 y4],'y')
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
set(gca,'Color','b')
end

figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RSME','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Misfit ');
set(gcf, 'Position', get(0, 'Screensize'));
grid on