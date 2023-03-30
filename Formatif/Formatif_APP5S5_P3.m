% APP5 S5 Formatif 
% Probleme 3
clc
close all
clear
clc

% Un système asservi à retour unitaire a pour fonction de transfert en boucle ouverte :
numGs = [1];
denGs = conv([1 0],conv([1 5],[1 11]));
FT = tf(numGs,denGs)

%% (a) Un asservissement de type proportionnel de gain 𝐾𝑝 = 218.7 a initialement été conçu. 
% Calculer les performances de cet asservissement : 
% - pôles en boucle fermée, 
% - dépassement maximum 𝑀𝑝,
% - temps de stabilisation à 2% 𝑡𝑠(mode dominant) 
% - erreur en régime permanent pour une entrée rampe unitaire.
Kp = 218.7;
p_FB = rlocus(FT,Kp);
p_des = [p_FB(2); p_FB(3)];
s_des = p_des(1);
wn = abs(s_des);
Zeta = -real(s_des)./wn;
Phi = acosd(Zeta);
Mp = 100*exp(-pi./tand(Phi));
ts2pc = 4./(Zeta.*wn);
tp = pi./imag(s_des);
Kvel = Kp.*(numGs(end))./(denGs(end-1));
Erp = 1./Kvel;
disp(['s* = ', num2str(real(s_des)), ' ± ', num2str(imag(s_des)), 'i (il y a un autre pôle à ', num2str(real(p_FB(1))),')'])
disp(['𝑀𝑝 = ', num2str(Mp), '%'])
disp(['𝑡𝑠 = ', num2str(ts2pc), ' 𝑠'])
disp(['𝑡𝑝 = ', num2str(tp), ' 𝑠'])
disp(' ')
disp(['𝐾𝑣𝑒𝑙 = ', num2str(Kvel)])
disp(['𝑒𝑅𝑃 = ', num2str(Erp)])
disp(' ')

figure('Name','Lieu des racines de la FT originale et pôles désirés')
hold on
rlocus(FT,'r')
plot(real(p_FB),imag(p_FB),'p','MarkerEdgeColor','b')
axis([-15 5 -15 15])
legend('FT Originale','Pôles désirés','Location','NorthWest')
hold off
%% (b) La performance obtenue avec ce compensateur P n’est pas acceptable. 
% Concevoir un compensateur AvPh et un compensateur RePh en cascade avec 𝐺(𝑠) pour réduire 
% le temps de stabilisation et le dépassement maximum d’un facteur 2 et réduire l’erreur en régime permanent d’un facteur 10.
% Utiliser un facteur 𝐹 = 10 pour le placement du zéro du RePh.
disp('===================(B)===================')
F = 10;
Facteur = 2;
ts2pc2 = ts2pc./2;
Mp2 = Mp./2;
phi2 = atand(-pi/log(Mp2/100));
Zeta2 = cosd(phi2);
wn2_des = 4.0/(ts2pc2*Zeta2); % fréquence des oscillations
wa2 = wn2_des.*sqrt(1-Zeta2.^2);
Erp2 = Erp./F;
p_des_2 = [(-(Zeta2.*wn2_des) + i*wa2); (-(Zeta2.*wn2_des) - i*wa2)];
s_des_2 = p_des_2(1);
disp(['𝑡𝑠* = ', num2str(ts2pc2), ' 𝑠']);
disp(['𝑀𝑝* = ', num2str(Mp2), ' %'])
disp(['𝑒𝑅𝑃* = ', num2str(Erp2)])
disp(['𝜙* = ', num2str(phi2), ' deg'])
disp(['𝜁* = cos𝜙* = ', num2str(Zeta2)])
disp(['𝜔𝑛* = ', num2str(wn2_des)])
disp(['Les pôles désirés sont à 𝑠* = ', num2str(real(s_des_2)), ' ± ', num2str(imag(s_des_2)), 'i'])
disp('======Compensateur AvPh======')
alpha = 180 - phi2;
delta_phi = -180 - rad2deg(angle(polyval(numGs,s_des_2)./polyval(denGs,s_des_2))) + 360;
Phi_Z = (alpha + delta_phi)./2;
Phi_P = (alpha - delta_phi)./2;
z_pos = real(s_des_2) - (imag(s_des_2)./tand(Phi_Z));
p_pos = real(s_des_2) - (imag(s_des_2)./tand(Phi_P));
Ka_AvPh = 1./abs(((s_des_2 - z_pos)./(s_des_2 - p_pos)).*(polyval(numGs,s_des_2)./polyval(denGs,s_des_2)));

disp(['Avance de phase requise selon les spécifications = ', num2str(delta_phi), ' deg'])
disp(' ')
disp(['La position du zero de l`AvPh est ', num2str(z_pos)])
disp(['La position du pôle de l`AvPh est ', num2str(p_pos)])
disp(['Le gain de l`AvPh est ', num2str(Ka_AvPh)])
disp(' ')

num_AvPh = Ka_AvPh .*[1 -z_pos];
den_AvPh = [1 -p_pos];
disp('Compensateur Avance de phase:')
AvPh = tf(num_AvPh,den_AvPh)
disp(' ')
AvPh_Gs = AvPh * FT;
[numa,dena] = tfdata(AvPh_Gs,'v');

figure('Name','FT originale et FT avec AvPh simple')
hold on
rlocus(FT,'r')
rlocus(AvPh_Gs,'b')
p_e = rlocus(AvPh_Gs,1);
plot(real(p_des_2), imag(p_des_2),'p','MarkerEdgeColor','b')
plot(real(p_e), imag(p_e),'s','MarkerEdgeColor','g')
axis([-16 2 -8 8])
legend('FT originale', 'FT avec AvPh simple', 'pôles désirés','pôles obtenus','Location','NorthWest')
title('FT originale et FT avec AvPh simple')
hold off
disp('=====Compensateur Retard de phase====')
Erp_RePh = Erp2;
Kvel_des = 1./Erp_RePh;
Kvel_now = numa(end)./dena(end-1);
Erp_now = 1./Kvel_now;

K_des = Kvel_des./Kvel_now;
disp(['eRP_des = ', num2str(Erp_RePh)])
disp(['Kvel_des = ', num2str(Kvel_des)])
disp([' '])
disp(['Kvel_now = ', num2str(Kvel_now)])
disp(['eRP_now = ', num2str(Erp_now)])
disp([' '])
disp(['K_des = ', num2str(K_des)])
disp([' '])

zr = real(p_des_2(1))/F;
pr = zr/K_des;

numr = [1 -zr];
denr = [1 -pr];
disp('Compensateur RePh')
FTr = tf(numr,denr)
AvPh_FT_RePh = AvPh_Gs*FTr;

figure('Name','FT originale,avec AvPh simple et RePh (F=10)')
hold on
rlocus(FT,'r')
rlocus(AvPh_Gs,'b')
rlocus(AvPh_FT_RePh,'g')
p_e = rlocus(AvPh_FT_RePh,1);
plot(real(p_des_2), imag(p_des_2),'p','MarkerEdgeColor','b')
plot(real(p_e), imag(p_e),'s','MarkerEdgeColor','g')
axis([-16 2 -8 8])
legend('FT originale', 'FT avec AvPh simple', 'FT avec AvPh et RePh simples', 'pôles désirés','pôles obtenus','Location','NorthWest')
title('FT originale,avec AvPh simple et RePh (F=10)')
hold off

%% (c) Refaire le RePh mais avec cette fois un facteur 𝐹 = 2 pour placer son zéro. 
% Comparer le graphique du lieu des racines et de l’erreur à la rampe pour 𝐹 = 10 et 𝐹 = 2 
% et observer la dégradation dans un cas et l’amélioration dans l’autre.
F_ramp2 = 2;

zr2 = real(p_des_2(1))/F_ramp2;
pr2 = zr2/K_des;

numr2 = [1 -zr2];
denr2 = [1 -pr2];
disp('Compensateur RePh (F=2)')
FTr2 = tf(numr2,denr2)
AvPh_FT_RePh2 = AvPh_Gs*FTr2;

figure('Name','FT originale,avec AvPh simple et RePh (F=2)')
hold on
rlocus(FT,'r')
rlocus(AvPh_Gs,'b')
rlocus(AvPh_FT_RePh2,'g')
p_e = rlocus(AvPh_FT_RePh2,1);
plot(real(p_des_2), imag(p_des_2),'p','MarkerEdgeColor','b')
plot(real(p_e), imag(p_e),'s','MarkerEdgeColor','g')
axis([-16 2 -8 8])
legend('FT originale', 'FT avec AvPh simple', 'FT avec AvPh et RePh simples', 'pôles désirés','pôles obtenus','Location','NorthWest')
title('FT originale,avec AvPh simple et RePh (F=2)')
hold off

