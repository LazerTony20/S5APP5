% APP5 S5 Laboratoire 1 
% Probleme 6
opengl software
close all
clear
clc

%Valeurs de depart
Mp = 6;              %En pourcentage
tr10_90 = 0.004;     %En secondes
tp = 0.008;          %En secondes
ts2pc = 0.010;       %En secondes
Erpr = 0.00005;      %Erreur Rampe plus petit que
numGs = [4500];      %Num de ma FTBO G(s)
denGs = [1 361.2 0]; %Den de ma FTBO G(s)
Gs = tf(numGs,denGs);%FTBO
% On peut rencontrer ces spécifications avec un AvPh suivi d’un RePh ou d’un PD suivi d’un PI (i.e. un PID).

%% A) Traduire les 4 premières spécifications en termes de pôles désirés.
Phi_A = atand(-pi./log((Mp./100)));
Zeta = cosd(Phi_A);  % Zeta = 0.6671
wn_tp = pi./(tp.*sqrt(1-Zeta.^2));
wn_ts2pc = 4./(Zeta.*ts2pc);
wn_tr_10_90 = (1 + 1.1.*Zeta + 1.4.*(Zeta.^2))./tr10_90;
PM_star = atand(2*Zeta/sqrt(sqrt(1+4*Zeta^4)-2*Zeta^2));
wg_star = 2*Zeta*wn_ts2pc/tand(PM_star);
disp(['Wn de tp = ',num2str(wn_tp)])
disp(['Wn de ts2% = ',num2str(wn_ts2pc)])
disp(['Wn de tr(10-90%) = ',num2str(wn_tr_10_90)])
disp(' ')
disp(['La valeur la plus haute (avoir le plus petit ts) est Wn de ts2%, soit : ',num2str(wn_ts2pc)])
wa_ts2pc = wn_ts2pc.*sqrt(1-Zeta.^2);
pstar = [(-Zeta.*wn_ts2pc + i.*(wa_ts2pc)); (-Zeta.*wn_ts2pc - i.*(wa_ts2pc))];
disp(['Pole desire = ',num2str(pstar(1))])
% Affichage des poles desires
figure('Name','Poles desires')
hold on
rlocus(Gs)
plot(real(pstar), imag(pstar),'p')
axis([-500 100 -500 500])
title('Poles desires')
hold off
%% Conception AvPh + RePh avec lieu des racines
%% B) Concevoir un compensateur AvPh

s_b = pstar(1);
Gspstar = polyval(numGs,s_b)./polyval(denGs,s_b);
mag_gspstar = abs(Gspstar);
pha_gspstar = angle(Gspstar).*180/pi -360;
disp(' ')
disp(['La valeur de la fonction à ', num2str(s_b), ' est: ', num2str(Gspstar)])
disp(['Le module de ', num2str(Gspstar), ' est = ', num2str(mag_gspstar)])
disp(['La phase de ', num2str(Gspstar), ' est = ', num2str(pha_gspstar), ' degrés'])

alpha_b = 180 - Phi_A;
DeltaPhi = -180 - pha_gspstar;
PhiZ = (alpha_b + DeltaPhi)./2;
PhiP = (alpha_b - DeltaPhi)./2;
Za = real(s_b)- imag(s_b)./tand(PhiZ);
Pa = real(s_b)- imag(s_b)./tand(PhiP);
disp(' ')
disp(['Alpha est = ', num2str(alpha_b), ' degrés'])
disp(['DeltaPhi est = ', num2str(DeltaPhi), ' degrés'])
disp(['Phi_z est = ', num2str(PhiZ), ' degrés'])
disp(['Phi_p est = ', num2str(PhiP), ' degrés'])
disp(['Za est = ', num2str(Za)])
disp(['Pa est = ', num2str(Pa)])
disp('')
disp(['Verif avec PhiZ - PhiP ', num2str((PhiZ-PhiP)), ' degrés'])
num_Gsa = [1 -Za];
den_Gsa = [1 -Pa];
Ka = abs( (polyval(denGs,s_b)*polyval(den_Gsa,s_b))./(polyval(num_Gsa,s_b)*polyval(numGs,s_b)));
Gsa = tf(Ka*num_Gsa,den_Gsa);
disp(['Compensateur Avance de phase: Ka = ', num2str(Ka), ' zéro = ',num2str(Za), ' pole = ', num2str(Pa)]);
disp(' ')
GsaGs = Gs*Gsa;
[numa, dena] = tfdata(GsaGs,'v');

figure('Name','FT originale, FT avec AvPh')
hold on
rlocus(Gs,'r')
rlocus(GsaGs, 'b')
plot(pstar,'p')
pos = rlocus(GsaGs,1);
plot(pos,'s')
title('FT originale, FT avec AvPh, pôles désirés * et pôles placés')
% 
% 
% rlocus(Gs)
% plot(real(pstar), imag(pstar),'p')
% rlocus()
% axis([-500 100 -500 500])
% title('Poles desires')
% hold off
%% C) En partant de l’AvPh en (b), faire la conception cette fois d’un RePh
Kvel_now = numa(end)/dena(end-1);
Kvel_des = 1./Erpr;
Kvel_fac_RePh = Kvel_des/Kvel_now;
F = 10;
disp(['Le Kvel actuel est de ', num2str(Kvel_now)])
disp(['Le Kvel désiré est de ', num2str(Kvel_des)])
disp(['L`augmentation de gain requise pour atteindre l`erreur en RP est de 𝐾^* = ', num2str(Kvel_fac_RePh)])
disp(' ')

[MAG_RePh, PHA_RePh] = bode(Kvel_fac_RePh.*GsaGs, wg_star);
beta_RePh = MAG_RePh

% Calcul du compensateur:
% On place le zéro a une décade avant wg_des
s_zer = real(s_b)/F;
% On place le pole a un facteur beta plus à gauche
s_pol = s_zer/Kvel_fac_RePh;
 
% L'amplification aux basses frequences est K = Kr*beta = Kvel_fac
% Donc, Kr = Kvel_fac/beta

[mag, pha] = bode([1 -s_zer],[1 -s_pol], wg_star);
% Kr = Kvel_fac/(beta*mag)


Kr = Kvel_fac_RePh/(MAG_RePh);

disp(['Le zéro est placé à 𝑟𝑒𝑎𝑙(𝑝∗)/𝐹 = ', num2str(real(s_b)/F), ' avec facteur 𝐹 = ', num2str(F)])
disp(['Le pôle est placé à zéro/𝐾∗ = ', num2str(s_pol)])
disp(['Le gain exact du RePh est Kr = ', num2str(Kr), ' mais on le laisse à ', num2str(1.0)])
disp(' ')

Kr = 1;
numr = [1 -s_zer]*Kr;
denr = [1 -s_zer/Kvel_fac_RePh];

disp('Compensateur RePh')
RePh = tf(numr,denr)

% Calcule la FTBO avec le compensateur 
GsAvPhRePh = GsaGs*RePh;
[numar, denar] = tfdata(GsAvPhRePh,'v');


figure('Name','Retard de phase')
hold on
rlocus(GsaGs,'b')
rlocus(GsAvPhRePh, 'r')
plot(pstar,'p')
pos = rlocus(GsaGs,1);
plot(pos,'s')
pos2 = rlocus(GsAvPhRePh,1);
plot(pos2,'bo')
legend('FT avec AvPh (bissectrice)','FT avec AvPh (bissectrice) et RePh','Pôles désirés','Pôles placés avec AvPh','Pôles avec AvPh (bissectrice) et RePh','Location','SouthWest')
hold off


%% Refaire (b) et (c) en prévoyant une surcompensation de 3 deg dans la conception de l’AvPh. 
% Vérifier que les pôles obtenus (les cercles) sont maintenant plus près des pôles désirés (les pentagrammes)




%% Conception PD + PI avec lieu des racines
%% F) À la place de l’AvPh en (b) faire un PD cette fois


%% G) À la place du RePh en (c) faire un PI en cascade au PD ci-dessus (créant ainsi un PID) 
% qui rencontre les mêmes spécifications en régime transitoire et permanent. 
% Commenter les différences. À noter qu’en pratique, on ne fait jamais un
% PI quand les spécifications ne requièrent pas un changement de classe.








