% APP5 S5 Formatif 
% Probleme 2
clc
close all
clear
clc


% Le système en boucle ouverte ci-dessous est asservi avec une rétroaction unitaire et un compensateur en cascade.
numGs = [1];
denGs = conv(conv([1 2],[1 2]),[1 3]);
disp('La fonction de transfert est :')
Gs = tf(numGs, denGs)
%% (a) Calculer la position des pôles dominants pour avoir temps de stabilisation à 2% de 1.6 seconde et un dépassement maximum de 25%.
ts2pc = 1.6;    % En secondes
Mp = 25;        % En Pourcent
% Conclusion, domaine temporel
Phi_A = atand(-pi./log(Mp./100));
Zeta_A = cosd(Phi_A);
wn_ts2pc = 4./(ts2pc*Zeta_A);
wa_A = wn_ts2pc.*sqrt(1 - Zeta_A.^2)
p_des = [(-(Zeta_A.*wn_ts2pc) + i*wa_A); (-(Zeta_A.*wn_ts2pc) - i*wa_A)]
s_des = p_des(1);
disp(['Φ = ', num2str(Phi_A), ' deg'])
disp(['ζ = ', num2str(Zeta_A)])
disp(['ωn = ', num2str(wn_ts2pc), ' rad/s'])
disp(' ')
disp('Les pôles désirés sont : ')
disp(num2str(p_des(1)))
disp(num2str(p_des(2)))
disp(' ')

figure('Name','Lieu des racines et pôles désirés')
hold on
rlocus(numGs,denGs)
plot(real(p_des), imag(p_des),'p')
title('Lieu des racines et pôles désirés')
axis([-6 1 -7 7])
hold off
%% (b) *** Si un compensateur AvPh avec un zéro à −1 est utilisé pour satisfaire les conditions (a), quelle devrait être la contribution angulaire du pôle du compensateur.
z_pos = -1;
Delta_Phi = (-180) - rad2deg(angle(polyval(numGs,s_des)./polyval(denGs,s_des))) + 360;
Phi_Z = 180 + atand((imag(s_des))./(real(s_des)-z_pos));
Phi_P = Phi_Z - Delta_Phi;
disp(['DeltaΦ = ', num2str(Delta_Phi), ' deg'])
disp(['ΦZ = ', num2str(Phi_Z), ' deg'])
disp(['ΦP = ', num2str(Phi_P), ' deg'])
disp(' ')
%% (c) Calculer la position du pôle du compensateur.
p_pos = real(s_des)-(imag(s_des)./tand(Phi_P));
disp(['La position du pôle de l`AvPh est ', num2str(p_pos)])
disp(' ')
%% (d) Calculer la valeur du gain du compensateur et donner la FT du compensateur.
Ka_AvPh = 1./abs(((s_des - z_pos)./(s_des - p_pos)).*(polyval(numGs,s_des)./polyval(denGs,s_des)));
disp(['Le gain de l`AvPh est ', num2str(Ka_AvPh)])
disp(' ')
num_AvPh = Ka_AvPh .*[1 -z_pos];
den_AvPh = [1 -p_pos];
disp('Compensateur Avance de phase:')
AvPh = tf(num_AvPh,den_AvPh)
disp(' ')
AvPh_Gs = AvPh * Gs;
%% (e) Avec rlocus, vérifier la position obtenue des pôles en boucle fermée du système compensé.
figure('Name','FT originale et FT avec AvPh simple')
hold on
rlocus(Gs,'r')
rlocus(AvPh_Gs,'b')
p_e = rlocus(AvPh_Gs,1);
plot(real(p_des), imag(p_des),'p')
plot(real(p_e), imag(p_e),'s')
axis([-7 2 -10 10])
legend('FT originale', 'FT avec AvPh simple', 'pôles désirés','pôles obtenus','Location','NorthWest')
title('FT originale et FT avec AvPh simple')
hold off
%% (f) Avec une simulation en boucle fermée, vérifier si les spécifications sont rencontrées. Est-ce que l’approximation de second ordre est valide?
AvPh_Gs_FB = feedback(AvPh_Gs,1);
t = [0:0.01:5];
u = ones(size(t));
[num_AvPh_Gs_FB, den_AvPh_Gs_FB] = tfdata(AvPh_Gs_FB,'v');
simu = lsim(num_AvPh_Gs_FB,den_AvPh_Gs_FB,u,t);
simuinfos = lsiminfo(simu,t);
Mp_simu = (max(simu)-simu(end))/simu(end);
figure('Name','Réponse à un échelon unitaire')
hold on
plot(t,simu,'b')
yline(0.98.*simu(end),'r--')
yline(1.02.*simu(end),'r--')
ylabel('Réponse')
title('Réponse à un échelon unitaire')
grid on
hold off
disp(['ts(2%) = ', num2str(simuinfos.SettlingTime), ' s', ' (vs ', num2str(ts2pc), ' s)'])
disp(['Mp = ' num2str(Mp_simu.*100), ' %', ' (vs ', num2str(Mp), ' %)'])
disp('Les spécifications ne sont pas rencontrées (mais pas trop loin).')
disp(' ')
%% (g) *** On remarque que l’avance de phase requise est plus élevée que 75 deg et le zéro de l’AvPh est à droite des pôles désirés, conditions qui ne sont pas acceptables. Refaire le problème avec un double  AvPh avec les deux zéros à −3. 
z_pos_AvPhdbl = -3;
Phi_Z_2AvPh = atand((imag(s_des))./(real(s_des)-z_pos_AvPhdbl));
Phi_P_2AvPh = Phi_Z_2AvPh - Delta_Phi./2;
p_pos_2AvPh = real(s_des)-(imag(s_des)./tand(Phi_P_2AvPh));
Ka_2AvPh = 1./abs((((s_des - z_pos_AvPhdbl)./(s_des - p_pos_2AvPh)).^2).*(polyval(numGs,s_des)./polyval(denGs,s_des)));
disp(['La contribution angulaire de chaque zéro à ', num2str(z_pos_AvPhdbl), ' est ', num2str(Phi_Z_2AvPh), ' deg (< 90 deg donc à gauche)'])
disp(['La contribution angulaire de chaque pôle est donc ', num2str(Phi_P_2AvPh) ' deg'])
disp(' ')
disp(['La position du zéro de chaque AvPh est ', num2str(z_pos_AvPhdbl)])
disp(['La position du pôle de chaque AvPh est ', num2str(p_pos_2AvPh)])
disp(['Le gain de l`AvPh est ', num2str(Ka_2AvPh)])
disp('FT pour double AvPh:')
num_2AvPh = conv([1 -z_pos_AvPhdbl],[1 -z_pos_AvPhdbl]);
den_2AvPh = conv([1 -p_pos_2AvPh],[1 -p_pos_2AvPh]);
AvPh2 = Ka_2AvPh.*tf(num_2AvPh,den_2AvPh)
AvPh2_Gs = AvPh2*Gs;
disp(' ')
figure('Name','FT originale et FT avec AvPh double')
hold on
rlocus(Gs,'r')
rlocus(AvPh2_Gs,'b')
p_e = rlocus(AvPh2_Gs,1);
plot(real(p_des), imag(p_des),'p')
plot(real(p_e), imag(p_e),'s')
axis([-10 4 -10 10])
legend('FT originale', 'FT avec AvPh double', 'pôles désirés','pôles obtenus','Location','NorthWest')
title('FT originale et FT avec AvPh double')
hold off
%% (h) Refaire (g) mais en ajoutant en cascade à 𝐺(𝑠) un capteur qui ajoute un retard pur de 0.05 s. Utiliser une approximation Pade d’ordre 7.
retard = 0.05;  % En secondes
ordrePade = 7;
[numEst, denEst] = pade(retard,ordrePade);
Est = tf(numEst, denEst)
GsEst = Gs * Est;
p_des_H = rlocus(GsEst,1)









