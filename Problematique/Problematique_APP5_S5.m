% APP5 S5 Problematique
% ROYA2019
opengl software
close all
clear
clc

show_intermediary_graphs = 1;

%% Valeurs de depart
% Modèle dynamique des deux télescopes en azimut (AZ)
numAZ = [1.59e09];
denAZ = [1 1020.51 25082.705 3102480.725 64155612.5 82700000 0];
% Modèle dynamique des deux télescopes en élévation (EL)
numEL = [7.95e09];
denEL = [1 1020.51 37082.705 15346520.725 320776412.5 413500000 0];
FT_AZ = tf(numAZ,denAZ);
FT_EL = tf(numEL,denEL);

util_cb = 0; %Utilisation du Coupe-Bande dans les calculs (VESTIGE NON UTILISÉ)



% figure('Name','Lieu des racines de FTAZ originale')
% hold on
% rlocus(FT_AZ,'b')
% hold off
% 
if show_intermediary_graphs == 1
    figure('Name','FT_AZ de Base')
    margin(FT_AZ)
    grid on
    figure('Name','FT_EL de Base')
    margin(FT_EL)
    grid on
end
%% =====Coupe bande======
% https://en.wikipedia.org/wiki/Band-stop_filter 
w0_cb_AZ = 54.8;    %Depuis bode
wc_cb_AZ = 20;
numBut_AZ = [1 1 w0_cb_AZ.^2];
denBut_AZ = [1 wc_cb_AZ w0_cb_AZ.^2];
Butt_AZ = tf(numBut_AZ,denBut_AZ);
%
w0_cb_EL = 122.5;    %Depuis bode
wc_cb_EL = 20;
numBut_EL = [1 1 w0_cb_EL.^2];
denBut_EL = [1 wc_cb_EL w0_cb_EL.^2];
Butt_EL = tf(numBut_EL,denBut_EL);

if util_cb == 1
    FT_AZ = tf(numAZ,denAZ)*Butt_AZ;
    FT_EL = tf(numEL,denEL)*Butt_EL;
    [numAZ,denAZ] = tfdata(FT_AZ,'v');
    [numEL,denEL] = tfdata(FT_EL,'v');
    figure(1)
    hold on
    margin(FT_AZ)
    legend('FTAZ','FTAZ (Avec CB)')
    grid on
    hold off
    figure(2)
    hold on
    margin(FT_EL)
    legend('FTEL','FTEL (Avec CB)')
    grid on
    hold off
end



% =======================================================
% ======================Telescope A======================
% =======================================================
disp(' ')
disp('=======================================================')
disp('======================Telescope A======================')
disp('=======================================================')
disp(' ')
% Reponse a un echelon unitaire
Erp_Ech_Ta = 0;
% Reponse a une rampe unitaire
Erp_Ramp_Ta_Az = 0.03;   %En degres
Erp_Ramp_Ta_El = 0;      %En degres
% Reponse a une parabole unitaire
Erp_Par_El = 0.08;    %En degres

% Traduction des spécifications
% En azimut et elevation
disp('Traduction des spécifications')
Mp_Ta = 25;          %En pourcentage
ts2pc_Ta = 1;        %En secondes
ts1090_Ta = 0.18;    %En secondes
phi_Ta = atand(-pi./(log(Mp_Ta./100)));
Zeta_Ta = cosd(phi_Ta); 
wn_ts2pc_Ta = 4./(Zeta_Ta.*ts2pc_Ta);
wn_ts1090_Ta = (1 + 1.1.*Zeta_Ta + 1.4.*(Zeta_Ta.^2))./ts1090_Ta;
disp(['wn_ts2pc_Ta = ', num2str(wn_ts2pc_Ta),' <---'])
disp(['wn_ts1090_Ta = ', num2str(wn_ts1090_Ta)])
disp(['Zeta_Ta = ', num2str(Zeta_Ta)])

% Calcul des pôles désirés selon les spécifications
poles_des_Ta = [-Zeta_Ta*wn_ts2pc_Ta+1i*wn_ts2pc_Ta*sqrt(1-Zeta_Ta^2) -Zeta_Ta*wn_ts2pc_Ta-1i*wn_ts2pc_Ta*sqrt(1-Zeta_Ta^2)];
disp(['Les pôles désires sont à '])
disp(num2str(poles_des_Ta(1)))
disp(num2str(poles_des_Ta(2)))
disp(' ')
% Calcul de la frequence de la traverse de gain
wg_des_Ta = wn_ts2pc_Ta.*sqrt(sqrt(1+4.*Zeta_Ta.^4)-2.*Zeta_Ta.^2);

% Concevoir un compensateur AvPh
disp('=========Concevoir un compensateur AvPh=========')
s_pchoix_Ta = poles_des_Ta(1);
alpha_Ta = 180 - phi_Ta;
% FTAZ
FT_AZ_des = polyval(numAZ,s_pchoix_Ta)./polyval(denAZ,s_pchoix_Ta);
Mag_des_Ta_AZ = abs(FT_AZ_des);
Pha_des_Ta_Az = angle(FT_AZ_des).*180/pi -360;
deltaPhi_des_Ta_Az = -180 - Pha_des_Ta_Az;

% FTEL
FT_EL_des = polyval(numEL,s_pchoix_Ta)./polyval(denEL,s_pchoix_Ta);
Mag_des_Ta_EL = abs(FT_EL_des);
Pha_des_Ta_EL = angle(FT_EL_des).*180/pi -360;
deltaPhi_des_Ta_EL = -180 - Pha_des_Ta_EL;

disp(' ')
disp(['La valeur de la fonction FTAZ_des à ', num2str(s_pchoix_Ta), ' est: ', num2str(FT_AZ_des)])
disp(['Le module de ', num2str(FT_AZ_des), ' est = ', num2str(Mag_des_Ta_AZ)])
disp(['La phase de ', num2str(FT_AZ_des), ' est = ', num2str(Pha_des_Ta_Az), ' degrés'])
disp(['deltaPhi_des_Ta_Az ', num2str(deltaPhi_des_Ta_Az), ' deg'])
disp(' ')
disp(['La valeur de la fonction FTEL_des à ', num2str(s_pchoix_Ta), ' est: ', num2str(FT_EL_des)])
disp(['Le module de ', num2str(FT_EL_des), ' est = ', num2str(Mag_des_Ta_EL)])
disp(['La phase de ', num2str(FT_EL_des), ' est = ', num2str(Pha_des_Ta_EL), ' degrés'])
disp(['deltaPhi_des_Ta_EL ', num2str(deltaPhi_des_Ta_EL), ' deg'])
disp(' ')

% Calcul du gain K* (Kdes) pour obtenir la fréquence de traverse à wg* (wg_des)
% C'est l'inverse de l'amplitude de la FTBO à wg*
K_des_Ta_FT_AZ = 1./Mag_des_Ta_AZ;
disp(['K_des_Ta_FT_AZ = ', num2str(K_des_Ta_FT_AZ), ' (K*)'])

K_des_Ta_FT_EL = 1./Mag_des_Ta_EL;
disp(['K_des_Ta_FT_EL = ', num2str(K_des_Ta_FT_EL), ' (K*)'])

PhiZ_AZ_des = (alpha_Ta + deltaPhi_des_Ta_Az)./2;
PhiP_AZ_des = (alpha_Ta - deltaPhi_des_Ta_Az)./2;
Za_AZ_des = real(s_pchoix_Ta)- imag(s_pchoix_Ta)./tand(PhiZ_AZ_des);
Pa_AZ_des = real(s_pchoix_Ta)- imag(s_pchoix_Ta)./tand(PhiP_AZ_des);

PhiZ_EL_des = (alpha_Ta + deltaPhi_des_Ta_EL)./2;
PhiP_EL_des = (alpha_Ta - deltaPhi_des_Ta_EL)./2;
Za_EL_des = real(s_pchoix_Ta)- imag(s_pchoix_Ta)./tand(PhiZ_EL_des);
Pa_EL_des = real(s_pchoix_Ta)- imag(s_pchoix_Ta)./tand(PhiP_EL_des);

disp(' ')
disp(['Alpha est = ', num2str(alpha_Ta), ' degrés'])
disp(' ')
disp(['DeltaPhi_AZ est = ', num2str(deltaPhi_des_Ta_Az), ' degrés'])
disp(['Phi_z_AZ est = ', num2str(PhiZ_AZ_des), ' degrés'])
disp(['Phi_p_AZ est = ', num2str(PhiP_AZ_des), ' degrés'])
disp(['Za_AZ est = ', num2str(Za_AZ_des)])
disp(['Pa_AZ est = ', num2str(Pa_AZ_des)])
disp(['Verif avec PhiZ - PhiP ', num2str((PhiZ_AZ_des-PhiP_AZ_des)), ' degrés'])
disp(' ')
disp(['DeltaPhi_EL est = ', num2str(deltaPhi_des_Ta_EL), ' degrés'])
disp(['Phi_z_EL est = ', num2str(PhiZ_EL_des), ' degrés'])
disp(['Phi_p_EL est = ', num2str(PhiP_EL_des), ' degrés'])
disp(['Za_EL est = ', num2str(Za_EL_des)])
disp(['Pa_EL est = ', num2str(Pa_EL_des)])
disp(['Verif avec PhiZ - PhiP ', num2str((PhiZ_EL_des-PhiP_EL_des)), ' degrés'])
disp(' ')

num_Gsa_AZ = [1 -Za_AZ_des];
den_Gsa_AZ = [1 -Pa_AZ_des];
Ka_AZ = abs( (polyval(denAZ,s_pchoix_Ta)*polyval(den_Gsa_AZ,s_pchoix_Ta))./(polyval(num_Gsa_AZ,s_pchoix_Ta)*polyval(numAZ,s_pchoix_Ta)));
disp('Fonction de transfer de l`AvPh (AZ)')
AvPh_AZ_Ta = tf(Ka_AZ*num_Gsa_AZ,den_Gsa_AZ)
disp(['Compensateur Avance de phase (AZ):']);
disp(['Ka = ', num2str(Ka_AZ)]);
disp([' zéro = ',num2str(Za_AZ_des)])
disp([' pole = ', num2str(Pa_AZ_des)])
AvPh_FT_AZ_Ta = AvPh_AZ_Ta*FT_AZ;
[numa_AZ_Ta, dena_AZ_Ta] = tfdata(AvPh_FT_AZ_Ta,'v');
disp(' ')
num_Gsa_EL = [1 -Za_EL_des];
den_Gsa_EL = [1 -Pa_EL_des];
Ka_EL = abs( (polyval(denEL,s_pchoix_Ta)*polyval(den_Gsa_EL,s_pchoix_Ta))./(polyval(num_Gsa_EL,s_pchoix_Ta)*polyval(numEL,s_pchoix_Ta)));
disp('Fonction de transfer de l`AvPh (EL)')
AvPh_EL_Ta = tf(Ka_EL*num_Gsa_EL,den_Gsa_EL)
disp(['Compensateur Avance de phase (EL):']);
disp(['Ka = ', num2str(Ka_EL)])
disp([' zéro = ',num2str(Za_EL_des)])
disp([' pole = ', num2str(Pa_EL_des)])
AvPh_FT_EL_Ta = AvPh_EL_Ta*FT_EL;
[numa_EL_Ta, dena_EL_Ta] = tfdata(AvPh_FT_EL_Ta,'v');
disp(' ')

% Comparaison Avant-Apres
% AZ
if show_intermediary_graphs == 1
    figure('Name','Lieu de bode de AZ avant et avec AvPh (Ta)')
    hold on
    margin(FT_AZ)
    margin(AvPh_FT_AZ_Ta)
    legend('FTAZ','FTAZ + AvPh','Location','NorthEast')
    grid on
    hold off
    % EL
    figure('Name','Lieu de bode de EL avant et avec AvPh (Ta)')
    hold on
    margin(FT_EL)
    margin(AvPh_FT_EL_Ta)
    legend('FTEL','FTEL + AvPh','Location','NorthEast')
    grid on
    hold off
end


% Concevoir un compensateur RePh
disp('=========Concevoir un compensateur RePh=========')
disp(' ')
Kvel_now_AZ_Ta = numa_AZ_Ta(end)/dena_AZ_Ta(end-1);

% Az
Kvel_des_AZ_Ta = 1./Erp_Ramp_Ta_Az;
K_des_AZ_RePh_Ta = Kvel_des_AZ_Ta./Kvel_now_AZ_Ta;
F_AZ = 10;

Zr_AZ_Ta = real(s_pchoix_Ta)./F_AZ;
Pr_AZ_Ta = Zr_AZ_Ta./K_des_AZ_RePh_Ta;
Kr_AZ_RePh_Ta = 1;

numRe_AZ_Ta = [1 -Zr_AZ_Ta];
denRe_AZ_Ta = [1 -Pr_AZ_Ta];
disp('Fonction de transfer du RePh (AZ)')
RePh_Ta = tf((Kr_AZ_RePh_Ta.*numRe_AZ_Ta),denRe_AZ_Ta)

% disp(['Fonction de transfert Azimut (Ta) avec AvPh et RePh:'])
AvPh_FTAZ_RePh_Ta = AvPh_FT_AZ_Ta * RePh_Ta;

% Comparaison Avant-Apres
% AZ
if show_intermediary_graphs == 1
    figure('Name','Lieu de bode de AZ original,avec AvPh et avec RePh (Ta)')
    hold on
    margin(FT_AZ)
    margin(AvPh_FT_AZ_Ta)
    margin(AvPh_FTAZ_RePh_Ta)
    legend('FTAZ','FTAZ + AvPh','FTAZ + AvPh + RePh','Location','NorthEast')
    grid on
    hold off 
end

% Concevoir un compensateur PI
disp('=========Concevoir un compensateur PI=========')
disp(' ')
Kvel_now_EL_Ta = numa_EL_Ta(end)/dena_EL_Ta(end-1);
% El
Kacc_des_EL_Ta = 1./Erp_Par_El;
Ki_EL_Ta = Kacc_des_EL_Ta./Kvel_now_EL_Ta;

F_EL_Ta = 10;
Zi_EL_Ta = real(s_pchoix_Ta)./F_EL_Ta;
Kp_EL_Ta = -Ki_EL_Ta./Zi_EL_Ta;

numPIEL = [Kp_EL_Ta Ki_EL_Ta];
denPIEL = [1 0];
disp('Fonction de transfer du PI (EL)')
PI_EL_Ta = tf(numPIEL,denPIEL)
% disp(['Fonction de transfert Élévation (Ta) avec AvPh et PI:'])
AvPh_FTEL_PI_Ta = AvPh_FT_EL_Ta*PI_EL_Ta;

% EL
if show_intermediary_graphs == 1
    figure('Name','Lieu de bode de EL original, avec AvPh et avec PI (Ta)')
    hold on
    margin(FT_EL)
    margin(AvPh_FT_EL_Ta)
    margin(AvPh_FTEL_PI_Ta)
    legend('FTEL','FTEL + AvPh','AvPh + FTEL + PI','Location','NorthEast')
    grid on
    hold off
end
% =======Coupe bande (Fin)========
if util_cb == 0
    AvPh_FTAZ_RePh_CB_Ta = AvPh_FTAZ_RePh_Ta*Butt_AZ;
    AvPh_FTEL_PI_CB_Ta = AvPh_FTEL_PI_Ta*Butt_EL;    
end


% =======================================================
% ======================Telescope B======================
% =======================================================
disp(' ')
disp('=======================================================')
disp('======================Telescope B======================')
disp('=======================================================')
disp(' ')

BW_Min_Tb = 10;         %En Rad/s
PM_Tb_avg = 50;         %En Degres
PM_Tb_min = PM_Tb_avg-1;%En Degres
PM_Tb_max = PM_Tb_avg+1;%En Degres

%Reponse a une rampe unitaire
Erp_Ramp_Tb = 0.005;    %En Degres
ts2pc_Erp_Ramp_Tb = 8;  %En secondes 

disp('Les spécifications du client pour le Téléscope B sont:')
disp(['PM désirée = ', num2str(PM_Tb_avg),' deg +/- 1 deg'])
disp(['BW minimale désirée = ', num2str(BW_Min_Tb),' rad/s'])
disp(['eRP désirée = ', num2str(Erp_Ramp_Tb)])
disp(' ')

% Traduction des spécifications
Zeta_des_Tb_avg = (1./2).*sqrt(tand(PM_Tb_avg).*sind(PM_Tb_avg));
Zeta_des_Tb_min = (1./2).*sqrt(tand(PM_Tb_min).*sind(PM_Tb_min));
Zeta_des_Tb_max = (1./2).*sqrt(tand(PM_Tb_max).*sind(PM_Tb_max));
Kvel_des_Tb = 1./Erp_Ramp_Tb;
wg_des_Tb = (BW_Min_Tb.*sqrt(sqrt(1 + 4.*(Zeta_des_Tb_avg.^4)) - (2.*(Zeta_des_Tb_avg.^2))))./(sqrt((1 - 2.*(Zeta_des_Tb_avg.^2)) + sqrt(4.*(Zeta_des_Tb_avg.^4) - 4.*(Zeta_des_Tb_avg.^2) + 2)));
disp('Les spécifications dérivées (paramètres standards) sont:')
disp(['Zeta désiré = ',num2str(Zeta_des_Tb_avg)])
disp(['Kvel désiré = ', num2str(Kvel_des_Tb)])
disp(['Fréquence de traverse en gain désirée (wg*) = ', num2str(wg_des_Tb),' rad/s'])
disp(' ')

% Lieu de Bode de la FT originale
if show_intermediary_graphs == 1
    figure('Name','Lieu de Bode de la FT AZ originale (Tb)')
    hold on
    margin(FT_AZ)
    xline(wg_des_Tb)
    grid on
    hold off

    figure('Name','Lieu de Bode de la FT EL originale (Tb)')
    hold on
    margin(FT_EL)
    xline(wg_des_Tb)
    grid on
    hold off
end
% Concevoir un compensateur AvPh
disp('=========Concevoir un compensateur AvPh (Bode)=========')
% Avec G(s) et wg* connus, on peut calculer ce gain K* qui donne l'amplitude
% du compensateur à cette fréquence de traverse. (7-29)
% disp('Compensateur avance-de-phase avec Bode')
[mag_AZ_Tb, pha_AZ_Tb] = bode(FT_AZ,wg_des_Tb); 
[mag_EL_Tb, pha_EL_Tb] = bode(FT_EL,wg_des_Tb); 
Kdes_AZ_Tb = 1/mag_AZ_Tb; 
Kdes_EL_Tb = 1/mag_EL_Tb; 
disp(['Kdes_AZ_Tb = ', num2str(Kdes_AZ_Tb), ' (K*)'])
disp(['Kdes_EL_Tb = ', num2str(Kdes_EL_Tb), ' (K*)'])

% Tracer le diagramme de Bode de la nouvelle FTBO, 
if show_intermediary_graphs == 1
    % AZ
    figure('Name',['FT AZ et Kdes*FT_AZ pour wg à ', num2str(wg_des_Tb), ' (Tb)'])
    hold on
    margin(FT_AZ)
    margin(Kdes_AZ_Tb*FT_AZ)
    legend('FTAZ','Kdes * FTAZ')
    grid on
    hold off
    % EL
    figure('Name',['FT EL et Kdes*FT_EL pour wg à ', num2str(wg_des_Tb), ' (Tb)'])
    hold on
    margin(FT_EL)
    margin(Kdes_EL_Tb*FT_EL)
    legend('FTEL','Kdes * FTEL')
    grid on
    hold off
end
% Calculer l'avance de phase requise à wg* (7-29)
warning('off', 'Control:analysis:MarginUnstable'); %Je sais que mon systeme est instable actuellement, je desactive temporairement 
[~, PM_des_AZ_Tb, ~, wgstarCheck_AZ] = margin(Kdes_AZ_Tb*FT_AZ);
[~, PM_des_EL_Tb, ~, wgstarCheck_EL] = margin(Kdes_EL_Tb*FT_EL);
warning('on', 'Control:analysis:MarginUnstable');
disp(' ')
deltaPhi_AZ_Tb = PM_Tb_avg - PM_des_AZ_Tb;
deltaPhi_EL_Tb = PM_Tb_avg - PM_des_EL_Tb;
disp(['L`avance de phase requise sans marge (AZ) est de ', num2str(deltaPhi_AZ_Tb), ' degres'])
disp(['L`avance de phase requise sans marge (EL) est de ', num2str(deltaPhi_EL_Tb), ' degres'])
disp(' ')
% AvPh AZ
disp('Caractéristiques du compensateur AvPh (AZ):')
alpha_AZ_Tb = (1 - sind(deltaPhi_AZ_Tb))./(1 + sind(deltaPhi_AZ_Tb));
TAvPh_AZ_Tb  = 1/(wg_des_Tb*sqrt(alpha_AZ_Tb));
KaAvPh_AZ_Tb = Kdes_AZ_Tb/sqrt(alpha_AZ_Tb);
numAvPh_AZ_Tb = [1 1/TAvPh_AZ_Tb]*KaAvPh_AZ_Tb;
denAvPh_AZ_Tb = [1 1/(alpha_AZ_Tb*TAvPh_AZ_Tb)];
disp(['paramètre alpha AZ = ',num2str(alpha_AZ_Tb)])
disp(['paramètre T (AZ) = ', num2str(TAvPh_AZ_Tb)])
disp(['gain Ka (EL) = ', num2str(KaAvPh_AZ_Tb)])
disp('Fonction du compensateur AvPh (AZ)')
AvPh_AZ_Tb = tf(numAvPh_AZ_Tb,denAvPh_AZ_Tb)
AvPh_FTAZ_Tb = AvPh_AZ_Tb * FT_AZ;
[numa_AZ_Tb,dena_AZ_Tb] = tfdata(AvPh_FTAZ_Tb,'v');
disp(' ')
% AvPh EL
disp('Caractéristiques du compensateur AvPh (EL):')
alpha_EL_Tb = (1 - sind(deltaPhi_EL_Tb))./(1 + sind(deltaPhi_EL_Tb));
TAvPh_EL_Tb  = 1/(wg_des_Tb*sqrt(alpha_EL_Tb));
KaAvPh_EL_Tb = Kdes_EL_Tb/sqrt(alpha_EL_Tb);
numAvPh_EL_Tb = [1 1/TAvPh_EL_Tb]*KaAvPh_EL_Tb;
denAvPh_EL_Tb = [1 1/(alpha_EL_Tb*TAvPh_EL_Tb)];
disp(['paramètre alpha EL = ',num2str(alpha_EL_Tb)])
disp(['paramètre T (EL) = ', num2str(TAvPh_EL_Tb)])
disp(['gain Ka (EL) = ', num2str(KaAvPh_EL_Tb)])
disp('Fonction du compensateur AvPh (EL)')
AvPh_EL_Tb = tf(numAvPh_EL_Tb,denAvPh_EL_Tb)
AvPh_FTEL_Tb = AvPh_EL_Tb * FT_EL;
[numa_EL_Tb,dena_EL_Tb] = tfdata(AvPh_FTEL_Tb,'v');
disp(' ')

% Lieu de Bode de la FT originale et FT avec AvPh
if show_intermediary_graphs == 1
    figure('Name','Lieu de Bode de la FT AZ originale et avec AvPh (Tb)')
    hold on
    margin(FT_AZ)
    margin(AvPh_FTAZ_Tb)
    legend('FTAZ','FTAZ + AvPh')
    grid on
    hold off

    figure('Name','Lieu de Bode de la FT EL originale et avec AvPh (Tb)')
    hold on
    margin(FT_EL)
    margin(AvPh_FTEL_Tb)
    legend('FTEL','FTEL + AvPh')
    grid on
    hold off
end

% Concevoir un compensateur RePh
disp('=========Concevoir un compensateur RePh=========')
% À partir des exigences sur l’erreur en régime permanent, déterminer le coefficient d'erreur en régime permanent désiré
Kr_RePh_Tb = 1; % Pour la conception, Kr = 1
F_RePh_Tb = 10;
T_RePh_Tb = F_RePh_Tb./wg_des_Tb;
Zr_RePh_Tb = -wg_des_Tb./F_RePh_Tb;

% Az
disp('Caractéristiques du compensateur RePh (AZ):')
Kvel_now_AZ_Tb = numa_AZ_Tb(end)/dena_AZ_Tb(end-1);
Kvel_des_AZ_Tb = 1./Erp_Ramp_Tb;
K_des_AZ_RePh_Tb = Kvel_des_AZ_Tb./Kvel_now_AZ_Tb;
beta_AZ_Tb = K_des_AZ_RePh_Tb./Kr_RePh_Tb;
Pr_RePh_AZ_Tb = (-1)./(beta_AZ_Tb.*T_RePh_Tb);
disp(['paramètre T = ', num2str(T_RePh_Tb)])
disp(['gain K* = ', num2str(K_des_AZ_RePh_Tb)])
numRe_AZ_Tb = [1 -Zr_RePh_Tb];
denRe_AZ_Tb = [1 -Pr_RePh_AZ_Tb];
disp('Fonction de transfer du RePh (AZ)')
RePh_AZ_Tb = tf(numRe_AZ_Tb,denRe_AZ_Tb)
AvPh_FTAZ_RePh_Tb = RePh_AZ_Tb*AvPh_FTAZ_Tb;

% EL
disp('Caractéristiques du compensateur RePh (EL):')
Kvel_now_EL_Tb = numa_EL_Tb(end)/dena_EL_Tb(end-1);
Kvel_des_EL_Tb = 1./Erp_Ramp_Tb;
K_des_EL_RePh_Tb = Kvel_des_EL_Tb./Kvel_now_EL_Tb;
beta_EL_Tb = K_des_EL_RePh_Tb./Kr_RePh_Tb;
Pr_RePh_EL_Tb = (-1)./(beta_EL_Tb.*T_RePh_Tb);
disp(['paramètre T = ', num2str(T_RePh_Tb)])
disp(['gain K* = ', num2str(K_des_EL_RePh_Tb)])
numRe_EL_Tb = [1 -Zr_RePh_Tb];
denRe_EL_Tb = [1 -Pr_RePh_EL_Tb];
disp('Fonction de transfer du RePh (EL)')
RePh_EL_Tb = tf(numRe_EL_Tb,denRe_EL_Tb)
AvPh_FTEL_RePh_Tb = RePh_EL_Tb*AvPh_FTEL_Tb;

% Lieu de Bode de la FT originale, avec AvPh et avec RePh
if show_intermediary_graphs == 1
    figure('Name','Lieu de Bode de la FT AZ originale avec AvPh et avec RePh (Tb)')
    hold on
    margin(FT_AZ)
    margin(AvPh_FTAZ_Tb)
    margin(AvPh_FTAZ_RePh_Tb)
    legend('FTAZ','FTAZ + AvPh','FTAZ + AvPh + RePh')
    grid on
    hold off

    figure('Name','Lieu de Bode de la FT EL originale avec AvPh et avec RePh (Tb)')
    hold on
    margin(FT_EL)
    margin(AvPh_FTEL_Tb)
    margin(AvPh_FTEL_RePh_Tb)
    legend('FTEL','FTEL + AvPh','FTEL + AvPh + RePh')
    grid on
    hold off
end
% =======Coupe bande (Fin)========
if util_cb == 0
    AvPh_FTAZ_RePh_CB_Tb = AvPh_FTAZ_RePh_Tb*Butt_AZ;
    AvPh_FTEL_RePh_CB_Tb = AvPh_FTEL_RePh_Tb*Butt_EL;    
end

%%
% =======================================================
% ===================Valeurs Initiales===================
% =======================================================
disp(' ')
disp('=======================================================')
disp('===================Valeurs Initiales===================')
disp('=======================================================')
disp(' ')

disp('======================Télescope A======================')
% disp(' ')
% ========AZ========
disp('========AZ========')
[num_Ini_AZ_Ta,den_Ini_AZ_Ta] = tfdata(AvPh_FTAZ_RePh_CB_Ta,'v');
t_AZ_Ta = [0:0.001:20]';
u_unit_AZ_Ta = ones(size(t_AZ_Ta)); % Échelon unitaire
u_ramp_AZ_Ta = t_AZ_Ta;
[num_FB_AZ_Ta, den_FB_AZ_Ta] = feedback(num_Ini_AZ_Ta,den_Ini_AZ_Ta,1,1);
% ou FTBF = feedback(FTBO,1)
ybf_unit_AZ_Ta = lsim(num_FB_AZ_Ta,den_FB_AZ_Ta,u_unit_AZ_Ta,t_AZ_Ta);
ybf_ramp_AZ_Ta = lsim(num_FB_AZ_Ta,den_FB_AZ_Ta,u_ramp_AZ_Ta,t_AZ_Ta);
figure('Name','Bode FTBF (AZ) (Ta)')
margin(tf(num_FB_AZ_Ta,den_FB_AZ_Ta))
grid on

figure('Name','Réponse a un echelon (AZ) (Ta)') % Cas de la réponse à l’échelon
hold on
plot(t_AZ_Ta,ybf_unit_AZ_Ta,'b', 'linewidth', 2)
plot([t_AZ_Ta(1); t_AZ_Ta(end)], 0.98*ybf_unit_AZ_Ta(end)*[1;1], 'r', 'linewidth', 1)
plot([t_AZ_Ta(1); t_AZ_Ta(end)], 1.02*ybf_unit_AZ_Ta(end)*[1;1], 'r', 'linewidth', 1)
grid on
hold off

% Retrouver mes paramètres
tmp_unit_AZ_Ta = lsiminfo(ybf_unit_AZ_Ta,t_AZ_Ta);
tmp_ramp_AZ_Ta = lsiminfo(ybf_ramp_AZ_Ta,t_AZ_Ta);
ts_unit_Fin_AZ_Ta = tmp_unit_AZ_Ta.SettlingTime;
ts_ramp_Fin_AZ_Ta = tmp_ramp_AZ_Ta.SettlingTime;
Mp_Fin_AZ_Ta = ((tmp_unit_AZ_Ta.Max - ybf_unit_AZ_Ta(end))/ybf_unit_AZ_Ta(end)).*100;
tp_Fin_AZ_Ta = tmp_unit_AZ_Ta.MaxTime;
Zeta_Fin_AZ_Ta = cosd(atand(-pi./(log(Mp_Fin_AZ_Ta./100))));
wa_Fin_AZ_Ta = pi./tp_Fin_AZ_Ta;
tr0_100_Fin_AZ = (pi - acos(Zeta_Fin_AZ_Ta))./wa_Fin_AZ_Ta;
Kvel_now_Fin_AZ_Ta = num_Ini_AZ_Ta(end)/den_Ini_AZ_Ta(end-1);
disp(['ts(2%) (ech) = ', num2str(ts_unit_Fin_AZ_Ta), ' s'])
disp(['Mp = ', num2str(Mp_Fin_AZ_Ta),' %'])
disp(['tr(0-100%) = ', num2str(tr0_100_Fin_AZ), ' s'])
disp(['ERP (Rampe) = ', num2str(1/Kvel_now_Fin_AZ_Ta), ' deg en Azimut'])
disp(['ts(2%) (ramp) = ', num2str(ts_ramp_Fin_AZ_Ta), ' s'])

figure('Name','Lieu de bode de AZ original, avec AvPh, avec RePh et avec CB (Ta)')
hold on
margin(FT_AZ)
margin(AvPh_FT_AZ_Ta)
margin(AvPh_FTAZ_RePh_Ta)
margin(AvPh_FTAZ_RePh_CB_Ta)
legend('FTAZ','FTAZ + AvPh','AvPh + FTAZ + RePh','AvPh + FTAZ + RePh + CB','Location','SouthWest')
grid on
hold off


% ========EL========
disp('========EL========')
[num_Ini_EL_Ta,den_Ini_EL_Ta] = tfdata(AvPh_FTEL_PI_CB_Ta,'v');
t_EL_Ta = [0:0.001:20]';
u_unit_EL_Ta = ones(size(t_EL_Ta)); % Échelon unitaire
u_para_EL_Ta = t_EL_Ta.^2/2; % Parabole unitaire
[num_FB_EL_Ta, den_FB_EL_Ta] = feedback(num_Ini_EL_Ta,den_Ini_EL_Ta,1,1);
% ou FTBF = feedback(FTBO,1)
ybf_unit_EL_Ta = lsim(num_FB_EL_Ta,den_FB_EL_Ta,u_unit_EL_Ta,t_EL_Ta);
ybf_para_EL_Ta = lsim(num_FB_EL_Ta,den_FB_EL_Ta,u_para_EL_Ta,t_EL_Ta);

figure('Name','Bode FTBF (EL) (Ta)')
margin(tf(num_FB_EL_Ta,den_FB_EL_Ta))
grid on

figure('Name','Réponse a un echelon (EL) (Ta)') % Cas de la réponse à l’échelon
hold on
plot(t_EL_Ta,ybf_unit_EL_Ta,'b', 'linewidth', 2)
plot([t_EL_Ta(1); t_EL_Ta(end)], 0.98*ybf_unit_EL_Ta(end)*[1;1], 'r', 'linewidth', 1)
plot([t_EL_Ta(1); t_EL_Ta(end)], 1.02*ybf_unit_EL_Ta(end)*[1;1], 'r', 'linewidth', 1)
grid on
hold off


% Retrouver mes paramètres
tmp_unit_EL_Ta = lsiminfo(ybf_unit_EL_Ta,t_EL_Ta);
tmp_para_EL_Ta = lsiminfo(ybf_para_EL_Ta,t_EL_Ta);
ts_unit_Fin_EL_Ta = tmp_unit_EL_Ta.SettlingTime;
ts_para_Fin_EL_Ta = tmp_para_EL_Ta.SettlingTime;
Mp_Fin_EL_Ta = ((tmp_unit_EL_Ta.Max - ybf_unit_EL_Ta(end))/ybf_unit_EL_Ta(end)).*100;
tp_Fin_EL_Ta = tmp_unit_EL_Ta.MaxTime;
Zeta_Fin_EL_Ta = cosd(atand(-pi./(log(Mp_Fin_EL_Ta./100))));
wa_Fin_EL_Ta = pi./tp_Fin_EL_Ta;
tr0_100_Fin_EL_Ta = (pi - acos(Zeta_Fin_EL_Ta))./wa_Fin_EL_Ta;
Kacc_now_Fin_EL_Ta = num_Ini_EL_Ta(end)/den_Ini_EL_Ta(end-2);
disp(['ts(2%) (ech) = ', num2str(ts_unit_Fin_EL_Ta), ' s'])
disp(['Mp = ', num2str(Mp_Fin_EL_Ta),' %'])
disp(['tr(0-100%) = ', num2str(tr0_100_Fin_EL_Ta), ' s'])
disp(['ERP (Parabole) = ', num2str(1/Kacc_now_Fin_EL_Ta), ' deg en Elevation'])
disp(['ts(2%) (para) = ', num2str(ts_para_Fin_EL_Ta), ' s'])

figure('Name','Lieu de bode de EL original, avec AvPh, avec PI et avec CB (Ta)')
hold on
margin(FT_EL)
margin(AvPh_FT_EL_Ta)
margin(AvPh_FTEL_PI_Ta)
margin(AvPh_FTEL_PI_CB_Ta)
legend('FTEL','FTEL + AvPh','AvPh + FTEL + PI','AvPh + FTEL + PI + CB','Location','SouthWest')
grid on
hold off


disp(' ')
disp('======================Télescope B======================')
% disp(' ')

% Retrouver mes paramètres
disp('========AZ========')
[num_Ini_AZ_Tb,den_Ini_AZ_Tb] = tfdata(AvPh_FTAZ_RePh_CB_Tb,'v');
[Gm_Ini_AZ_Tb,Pm_Ini_AZ_Tb,Wcp_Ini_AZ_Tb,Wg_Ini_AZ_Tb] = margin(AvPh_FTAZ_RePh_CB_Tb);
Zeta_Ini_AZ_Tb = (1./2).*sqrt(tand(Pm_Ini_AZ_Tb).*sind(Pm_Ini_AZ_Tb));
FTBF_AZ_Tb = feedback(AvPh_FTAZ_RePh_CB_Tb,1);
BW_Ini_AZ_Tb = bandwidth(FTBF_AZ_Tb);
delay_Ini_AZ_Tb = (Pm_Ini_AZ_Tb*pi)/(Wg_Ini_AZ_Tb*180);
Kvel_now_Fin_AZ_Tb = num_Ini_AZ_Tb(end)/den_Ini_AZ_Tb(end-1);
disp(['Le BW initial du systeme complet (AZ) est : ', num2str(BW_Ini_AZ_Tb), ' rad/s'])
disp(['La PM initial du systeme complet (AZ) est : ', num2str(Pm_Ini_AZ_Tb), ' deg'])
disp(['La GM initial du systeme complet (AZ) est : ', num2str(20*log10(Gm_Ini_AZ_Tb)), ' db'])
% disp(['La DM initial du systeme complet (AZ) est : ', num2str(delay_Ini_AZ_Tb), ' s'])
disp(['ERP (Rampe) = ', num2str(1/Kvel_now_Fin_AZ_Tb), ' deg en Azimut'])


t_AZ_Tb = [0:0.001:20]';
u_unit_AZ_Tb = ones(size(t_AZ_Tb)); % Échelon unitaire
[num_FB_AZ_Tb, den_FB_AZ_Tb] = feedback(num_Ini_AZ_Tb,den_Ini_AZ_Tb,1,1);
% ou FTBF = feedback(FTBO,1)
ybf_AZ_Tb = lsim(num_FB_AZ_Tb,den_FB_AZ_Tb,u_unit_AZ_Tb,t_AZ_Tb);
tmp_AZ_Tb = lsiminfo(ybf_AZ_Tb,t_AZ_Tb);
ts_Fin_AZ_Tb = tmp_AZ_Tb.SettlingTime;
disp(['ts(2%) (AZ) = ', num2str(ts_Fin_AZ_Tb), ' s'])

figure('Name','Réponse a un echelon (AZ) (Tb)') % Cas de la réponse à l’échelon
hold on
plot(t_AZ_Tb, ybf_AZ_Tb,'b', 'linewidth', 2)
plot([t_AZ_Tb(1); t_AZ_Tb(end)], 0.98*ybf_AZ_Tb(end)*[1;1], 'r', 'linewidth', 1)
plot([t_AZ_Tb(1); t_AZ_Tb(end)], 1.02*ybf_AZ_Tb(end)*[1;1], 'r', 'linewidth', 1)
grid on
hold off


% Lieu de Bode de la FT originale, avec AvPh et avec RePh
figure('Name','Lieu de Bode de la FT AZ originale avec AvPh,avec RePh et avec CB (Tb)')
hold on
margin(FT_AZ)
margin(AvPh_FTAZ_Tb)
margin(AvPh_FTAZ_RePh_Tb)
margin(AvPh_FTAZ_RePh_CB_Tb)
legend('FTAZ','FTAZ + AvPh','FTAZ + AvPh + RePh','FTAZ + AvPh + RePh + CB','Location','SouthWest')
grid on
hold off


disp('========EL========')
[num_Ini_EL_Tb,den_Ini_EL_Tb] = tfdata(AvPh_FTEL_RePh_CB_Tb,'v');
[Gm_Ini_EL_Tb,Pm_Ini_EL_Tb,Wcp_Ini_EL_Tb,Wg_Ini_EL_Tb] = margin(AvPh_FTEL_RePh_CB_Tb);
Zeta_Ini_EL_Tb = (1./2).*sqrt(tand(Pm_Ini_EL_Tb).*sind(Pm_Ini_EL_Tb));
FTBF_EL_Tb = feedback(AvPh_FTEL_RePh_CB_Tb,1);
BW_Ini_EL_Tb = bandwidth(FTBF_EL_Tb);
delay_Ini_EL_Tb = (Pm_Ini_EL_Tb*pi)/(Wg_Ini_EL_Tb*180);
Kvel_now_Fin_EL_Tb = num_Ini_EL_Tb(end)/den_Ini_EL_Tb(end-1);
disp(['Le BW initial du systeme complet (EL) est : ', num2str(BW_Ini_EL_Tb), ' rad/s'])
disp(['La PM initial du systeme complet (EL) est : ', num2str(Pm_Ini_EL_Tb), ' deg'])
disp(['La GM initial du systeme complet (EL) est : ', num2str(20*log10(Gm_Ini_EL_Tb)), ' db'])
% disp(['La DM initial du systeme complet (EL) est : ', num2str(delay_Ini_EL_Tb), ' s'])
disp(['ERP (Rampe) = ', num2str(1/Kvel_now_Fin_EL_Tb), ' deg en Élévation'])

t_EL_Tb = [0:0.001:20]';
u_unit_EL_Tb = ones(size(t_EL_Tb)); % Échelon unitaire
[num_FB_EL_Tb, den_FB_EL_Tb] = feedback(num_Ini_EL_Tb,den_Ini_EL_Tb,1,1);
% ou FTBF = feedback(FTBO,1)
ybf_EL_Tb = lsim(num_FB_EL_Tb,den_FB_EL_Tb,u_unit_EL_Tb,t_EL_Tb);
tmp_EL_Tb = lsiminfo(ybf_EL_Tb,t_EL_Tb);
ts_Fin_EL_Tb = tmp_EL_Tb.SettlingTime;
disp(['ts(2%) (EL) = ', num2str(ts_Fin_EL_Tb), ' s'])

figure('Name','Réponse a un echelon (EL) (Tb)') % Cas de la réponse à l’échelon
hold on
plot(t_EL_Tb, ybf_EL_Tb,'b', 'linewidth', 2)
plot([t_EL_Tb(1); t_EL_Tb(end)], 0.98*ybf_EL_Tb(end)*[1;1], 'r', 'linewidth', 1)
plot([t_EL_Tb(1); t_EL_Tb(end)], 1.02*ybf_EL_Tb(end)*[1;1], 'r', 'linewidth', 1)
grid on
hold off

% Lieu de Bode de la FT originale, avec AvPh et avec RePh
figure('Name','Lieu de Bode de la FT EL originale avec AvPh, avec RePh et avec CB (Tb)')
hold on
margin(FT_EL)
margin(AvPh_FTEL_Tb)
margin(AvPh_FTEL_RePh_Tb)
margin(AvPh_FTEL_RePh_CB_Tb)
legend('FTEL','FTEL + AvPh','FTEL + AvPh + RePh','FTEL + AvPh + RePh + CB','Location','SouthWest')
grid on
hold off

% 
% %%
% % =======================================================
% % ========================Itérations=====================
% % =======================================================
% disp(' ')
% disp('=======================================================')
% disp('=======================Itérations======================')
% disp('=======================================================')
% disp(' ')
% %Criteres Maximaux Ta
% Mp_Ta_Max_Az = 30;
% Mp_Ta_Max_El = 35;
% ts2pc_Ta_Max = 1.20;   %En secondes
% AttVibFB_Ta = -15;     %En dB
% % Traduction des spécifications
% Zeta_Ta_Max_Az = cosd(atand(-pi./(log(Mp_Ta_Max_Az./100))));
% Zeta_Ta_Max_El = cosd(atand(-pi./(log(Mp_Ta_Max_El./100))));
% wn_ts2pc_Ta_Max = 4./(Zeta_Ta.*ts2pc_Ta_Max);
% 
% 
% 
% 
% 
% 
% 
% disp(' ')
% disp('======================Télescope B======================')
% %Criteres 
% BW_Ite_Tb = 10;         %En Rad/s
% PM_Ite_AZ_Tb = 60;
% Zeta_des_Ite_Tb = (1./2).*sqrt(tand(PM_Ite_AZ_Tb).*sind(PM_Ite_AZ_Tb))
% wg_des_Ite_Tb = (BW_Ite_Tb.*sqrt(sqrt(1 + 4.*(Zeta_des_Ite_Tb.^4)) - (2.*(Zeta_des_Ite_Tb.^2))))./(sqrt((1 - 2.*(Zeta_des_Ite_Tb.^2)) + sqrt(4.*(Zeta_des_Ite_Tb.^4) - 4.*(Zeta_des_Ite_Tb.^2) + 2)))
% 
% % Traduction des spécifications
% [mag_AvPh_Ite_AZ_Tb, pha_AvPh_Ite_AZ_Tb] = bode(FT_AZ,wg_des_Ite_Tb); 
% Kdes_AvPh_Ite_AZ_Tb = 1/mag_AvPh_Ite_AZ_Tb; 
% disp(' ')
% 
% disp('========AZ========')
% disp(['Kdes_AvPh_Ite_AZ_Tb = ', num2str(Kdes_AvPh_Ite_AZ_Tb), ' (K*)'])
% [~, PM_des_Ite_AZ_Tb, ~, wgstarCheck_Ite_AZ] = margin(Kdes_AvPh_Ite_AZ_Tb*FT_AZ)
% deltaPhi_Ite_AZ_Tb = PM_Ite_AZ_Tb - PM_des_Ite_AZ_Tb
% alpha_Ite_AZ_Tb = (1 - sind(deltaPhi_Ite_AZ_Tb))./(1 + sind(deltaPhi_Ite_AZ_Tb));
% TAvPh_Ite_AZ_Tb  = 1/(wg_des_Ite_Tb*sqrt(alpha_Ite_AZ_Tb));
% KaAvPh_AZ_Tb = Kdes_AvPh_Ite_AZ_Tb/sqrt(alpha_Ite_AZ_Tb);
% numAvPh_Ite_AZ_Tb = [1 1/TAvPh_Ite_AZ_Tb]*Kdes_AvPh_Ite_AZ_Tb;
% denAvPh_Ite_AZ_Tb = [1 1/(alpha_Ite_AZ_Tb*TAvPh_Ite_AZ_Tb)];
% AvPh_Ite_AZ_Tb = tf(numAvPh_Ite_AZ_Tb,denAvPh_Ite_AZ_Tb)
% % AvPh_FTAZ_Tb = AvPh_Ite_AZ_Tb * FT_AZ;
% % [numa_AZ_Tb,dena_AZ_Tb] = tfdata(AvPh_FTAZ_Tb,'v')
% 
% 
% figure('Name','Nouvelle AvPh AZ (TB)')
% margin((AvPh_Ite_AZ_Tb*FT_AZ*RePh_AZ_Tb*Butt_AZ))
% grid on
% 
% disp('========EL========')
% [mag_AvPh_Ite_EL_Tb, pha_AvPh_Ite_EL_Tb] = bode(FT_EL,wg_des_Ite_Tb); 
% Kdes_AvPh_Ite_EL_Tb = 1/mag_AvPh_Ite_EL_Tb; 
% disp(['Kdes_AvPh_Ite_EL_Tb = ', num2str(Kdes_AvPh_Ite_EL_Tb), ' (K*)'])
% [~, PM_des_Ite_EL_Tb, ~, wgstarCheck_Ite_EL] = margin(Kdes_AvPh_Ite_EL_Tb*FT_EL)
% deltaPhi_Ite_EL_Tb = PM_Tb_avg - PM_des_Ite_EL_Tb
% 
% 
% 
% 
% 
