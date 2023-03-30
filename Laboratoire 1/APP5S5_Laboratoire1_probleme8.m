% APP5 S5 Laboratoire 1 
% Probleme 8 : Conception de AvPh et RePh avec le diagramme de Bode (Méthode 1) (7.5.3, 7.6.3)
close all
clear
clc

%Valeurs de depart
numGs = [1 4];
denGs = [1 8 3 0];
disp('La fonction de transfert originale est:')
FTBO = tf(numGs,denGs)

%Le système en boucle ouverte ci-dessous est asservi avec une rétroaction de gain unitaire. 
%Les critères de performance à atteindre sont :
PM = 50;            %degres
BWbf = 5.2;         %rad/s
ErpRampe = 0.005;
disp('Les spécifications du client sont:')
disp(['PM désirée = ', num2str(PM),' deg'])
disp(['BW désirée = ', num2str(BWbf),' rad/s'])
disp(['eRP désirée = ', num2str(ErpRampe)])
disp(' ')

% spécifications dérivées
zetastar = (1./2).*sqrt(tand(PM).*sind(PM));
Kvelstar = 1./ErpRampe;
wgstar = (BWbf.*sqrt(sqrt(1 + 4.*(zetastar.^4)) - (2.*(zetastar.^2))))./(sqrt((1 - 2.*(zetastar.^2)) + sqrt(4.*(zetastar.^4) - 4.*(zetastar.^2) + 2)));
disp('Les spécifications dérivées (paramètres standards) sont:')
disp(['zeta désiré = ',num2str(zetastar)])
disp(['Kvel désiré = ', num2str(Kvelstar)])
disp(['Fréquence de traverse en gain désirée (wg*) = ', num2str(wgstar),' rad/s'])
disp(' ')

% Lieu de Bode de la FT originale
figure('Name','Lieu de Bode de la FT originale')
hold on
margin(FTBO)
xline(wgstar)
grid on
hold off
% Plusieurs types de compensateurs en cascade sont possibles en utilisant le diagramme de Bode. 
% La methode successive de faire un AvPh (PD) suivi d’un RePh (PI) est appelée Méthode 1 dans les notes de cours.
% Dans le présent problème, on utilise cette approche avec les Méthodes 1 de l’AvPh et du RePh. 
% Cette approche est identique à la conception dans le lieu des racines : on règle le régime transitoire en 
% premier avec un AvPh (ou un PD) et ensuite on règle le régime permanent avec un RePh (ou un PI)


% Calcul du gain K* (Kdes) pour obtenir la fréquence de traverse à wg* (wg_des)
% C'est l'inverse de l'amplitude de la FTBO à wg*
[MAG, PHA] = bode(FTBO,wgstar);
Kstar = 1./MAG;
% Kstar = abs( polyval(denGs,j*wgstar)/polyval(numGs,j*wgstar));
disp(['Kdes = ', num2str(Kstar), ' (K*)'])

% Calcul de la phase requise pour obtenir la PM désirée avec marge (5 deg) pour RePh
[~, PMstar, ~, wgstarCheck] = margin(Kstar*FTBO);
deltaPhi = PM - PMstar;
disp(['L`avance de phase requise sans marge est de ', num2str(deltaPhi), ' rad/s'])
disp(['L`avance de phase requise avec marge est de ', num2str(deltaPhi+5), ' rad/s']) 
disp('On ajoute 5 deg de marge pour prévoir la perte de phase causée par le RePh')
disp(' ')

figure('Name',['FT et Kdes*FT pour wg à', num2str(wgstar)])
margin(numGs,denGs)
hold on
title(['FT et Kstar*FT pour wg à', num2str(wgstar)])
margin(Kstar*FTBO)
% title(['FT et Kdes*FT pour wg à', num2str(wgstar)])
legend('FT originale','Kstar * FT')
grid on
hold off

Phi_m = deltaPhi+5;
alphaAvPh = (1-sind(Phi_m))/(1+sind(Phi_m));
TAvPh  = 1/(wgstar*sqrt(alphaAvPh));
KaAvPh = Kstar/sqrt(alphaAvPh);

disp('Caractéristiques du compensateur AvPh:')
disp(['paramètre alpha = ',num2str(alphaAvPh)])
disp(['paramètre T = ', num2str(TAvPh)])
disp(['gain Ka = ', num2str(KaAvPh)])
disp('Compensateur avance-de-phase avec Bode')
numAvPh = [1 1/TAvPh]*KaAvPh;
denAvPh = [1 1/(alphaAvPh*TAvPh)];
FTa = tf(numAvPh,denAvPh)
disp('Compensateur avance-de-phase x FT originale avec Bode (méthode 1)')
FTaFTBO = FTa*FTBO

[GMftfta,PMftfta,wss,wgftfta] = margin(FTaFTBO)

% Faire la conception d’un AvPh qui rencontre les critères PM et BW (avec la Méthode 1, Notes JdeL section 
% 7.5.3) suivi d’un RePh en cascade (avec la Methode 1, Notes JdeL section 7.6.3), qui maintient la BW et 
% donc 𝜔𝑔, pour rencontrer le critère d’erreur en RP. Prévoir une surcompensation de 5 deg dans la conception 
% de l’AvPh en prévision du RePh. A noter qu’on veut maintenir 𝜔𝑔 (qui dépend de 𝜔𝑛) et PM (qui dépend de 
% 𝜁) deja obtenus avec l’AvPh tout comme avec le lieu des racines on voulait conserver les pôles désirés (𝜁, 𝜔𝑛).








% 
% Notes de cours : La Méthode 1 est utilisée quand les performances demandées sont exprimées en termes de marge de phase 
% et de bande passante BW (ou équivalent : fréquence de traverse en gain wg,fréquence naturelle wn).
%
% 7.5.3 Conception de l’AvPh simple et double par le diagramme de Bode (Méthode 1) 
% La Méthode 1 est utilisée quand les performances demandées sont exprimées en termes de marge de phase 
% et de bande passante BW (ou équivalent : fréquence de traverse en gain wg*, fréquence naturelle wn)
% 
% * Il a été mentionné (à la Section 7.3) que l’utilisation du diagramme de Bode pour la conception de 
% compensateurs peut se faire avec l’approche nominale de ce document, c’est-à-dire (1) avec deux 
% compensateurs, un AvPh (ou PD) pour régler le régime transitoire et ensuite un RePh (ou PI) en cascade 
% pour régler le régime permanent ou (2) un seul compensateur qui fait les deux en même temps. L’approche 
% (1) est appelée la Méthode 1 et est présentée dans cette section.
% 
% * La Méthode 1 est en parallèle exact avec la conception avec le lieu des racines telle que présentée à la 
% Section 7.4 : on traite les spécifications sur le régime transitoire en premier avec un AvPh (ou PD) et les 
% spécifications sur le régime permanent ensuite avec un RePh (ou PI).
% 
% *Le principe de la compensation AvPh avec le diagramme de Bode consiste à modifier la réponse en 
% fréquence d’une fonction de transfert pour que sa marge de phase (PM) et sa bande passante (BW) 
% rencontre les performances désirées du client. Les spécifications sur la BW sont parfois exprimées de 
% façon équivalente par des spécifications sur la fréquence de traverse en
% gain wg* ou la frequence naturelle wn).
%
% À partir de la fonction de transfert d’un ordre 2 standard, la Section 7.1.3 démontre qu’il y a un lien 
% analytique direct entre les paramètres de performance PM, BW (ou wg*)et les paramètres standards Zeta et wn. 
% En effet, PM est relié par une seule équation a Zeta et BW (ou wg*) est lié aux deux Zeta et wn).
% On peut donc dire que Mp et PM sont dans la meme categorie : ils ne dependent que de Zeta.
% On peut aussi dire que tp,ts,tr et BW sont dans la meme categorie : ils dependent de zeta et wn.
%
% Principe de base : On utilise le fait que le diagramme de Bode est le tracé de l’amplitude et de la phase 
% de la FTBO à partir duquel la marge de phase PM peut être calculée à la fréquence de traverse en gain wg. 
%









