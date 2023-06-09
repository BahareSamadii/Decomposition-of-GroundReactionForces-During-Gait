% % Decomposition of ground reaction forces (GRF) during gait using parametric
% curve fitting model
% Phase detection from single to double stance phase : Written by Prof.
% Maxime Raison
% % % % 
%  Decomposition of GRF : Written by Bahare Samadi (M.Sc. Student)
% % %-----------------------------------------------------
%
%  Rehabilitation Chair Applied to Pediatrics (RECAP)
%
%      Polytechnique Montreal and Ste-Justine UHC
% 
%	     %
%        Evaluation laboratory: The motion laboratory of Marie-enfant
%        rehabilitation centre

%                  ----------------
% 
%   Informations:
%
%	- Echantillonnage: 100 Hz - filtre adaptatif
%   - Using two separate force platforms
%   - Fx = Medio-lateral forces
%   - Fy = Antero-posterior forces
%   - Fz = Vertical forces


% % % Colors of the figures:
%R : Recorded GRF
%C : Sine
%m : Spline
%g : Poly
%blue : Sine-Sigmoid
% Orange : sigmoid


% Decomposition of Fx
close all;
clear all;
clc



% % % % % % % % % % % % % % % 
% % % % %  
% % % % % % % % % % % % % % % 

% % % Reading recorded GRF
PTF1 = xlsread('26_03_F.xls');
R_Fx1 = PTF1 (:,1);
R_Fy1 = PTF1 (:,2);
R_Fz1 = PTF1 (:,3);

R_Fx2 = PTF1 (:,7);
R_Fy2 = PTF1 (:,8); 
R_Fz2 = PTF1 (:,9);


R_Fx3 = PTF1 (:,13);
R_Fy3 = PTF1 (:,14);
R_Fz3 = PTF1 (:,15);


R_Mx1 = PTF1 (:,4);
R_My1 = PTF1 (:,5);
R_Mz1 = PTF1 (:,6);

R_Mx2 = PTF1 (:,10);
R_My2 = PTF1 (:,11);
R_Mz2 = PTF1 (:,12);

R_Mx3 = PTF1 (:,16);
R_My3 = PTF1 (:,17);
R_Mz3 = PTF1 (:,18);

%==========================================================================
%
% 1. Interface utilisateur :
%
%==========================================================================

echantillons = (size(PTF1,1));

%-----------------------------------------------------
% choix pour l'affichage :
%-----------------------------------------------------
% % Utile pour calculs détection sauts Delta COP :

affichage_inf = 1;
affichage_sup = echantillons;


%==========================================================================
%
% 2. Calcul CoP Pour chaque platform de force et le CoP globale:
%
%==========================================================================
Fztot=R_Fz1+R_Fz2+R_Fz3;
Fxtot=R_Fx1+R_Fx2+R_Fx3;

COPxp1=(R_My1)./Fztot;
COPyp1=(R_Mx1)./Fztot;
% 
COPxp2=(R_My2)./Fztot;
COPyp2=(R_Mx2)./Fztot;
% 
COPxp3=(R_My3)./Fztot;
COPyp3=(R_Mx3)./Fztot;
% 
COPxp1_adap=COPxp1+333.5;
COPyp1_adap=COPyp1+352.5;
% 
COPxp2_adap=COPxp2+106;
COPyp2_adap=COPyp2+858.5;
% 
COPxp3_adap=COPxp3+559;
COPyp3_adap=COPyp3+956.5;



COPx = ((COPxp1_adap .* R_Fz1 + COPxp2_adap.*R_Fz2 + COPxp3_adap.*R_Fz3)./ (R_Fz1+R_Fz2+R_Fz3))./1000; 
COPy = ((COPyp1_adap .* R_Fz1 + COPyp2_adap.*R_Fz2 + COPyp3_adap.*R_Fz3)./ (R_Fz1+R_Fz2+R_Fz3))./1000;


%-----------------------------------------------------
% Définition DeltaCoP : 
%-----------------------------------------------------

Delta_square_COP1 = zeros(1,echantillons);
Delta_COP1 = zeros(1,echantillons);

for i = 1:size(COPy,1)-1,
    Delta_COPy1(i) = COPy(i+1)-COPy(i);
    Delta_COPx1(i) = COPx(i+1)-COPx(i);
    Delta_square_COP1(i) = Delta_COPx1(i)^2 + Delta_COPy1(i)^2;
    Delta_COP1 = sqrt(Delta_square_COP1);
end

Delta_square_COP5 = zeros(1,echantillons);
Delta_COP5 = zeros(1,echantillons);
for i = 1:size(COPy,1)-5,  % what is the 5, and why do we consider the 5th?
    Delta_COPy5(i) = COPy(i+5)-COPy(i);
    Delta_COPx5(i) = COPx(i+5)-COPx(i);
    Delta_square_COP5(i) = Delta_COPx5(i)^2 + Delta_COPy5(i)^2;
    Delta_COP5 = sqrt(Delta_square_COP5);
end


%--------------------------------------------------------------------------------------------
% Détection du changement de DeltaCoP pour detection de phases Single stance et Double Stance: 
%--------------------------------------------------------------------------------------------

clear Lim;
Lim =ones(1,3)*affichage_sup;  % We need Lim(1) to Lim(6) for spline interpolation, so we put 5 to repeat the last index three times
ind_lim = 1;
stop_iter = 0;

Cf = 0.003;
Cl = Cf;
Cp = 5*Cf;

j_ini=400; %4s
% j = j_ini; % au lieu de 1 (0.2 s) % You can get rid of this if everything is ok. In case that the first samples are not appropriate(noise or ...) we can remove them and start with the sample number 20, for example
j=j_ini;
while stop_iter < 1
%     COP5(j) is T1 in paper
    while (j< affichage_sup)& (Delta_COP5(j)< Cp); % (j < affichage_sup-constr)&(Delta_COP5(j) < Cp) %0.02 %0.0387 %0.02 % saut: seuil à fixer! sqrt(0.0015); 5 * 3 * std(COP at rest, so the noise of COP)
    j = j+1;
    end
    %     COP1(j) is T2 in paper
    while (j< affichage_sup)& (Delta_COP1(j)< Cf); %(j < affichage_sup-constr)&(Delta_COP1(j)< Cf) %/Delta_COP1(j-1) < 3 %0.0173 % 0.001 - 0.05 % should be higher if detected more S/D stances than expected, due to noise; should be lower if no noise, but you did not detect too late % Delta_COP1(j)/Delta_COP1(j-1) < 3 %(3=CF),CL pas: seuil à fixer!   < 0.015 % 3 * std(COP at rest, so the noise of COP)
        j = j+1;
    end
%     j = j -1;
%     j
%     plot(j/100-deb_aff,COPy(j),'k.');
    Lim(ind_lim) = j-1;%+1 % attention: point precedent retenu!
    ind_lim = ind_lim +1;
    stop_iter = 1;
end
    
j_fin=800; %8s
% j = j_fin; % au lieu de 1 (0.2 s) % You can get rid of this if everything is ok. In case that the first samples are not appropriate(noise or ...) we can remove them and start with the sample number 20, for example
j=j_fin;
stop_iter=0;
while stop_iter < 1
%     COP5(j) is T1 in paper
    while (j> affichage_inf)& (Delta_COP5(j)< Cp); % (j < affichage_sup-constr)&(Delta_COP5(j) < Cp) %0.02 %0.0387 %0.02 % saut: seuil à fixer! sqrt(0.0015); 5 * 3 * std(COP at rest, so the noise of COP)
    j = j-1;
    end
    %     COP1(j) is T2 in paper
    while (j> affichage_inf)& (Delta_COP1(j)< Cf); %(j < affichage_sup-constr)&(Delta_COP1(j)< Cf) %/Delta_COP1(j-1) < 3 %0.0173 % 0.001 - 0.05 % should be higher if detected more S/D stances than expected, due to noise; should be lower if no noise, but you did not detect too late % Delta_COP1(j)/Delta_COP1(j-1) < 3 %(3=CF),CL pas: seuil à fixer!   < 0.015 % 3 * std(COP at rest, so the noise of COP)
        j = j-1;
    end
%     j = j -1;
%     j
%     plot(j/100-deb_aff,COPy(j),'k.');
    Lim(ind_lim) = j+2;%+1 % attention: point suivant retenu!
    ind_lim = ind_lim +1;
    stop_iter = 1;
end
%     
%     
Lim; % saut, pas, saut, pas, saut ...

Lim = [1 Lim]; % echantillons];
sizelim = size(Lim,2);

frequency=100;

%--------------------------------------------------------------------------------------------
% La détection de la plate-forme de force engagée par les participants, deux des trois plates-formes sont embauchés à chaque cycle
%--------------------------------------------------------------------------------------------
    if R_Fx1(1)~=0
       Fx_gauche=R_Fx1;
        if R_Fx2(Lim(2))~=0
            Fx_droite=R_Fx2;
        else
            Fx_droite=R_Fx3;
        end
    elseif R_Fx3(1)~=0
            Fx_gauche=R_Fx3;
            Fx_droite=R_Fx1;
    else Fx_gauche=R_Fx2;
            Fx_droite=R_Fx1;
    end           

    
    if R_Fy1(1)~=0
       Fy_gauche=R_Fy1;
        if R_Fy2(Lim(2))~=0
            Fy_droite=R_Fy2;
        else
            Fy_droite=R_Fy3;
        end
    elseif R_Fy3(1)~=0
            Fy_gauche=R_Fy3;
            Fy_droite=R_Fy1;
    else Fy_gauche=R_Fy2;
            Fy_droite=R_Fy1;
    end           

    
    if R_Fz1(1)~=0
       Fz_gauche=R_Fz1;
        if R_Fz2(Lim(2))~=0
            Fz_droite=R_Fz2;
        else
            Fz_droite=R_Fz3;
        end
    elseif R_Fz3(1)~=0
            Fz_gauche=R_Fz3;
            Fz_droite=R_Fz1;
    else Fz_gauche=R_Fz2;
            Fz_droite=R_Fz1;
    end           

%--------------------------------------------------------------------------------------------
% Appliquer butterworth low-pass filtre afin d'enlever les bruits, cut-off frequency = 4 Hz
%--------------------------------------------------------------------------------------------
    
Wn=0.08;
Order = 4;
[bb,aa] = butter(Order,Wn,'low');

Fx_gauche1=filtfilt(bb,aa,Fx_gauche(Lim(1):Lim(3)-1)); 
Fx_gauche2=filtfilt(bb,aa,Fx_gauche(Lim(3):Lim(4))); 
Fx_gauche=[(Fx_gauche1);(Fx_gauche2)];

Fx_droite1=filtfilt(bb,aa,Fx_droite(Lim(1):Lim(2))); 
Fx_droite2=filtfilt(bb,aa,Fx_droite(Lim(2)+1:Lim(4)));
Fx_droite=[(Fx_droite1);(Fx_droite2)];


Fy_gauche1=filtfilt(bb,aa,Fy_gauche(Lim(1):Lim(3)-1)); 
Fy_gauche2=filtfilt(bb,aa,Fy_gauche(Lim(3):Lim(4))); 
Fy_gauche=[(Fy_gauche1);(Fy_gauche2)];

Fy_droite1=filtfilt(bb,aa,Fy_droite(Lim(1):Lim(2))); 
Fy_droite2=filtfilt(bb,aa,Fy_droite(Lim(2)+1:Lim(4))); 
Fy_droite=[(Fy_droite1);(Fy_droite2)];


Fz_gauche1=filtfilt(bb,aa,Fz_gauche(Lim(1):Lim(3)-1)); 
Fz_gauche2=filtfilt(bb,aa,Fz_gauche(Lim(3):Lim(4))); 
Fz_gauche=[(Fz_gauche1);(Fz_gauche2)];

Fz_droite1=filtfilt(bb,aa,Fz_droite(Lim(1):Lim(2))); 
Fz_droite2=filtfilt(bb,aa,Fz_droite(Lim(2)+1:Lim(4))); 
Fz_droite=[(Fz_droite1);(Fz_droite2)];

% Forces totales :
sizePTF = echantillons;            


% Forces totales :
R_FX = Fx_gauche + Fx_droite;
R_FY = Fy_gauche + Fy_droite;
R_FZ = Fz_gauche + Fz_droite;

% 
% Forces totales :
sizePTF = echantillons; 

R_F = sqrt(R_FX.^2 + R_FY.^2 + R_FZ.^2);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % on impose le choix pied gauche ou droit (pour calcul d'erreur):
PiedGauche = 1;  % 1er pied : gauche OK pour 0002~aa~Walking 01.tdf.
PiedDroit = 0;

%==========================================================================
%
% 1. Interface utilisateur :

%==========================================================================
% supprimé...
% fprintf('\n...\n')

%-----------------------------------------------------
% choix pour l'affichage :
%-----------------------------------------------------
% Utile pour calculs détection sauts Delta COP :

sampling=frequency;
% A changer si on veut particulariser abscisse (temps) figures :
deb_aff =( Lim(2)-20)/frequency;
fin_aff = (Lim(3)+20)/100; 

%==========================================================================
%
% Découplage des forces :
%
%==========================================================================
% % % % % Force X


%-----------------------------------------------------
% Interpolation linéaire :
%-----------------------------------------------------

clear R_FXG
clear R_FXD
R_FXG = zeros(sizePTF,1);
R_FXD = zeros(sizePTF,1);

if PiedGauche > PiedDroit % détection début (pied G ou D)
  for var1 =1:4:sizelim-1,
    R_FXG(Lim(var1):Lim(var1+1)) = R_FX(Lim(var1):Lim(var1+1),1);
    R_FXD(Lim(var1):Lim(var1+1)) = zeros(Lim(var1+1)-Lim(var1)+1,1);
  end
  for var1 =2:4:sizelim-1,
    R_FXG(Lim(var1)+0:Lim(var1+1)-1) = (-R_FX(Lim(var1)+0,1))*([Lim(var1)+0:Lim(var1+1)-1]-Lim(var1+1)) / (Lim(var1+1)-(Lim(var1)+0));
    R_FXD(Lim(var1)+0:Lim(var1+1)-1) = (R_FX(Lim(var1+1),1))*([Lim(var1)+0:Lim(var1+1)-1]-(Lim(var1)+0)) / (Lim(var1+1)-(Lim(var1)+0));  
  end
  for var1 =3:4:sizelim-1,
    R_FXD(Lim(var1):Lim(var1+1)) = R_FX(Lim(var1):Lim(var1+1),1);
    R_FXG(Lim(var1):Lim(var1+1)) = zeros(Lim(var1+1)-Lim(var1)+1,1);
  end
  for var1 =4:4:sizelim-1,
    if var1 < 8,
    R_FXD(Lim(var1)+0:Lim(var1+1)-1) = (-R_FX(Lim(var1)+0,1))*([Lim(var1)+0:Lim(var1+1)-1]-Lim(var1+1)) / (Lim(var1+1)-(Lim(var1)+0)); 
    R_FXG(Lim(var1)+0:Lim(var1+1)-1) = (R_FX(Lim(var1+1),1))*([Lim(var1)+0:Lim(var1+1)-1]-(Lim(var1)+0)) / (Lim(var1+1)-(Lim(var1)+0));
    end

  end
else
  for var1 = 1:4:sizelim-1, %(si pied D d'abord)
    R_FXD(Lim(var1):Lim(var1+1)) = R_FX(Lim(var1):Lim(var1+1),1);
    R_FXG(Lim(var1):Lim(var1+1)) = zeros(Lim(var1+1)-Lim(var1)+1,1);
  end
  for var1 = 2:4:sizelim-1,
    R_FXD(Lim(var1)+0:Lim(var1+1)-1) = (-R_FX(Lim(var1)+0,1))*([Lim(var1)+0:Lim(var1+1)-1]-Lim(var1+1)) / (Lim(var1+1)-(Lim(var1)+0));  
    R_FXG(Lim(var1)+0:Lim(var1+1)-1) = (R_FX(Lim(var1+1),1))*([Lim(var1)+0:Lim(var1+1)-1]-(Lim(var1)+0)) / (Lim(var1+1)-(Lim(var1)+0));  
  end
  for var1 = 3:4:sizelim-1,
    R_FXG(Lim(var1):Lim(var1+1)) = R_FX(Lim(var1):Lim(var1+1),1);
    R_FXD(Lim(var1):Lim(var1+1)) = zeros(Lim(var1+1)-Lim(var1)+1,1);
  end
  for var1 = 4:4:sizelim-1,
    if var1 < 8,
    R_FXG(Lim(var1)+0:Lim(var1+1)-1) = (-R_FX(Lim(var1)+0,1))*([Lim(var1)+0:Lim(var1+1)-1]-Lim(var1+1)) / (Lim(var1+1)-(Lim(var1)+0)); 
    R_FXD(Lim(var1)+0:Lim(var1+1)-1) = (R_FX(Lim(var1+1),1))*([Lim(var1)+0:Lim(var1+1)-1]-(Lim(var1)+0)) / (Lim(var1+1)-(Lim(var1)+0));
    end
  end 

end 

sizeR_FXG = size(R_FXG,1);
sizeR_FXD = sizeR_FXG; 


% % % % % % L'imprimer les figures

figure

% Single et double stance phase
hauteur_fleche = max(Fx_gauche(Lim(2)-20:Lim(3)+20))+3;
hauteur_texte = max(Fx_gauche(Lim(2)-20:Lim(3)+20))+5;
distance_fleche = 0.9;

plot((Lim(1)+distance_fleche)/frequency,hauteur_fleche,'color','k','Marker','<','MarkerFaceColor','black')
hold on
plot((Lim(2)-distance_fleche)/frequency,hauteur_fleche,'color','k','Marker','>','MarkerFaceColor','black')

line([(Lim(2)-10)/frequency,(Lim(3)+5)/frequency],[hauteur_fleche,hauteur_fleche],'color','k')
line([(Lim(2)-20)/frequency,(Lim(2)-10)/frequency],[hauteur_fleche,hauteur_fleche],'color','k','LineStyle','- -')
line([(Lim(3)+5)/frequency,(Lim(3)+20)/frequency],[hauteur_fleche,hauteur_fleche],'color','k','LineStyle','- -')



text((Lim(2)/100)-0.18,hauteur_texte,sprintf('SS1'),'FontName','arial','Fontsize',9);
plot((Lim(2)+distance_fleche)/frequency,hauteur_fleche,'color','k','Marker','<','MarkerFaceColor','black')
plot((Lim(3)-distance_fleche)/frequency,hauteur_fleche,'color','k','Marker','>','MarkerFaceColor','black')



text((Lim(2)+Lim(3)-4)/200,hauteur_texte,sprintf('DS'),'FontName','arial','Fontsize',9);
plot((Lim(3)+distance_fleche)/frequency,hauteur_fleche,'color','k','Marker','<','MarkerFaceColor','black')
plot((Lim(4)-distance_fleche)/frequency,hauteur_fleche,'color','k','Marker','>','MarkerFaceColor','black')

text((Lim(3)/100)+0.02,hauteur_texte,sprintf('SS2'),'FontName','Arial','Fontsize',9);
plot((Lim(4)+distance_fleche)/frequency,hauteur_fleche,'color','k','Marker','<','MarkerFaceColor','black')

for var1 = 2:4:sizelim-1,
%     plot([Lim(var1):Lim(var1+1)]/frequency,R_FX([Lim(var1):Lim(var1+1)]),'k:');
    line([Lim(var1)/frequency,Lim(var1)/frequency],[-2000,2000],'color','k','LineStyle',':');
end
for var1 = 3:4:sizelim,
    line([Lim(var1)/frequency,Lim(var1)/frequency],[-2000,2000],'color','k','LineStyle',':');
end

% % % Affichage des forces enregistrés par platformes des forces (Double stnace et 20 points de single stance)

for var1 = Lim(2)-20:1:Lim(3)-1, 
    line([var1,var1+1]/frequency,[Fx_gauche(var1),Fx_gauche(var1+1)],'color','k','LineStyle','-','Linewidth',1.35); %,'Linewidth',5);
    hold on;
end

for var1 = Lim(3):1:Lim(3)+20, 
     line([var1,var1+1]/frequency,[Fx_gauche(var1),Fx_gauche(var1+1)],'color','k','LineStyle','-','Linewidth',1.35);
    hold on;
end


% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Décomposition des forces
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
%-----------------------------------------------------
% Interpolation Spline :
%-----------------------------------------------------
tic
% C.I. : 
spline_x1 = Lim(1):Lim(2);
spline_x3 = Lim(3):Lim(4);

spline_xa = [spline_x1, spline_x3];

spline_R_FXD = zeros(1,echantillons);
spline_R_FXG = zeros(1,echantillons);

clear spli_R_FXG
for var1 = 1:size(spline_xa,2),
    spli_R_FXG(var1) = R_FXG(spline_xa(var1));
end

spline_R_FXG(Lim(1):Lim(1+3)) = spline(spline_xa,spli_R_FXG,Lim(1):Lim(1+3));
spline_R_FXD(Lim(1):Lim(1+3)) = R_FX(Lim(1):Lim(1+3))'- spline_R_FXG(Lim(1):Lim(1+3));

% affichage interpolation spline :
for var1 = Lim(2):1:Lim(3)-1, 
    line([var1,var1+1]/frequency,[spline_R_FXG(var1),spline_R_FXG(var1+1)],'color','m','LineStyle',':','Linewidth',1.2);
    hold on;
end

%-----------------------------------------------------
% Calculs erreurs et covariances d'interpolation Spline
% pendant phases de double appui :
%-----------------------------------------------------

if PiedDroit == 1, % choisir bonne courbe !
    spline_R_FXG2 = spline_R_FXD;
    spline_R_FXD2 = spline_R_FXG;
    spline_R_FXD = spline_R_FXD2;
    spline_R_FXG = spline_R_FXG2;
end
% Lim(1:7) = Lim(2:8);
    transition = 1;
    % Transition 1 :
    errorabs_G1 = abs (Fx_gauche(Lim(2):Lim(3)) - spline_R_FXG(Lim(2):Lim(3))');
    max_force1 = max( Fx_gauche(Lim(2):Lim(3)) + Fx_droite(Lim(2):Lim(3))) ; % JMNI
    error_G1 = abs ( errorabs_G1 / max_force1 ); % JMNI
    mean_errorG1 = mean(error_G1);
    cov_errorG1 = cov(error_G1);
    errorabs_D1 = abs (Fx_droite(Lim(2):Lim(3)) - spline_R_FXD(Lim(2):Lim(3))');
    error_D1 = abs ( errorabs_D1 / max_force1 ); % JMNI
    mean_errorD1 = mean(error_D1);
    cov_errorD1 = cov(error_D1);
    erreur_moy(transition) = mean_errorG1;
    erreur_moy(transition+1) = mean_errorD1;
    erreur_cov(transition) = cov_errorG1;
    erreur_cov(transition+1) = cov_errorD1;
    transition = transition + 2;
    
    max_forceX=(max_force1); %+max_force2)/2;
    errorabs_X_Spline = (mean(errorabs_G1)+mean(errorabs_D1))/2; %+mean(errorabs_G2)+mean(errorabs_D2))/4; % err. totale approx.
    errorrel_X_Spline = errorabs_X_Spline/max_forceX*100 % err. rel.
toc

% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Décomposition des forces
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
%-----------------------------------------------------
% Polynome d'ordre 3 : f(x) = Ax^3 + bx^2+ cx + d
%-----------------------------------------------------
tic

parD1 = Fx_gauche(Lim(2)); % Constraint 1: 1st point of Poly3_R_FXG is last point of the SS1
ndiff = 5; % number of points to calculate the slope on the ndiff last points of the SS1
parC1 = (Fx_gauche(Lim(2))-Fx_gauche(Lim(2)-ndiff-10))/(ndiff); % Constraint 2: slope of the 1st point of Poly3_R_FXG is slope of ndiff (5) last points of the SS1
LimEnd = Lim(3)-Lim(2); % Number of points of the DS1
parE1 = 0; % Arbitrary slope at the end of DS1
parB1 = -(3*parD1+2*parC1*LimEnd+parE1*LimEnd)/LimEnd^2; % Constraint 3 and 4: last point of the DS1 is at LimEnd and the slope at this last point is arbitrary (parE1), close to 0
parA1 = -(parD1+parC1*LimEnd+parB1*LimEnd^2)/LimEnd^3; % Constraint 3 and 4

Poly3_R_FXD = spline_R_FXD; % Initial condition. Poly3_R_FXD will be adapted further in the DS1, the rest being like spline_R_FXD.
Poly3_R_FXG = spline_R_FXG;

for Xpoly=1:LimEnd-10, % Creates the Poly3 during the DS1
    Poly3_R_FXG(Lim(2)+Xpoly) = parA1*Xpoly^3 + parB1 * Xpoly^2 + parC1 * Xpoly + parD1;
end
Poly3_R_FXD(Lim(2)+1:LimEnd-2) = R_FX(Lim(2)+1:LimEnd-2)'- Poly3_R_FXG(Lim(2)+1:LimEnd-2); % Force on other foot = Ftot - force on this foot.

% affichage Poly3 :
for var1 = Lim(2):1:Lim(3)-1, 
    line([var1,var1+1]/frequency,[Poly3_R_FXG(var1),Poly3_R_FXG(var1+1)],'color','r','LineStyle',':','Linewidth',1.2);
    hold on;
end

spline_R_FXG = Poly3_R_FXG;
spline_R_FXD = Poly3_R_FXD;
%-----------------------------------------------------
% Calculs erreurs et covariances de Polynomial d'ordre 3
% pendant phases de double appui :
%-----------------------------------------------------

if PiedDroit == 1, % choisir bonne courbe !
    spline_R_FXG2 = spline_R_FXD;
    spline_R_FXD2 = spline_R_FXG;
    spline_R_FXD = spline_R_FXD2;
    spline_R_FXG = spline_R_FXG2;
end

    transition = 1;
 
    errorabs_G1 = abs (Fx_gauche(Lim(2):Lim(3)) - spline_R_FXG(Lim(2):Lim(3))');
    max_force1 = max( Fx_gauche(Lim(2):Lim(3)) + Fx_droite(Lim(2):Lim(3))) ; % JMNI
    error_G1 = abs ( errorabs_G1 / max_force1 ); % JMNI
    mean_errorG1 = mean(error_G1);
    cov_errorG1 = cov(error_G1);
    errorabs_D1 = abs (Fx_droite(Lim(2):Lim(3)) - spline_R_FXD(Lim(2):Lim(3))');
    error_D1 = abs ( errorabs_D1 / max_force1 ); % JMNI
    mean_errorD1 = mean(error_D1);
    cov_errorD1 = cov(error_D1);
    erreur_moy(transition) = mean_errorG1;
    erreur_moy(transition+1) = mean_errorD1;
    erreur_cov(transition) = cov_errorG1;
    erreur_cov(transition+1) = cov_errorD1;
    transition = transition + 2;
    
    max_forceX=(max_force1); 
    errorabs_X_Poly3 = (mean(errorabs_G1)+mean(errorabs_D1))/2; 
    errorrel_X_Poly3 = errorabs_X_Poly3/max_forceX*100 % err. rel.
toc

% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Décomposition des forces
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
%-----------------------------------------------------
% Sinus : f(x) = a * sin ((x - c)/b) + d 
%-----------------------------------------------------    
tic
FSS=Fx_gauche(Lim(2));          %Last value of GRF in single stance before doublestance     
TDS=(Lim(3)-Lim(2))/100;        %Duration of double stance DS 


a=FSS*0.6;
b=2*TDS;
c=((Lim(2)-20)/100)-(TDS/2);
d=FSS*0.5; 


for var1 = Lim(2):1:Lim(3), 
    Sine_Fx_G(var1) = a*sin(((var1/100-c)*2*pi)/b)+d;      
end


Sine_Fx_D(Lim(2):1:Lim(3)) = R_FX(Lim(2):Lim(3))'- Sine_Fx_G(Lim(2):Lim(3));

% % % % % % % % % % % % 
% % % Affichage Sinus
% % % % % % % % % % % % 
for var1 = Lim(2):1:Lim(3)-1,
    line([var1,var1+1]/frequency,[Sine_Fx_G(var1),Sine_Fx_G(var1+1)],'color','[0 0.8 0]','LineStyle',':','Linewidth',1.2);
    hold on;
end

%-----------------------------------------------------
% Calculs erreurs et covariances de Sinus
% pendant phases de double appui :
%-----------------------------------------------------
transition=1;
 errorabs_G1 = abs (Fx_gauche(Lim(2):Lim(3)) - Sine_Fx_G(Lim(2):Lim(3))');
            max_force1 = max( Fx_gauche(Lim(2):Lim(3)) + Fx_droite(Lim(2):Lim(3))) ; % JMNI
            error_G1 = abs ( errorabs_G1 / max_force1 ); % JMNI
            mean_errorG1 = mean(error_G1);
            cov_errorG1 = cov(error_G1);
            errorabs_D1 = abs (Fx_droite(Lim(2):Lim(3)) - Sine_Fx_D(Lim(2):Lim(3))');
            error_D1 = abs ( errorabs_D1 / max_force1 ); % JMNI
            mean_errorD1 = mean(error_D1);
            cov_errorD1 = cov(error_D1);
            erreur_moy(transition) = mean_errorG1;
            erreur_moy(transition+1) = mean_errorD1;
            erreur_cov(transition) = cov_errorG1;
            erreur_cov(transition+1) = cov_errorD1;
            transition = transition + 2;
    
            max_forceX=(max_force1); %+max_force2)/2;
            errorabs_X_Sine = (mean(errorabs_G1)+mean(errorabs_D1))/2; %+mean(errorabs_G2)+mean(errorabs_D2))/4; % err. totale approx.
            errorrel_X_Sine = errorabs_X_Sine/max_forceX*100 % err. rel. 
toc
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Décomposition des forces
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
%-----------------------------------------------------
% Sigmoid f(x) = a / (1 + e ^( b (x - c)))
%-----------------------------------------------------
tic

FSS=Fx_gauche(Lim(2));
TDS=(Lim(3)-Lim(2))/100;
%%%% parameter starts%%%%
dFSS1=(Fx_gauche(Lim(2))-Fx_gauche(Lim(2)-1));


b=11.5;

vectordFSS=dFSS1;
for var1 = Lim(2):1:Lim(3), 
    Sine_Fx_G(var1) = FSS/(1+exp(b*(var1-Lim(2)-100*(0.5*TDS)+100*(-(dFSS1)/4))/100));  %                              
end
%%% parameter ends%%

Sine_Fx_D(Lim(2):1:Lim(3)) = R_FX(Lim(2):Lim(3))'- Sine_Fx_G(Lim(2):Lim(3));
   
% % % % % % % % % % % % % % % % % % % % % 
% % % Affichage Sigmoid function
% % % % % % % % % % % % % % % % % % % % %
for var1 = Lim(2):1:Lim(3)-1,
    line([var1,var1+1]/frequency,[Sine_Fx_G(var1),Sine_Fx_G(var1+1)],'color','[1 0.5 0]','LineStyle',':','Linewidth',1.2);
    hold on;
end


%-----------------------------------------------------
% Calculs erreurs et covariances de Sigmoid
% pendant phases de double appui :
%-----------------------------------------------------
transition=1;
 errorabs_G1 = abs (Fx_gauche(Lim(2):Lim(3)) - Sine_Fx_G(Lim(2):Lim(3))');
            max_force1 = max( Fx_gauche(Lim(2):Lim(3)) + Fx_droite(Lim(2):Lim(3))) ; % JMNI
            error_G1 = abs ( errorabs_G1 / max_force1 ); % JMNI
            mean_errorG1 = mean(error_G1);
            cov_errorG1 = cov(error_G1);
            errorabs_D1 = abs (Fx_droite(Lim(2):Lim(3)) - Sine_Fx_D(Lim(2):Lim(3))');
            error_D1 = abs ( errorabs_D1 / max_force1 ); % JMNI
            mean_errorD1 = mean(error_D1);
            cov_errorD1 = cov(error_D1);
            erreur_moy(transition) = mean_errorG1;
            erreur_moy(transition+1) = mean_errorD1;
            erreur_cov(transition) = cov_errorG1;
            erreur_cov(transition+1) = cov_errorD1;
            transition = transition + 2;
    
            max_forceX=(max_force1); 
            errorabs_X_Sigmoid = (mean(errorabs_G1)+mean(errorabs_D1))/2;
            errorrel_X_Sigmoid = errorabs_X_Sigmoid/max_forceX*100 % err. rel.

toc


            
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
% Décomposition des forces
% % % % % % % % % % % % % % % % % % % % % % % % % % % 
%-----------------------------------------------------
% Sinus + Sigmoid function
%-----------------------------------------------------           
tic
FSS=Fx_gauche(Lim(2));
TDS=(Lim(3)-Lim(2))/100;


%%%% parameter starts%%%%
dFSS1=(Fx_gauche(Lim(2))-Fx_gauche(Lim(2)-1));
dE=(Fx_gauche(Lim(3))-Fx_gauche(Lim(3)-3)/3);

b=11.5;

vectorb=-dFSS1/(3.5/1.7);vectordFSS=dFSS1;
flat=[];
for var1 = Lim(2):1:Lim(3), 
    Sine_Fx_G(var1) = FSS/(1+exp(b*(var1-Lim(2)-100*(0.5*TDS)+100*(-abs(dFSS1)/4))/100));  % 
    if (abs(Sine_Fx_G(var1)-Sine_Fx_G(var1-1))>(0.01))
        flat=[flat,var1];
    end
end

for var1=Lim(2):1:flat(length(flat)),
    Sine_Fx_G(var1)=dFSS1*((flat(length(flat))-Lim(2))/3)*sin(pi*(var1-Lim(2))/((flat(length(flat))-Lim(2))))+Sine_Fx_G(var1);
end
display (flat(length(flat))); display(Lim(2));
%%% parameter ends%%

Sine_Fx_D(Lim(2):1:Lim(3)) = R_FX(Lim(2):Lim(3))'- Sine_Fx_G(Lim(2):Lim(3));


% % % % % % % % % % % % % % % % % % % % % 
% % % Affichage Sinus + Sigmoid function
% % % % % % % % % % % % % % % % % % % % %
for var1 = Lim(2):1:Lim(3)-1,
    line([var1,var1+1]/frequency,[Sine_Fx_G(var1),Sine_Fx_G(var1+1)],'color','b','LineStyle','-','Linewidth',1.35);
    hold on;
end


%-----------------------------------------------------
% Calculs erreurs et covariances de Sinus + Sigmoid
% pendant phases de double appui :
%-----------------------------------------------------
transition=1;
 errorabs_G1 = abs (Fx_gauche(Lim(2):Lim(3)) - Sine_Fx_G(Lim(2):Lim(3))');
            max_force1 = max( Fx_gauche(Lim(2):Lim(3)) + Fx_droite(Lim(2):Lim(3))) ; % JMNI
            error_G1 = abs ( errorabs_G1 / max_force1 ); % JMNI
            mean_errorG1 = mean(error_G1);
            cov_errorG1 = cov(error_G1);
            errorabs_D1 = abs (Fx_droite(Lim(2):Lim(3)) - Sine_Fx_D(Lim(2):Lim(3))');
            error_D1 = abs ( errorabs_D1 / max_force1 ); % JMNI
            mean_errorD1 = mean(error_D1);
            cov_errorD1 = cov(error_D1);
            erreur_moy(transition) = mean_errorG1;
            erreur_moy(transition+1) = mean_errorD1;
            erreur_cov(transition) = cov_errorG1;
            erreur_cov(transition+1) = cov_errorD1;
            transition = transition + 2;
    
            max_forceX=(max_force1); 
            errorabs_X_Sine_Sigmoid = (mean(errorabs_G1)+mean(errorabs_D1))/2;
            errorrel_X_Sine_Sigmoid = errorabs_X_Sine_Sigmoid/max_forceX*100 % err. rel.

toc
% Write error
filename = 'Fx_Error.xlsx';
A = {'errorrel_X_Spline','errorrel_X_Poly3','errorrel_X_Sine','errorrel_X_Sigmoid','errorrel_X_Sine_Sigmoid'; errorrel_X_Spline, errorrel_X_Poly3,errorrel_X_Sine,errorrel_X_Sigmoid,errorrel_X_Sine_Sigmoid};
sheet = 1;
xlRange = 'A1';
xlswrite(filename,A,sheet,xlRange)

B = {'errorabs_X_Spline','errorabs_X_Poly3','errorabs_X_Sine','errorabs_X_Sigmoid','errorabs_X_Sine_Sigmoid'; errorabs_X_Spline, errorabs_X_Poly3,errorabs_X_Sine,errorabs_X_Sigmoid,errorabs_X_Sine_Sigmoid};
sheet = 2;
xlRange = 'A1';
xlswrite(filename,B,sheet,xlRange)

% % % Legend of the figure
title('Medio-Lateral forces (F_X)','FontName','arial','Fontsize',12)
xlabel('Time (s)','FontName','arial','Fontsize',10)
ylabel('Force (N)','FontName','arial','Fontsize',10)


axis([deb_aff fin_aff min(Fx_gauche(Lim(2)-20:Lim(3)+20))-10 max(Fx_gauche(Lim(2)-20:Lim(3)+20))+10]); %equal; 

h1 = plot(nan, nan, 'k-','Linewidth',1.35);
h2 = plot(nan, nan, 'm:','Linewidth',1.2);
h3 = plot(nan, nan,'r:','Linewidth',1.2); 
h4 = plot(nan, nan,'color',[0 0.8 0],'LineStyle',':','Linewidth',1.2); 
h5 = plot(nan, nan,'color',[1 0.5 0],'LineStyle',':','Linewidth',1.2);
h6 = plot(nan, nan, 'b-','Linewidth',1.35);



legend([h1, h2, h3, h4, h5 h6], 'Recorded GRF on the FPFs', 'Spline interpolation','3^r^d order polynomial function','Sine function','Sigmoid function','Sine-Sigmoid function');


