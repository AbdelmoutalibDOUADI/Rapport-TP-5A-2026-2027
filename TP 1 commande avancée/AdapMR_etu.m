% Script de simulation du système pendulaire - Commande Adaptative à Modèle de Référence
clear all; close all; clc;

%% ===== Choix de l'essai =====================================
% Tâche 16 : mode_BO=1, A0=0.5, A1_fin=1,  m0_fin=0.25
% Tâche 18 : mode_BO=0, A0=0.1, A1_fin=50, m0_fin=0.25   (amplitude 0.1 -> 5)
% Tâche 19 : mode_BO=0, A0=0.1, A1_fin=20, m0_fin=0.25   (amplitude 0.1 -> 2.0)
% Tâche 20 : mode_BO=0, A0=0.1, A1_fin=1,  m0_fin=0.5    (masse 0.25 -> 0.5)
mode_BO = 1;        % 1 = test de l'estimateur en boucle ouverte (u = consigne)
A0      = 0.5;      % amplitude de départ de la consigne
A1_fin  = 1;        % gain d'amplification de la consigne à t = 100 s
m0_fin  = 0.25;     % masse après t = 100 s
c       = 0;        % polynôme C(q^-1) = 1 + c*q^-1   (essayer 0 puis -0.5)
prefixe = 'CAMR_tache16';   % nom de base des images sauvegardées
%% ============================================================

%% Paramètres du système
m0 = 0.25;      % masse (kg)
r = 1;          % longueur du bras (m)
K = 0.1;        % coefficient de frottement visqueux
g = 10;         % gravité (m/s²)
J = m0*(r^2);   % moment d'inertie

%%  Période d'échantillonnage et Pas de calcul de la méthode de RK
Te = 0.05;  % période d'échantillonnage
Tc = 0.01;  % Pas de calcul
m=Te/Tc;


%% Paramètres de simulation
t = 0:Te:4000*Te;   % grille de temps de simulation (s)  -> [0, 200 s]
n=length(t);   %  durée de simulation avec Tc comme unité de temps  

%% Fonction utilisée dans la simulation du pendule
f = @(x, u,m0) [x(2); (u - K*x(2) - m0*g*r*sin(x(1)))/(m0*r^2)];  %modèle NL
% f = @(x, u,m0) [x(2); (u - K*x(2) - m0*g*r*x(1))/(m0*r^2)];   % modèle linéaire


%% Initialisations des signaux du régulateur et vecteurs/matrices de sauvegarde
theta0=1; thetadot0=-1;  %initialisation position et vitesse initiales
x=[theta0 thetadot0]';      %  initialisation du vecteur d'état [position, vitesse]'
y=x(1);       % denotera la sortie y(k) à l'instant présent k
y1=1;         %  y1 denotera y(k-1)
y2=1;         %  y2 denotera y(k-2)
X=x;           % X matrice desauvegarde des valeurs de x(k) à différents instants k
u = 0;   % u commande à l'instant présent k en boucle fermée
u1=0;    % u1 commande à l'instant présent k-1
u2=0;     % u2 commande à l'instant présent k-2
U=[u];     % Sauvegarde des valeurs de la commande u(k) 
I=0;       % integrale de u(k)
Y=[y];     % sauvegarde des sorties antérieures

 %% Début de calcul et discrétisation de la fonction de transfert (unité temps horloge = Te)
 %%
% coefficient fonction de transfert
a1c=K/(m0*(r^2));
a0c=g/r;
b0c=1/(m0*(r^2));

 % Fonctions de transfert continue du système
G=tf(b0c,[1 a1c a0c]);
%disp('FT continue G(p)=')
%G

 % Fonctions de transfert échantillonnée du système
Ge=c2d(G,Te,"zoh");
%disp('FT continue G(p)=')
%Ge

% Création de A(q^-1) et B(q^-1) tel que Ge(z)=(z^-d)*B(z^-1)/A(z^-1)
% Ici  d=1
% Récupération des coefficients B(z^-1) et A(z^-1)
B=Ge.num{1};
A=Ge.den{1};
b0=B(2);  b1=B(3);
a0=A(1); a1=A(2); a2=A(3);

%disp('B(q^-1)=')
B=tf([b0 b1],[1 0],Te,'variable','q^-1');

%disp('A(q^-1)=')
A=tf([a0 a1 a2],[1 0],Te,'variable','q^-1');

%% Formulation des performances désirées
%% Définition du modèle de référence
xi=1;  % coefficient d'amortissement
w0=1;   % pulsation propre

%disp('Modèle de référence continu Gm(p)=')
 Gm=tf(1, [1/w0^2  2*xi/w0  1]);

%disp('Modèle de référence discret  Gme(z)=')
Gme=c2d(Gm,Te,"zoh");

% signal de consigne
% consigne= @(t) 0.1 * (t >= 0); % Échelon de consigne A0 (N.m)
consigne = @(t) A0 * (mod(t, 30) < 15); % Consigne carrée: A0 pendant 15s, 0 pendant 15s

 % Définition du signal de référence

 %disp('Signal de référence yref')
 Consigne=[]; 
 Yref=[];   % sauvegarde signal référence
 y=theta0;  % initialisation sortie
 y1=theta0;  % initialisation sortie
 e=[y-0 y-0]';  % erreurs de poursuite initiales
 E=[e];      % sauvegarde erreurs de poursuite

 % changement de gain d'amplification à mi-parcours de la consigne
    A1=1;  % gain d'amplification de la consigne
    n1=n/2;   % instant de changement de l'amplitude de la consigne.
  for k=1:n
  consigne_k=A1*( 2*consigne(t(k))-A0 );

    Consigne=[Consigne consigne_k];
    
   % changement d'amplitude de la consigne
    if k==floor(n1)
        A1=A1_fin;  
    end
  end
Yref=lsim(Gme,Consigne,t);

% Création du polynôme C(q^-1)=1+c*q^-1
%disp('Polynomial C(z^-1)=')
C=tf([1 c],[1 0],Te);

%% Synthèse de régulateur à Modèle de Référence
%%
% Calcul des coefficients des  polynômes R(q^-1) et S(q^-1) solution de A*R+(q^-1)S=C
% Loi de commande R*B*u(t)+S*y(t)=C*yref(t+d)
%R=1;
 s0=c-a1;  s1=-a2;

%% Initialisation relatives à l'estimateur paramétrique
%%
Theta=[0.01; 0; 0; 0];   %vecteur des paramètres estiméz [b0 b1 s0 s1]
Thetav=[b0; b1; s0; s1];   %vecteur des vrais paramètres [b0 b1 s0 s1]
THETA=Theta;  THETAv=Thetav;
phi=[u1; u2; y1; y2];
epsilon=(y+c*y1)-phi'*Theta;
Eps=[epsilon];
p0=1000;
P=p0*eye(4);
vp=min(eig(P));
Vp=[vp];

%% Début de la boucle de commande à MR (unité temps horloge = Te)
   %% Acquisition de la mesure de sortie
   % Résolution de l'équation différentielle du pendule 
   % par la méthode de Runge-Kutta d'ordre 2 (méthode de Heun)
   % pas de calcul=Tc

  for k=1:n-1
        y2=y1;  y1=y;   % sauvegarde de y(k-1) et y(k-2)
        u2=u1;  % sauvegarde de u(k-2)
        u1=u;   % sauvegarde de u(k-1)

     for i = 1:m
        x1=x+Tc*f(x,u1,m0);
        x = x+(Tc/2)*( f(x1,u,m0) + f(x,u1,m0) );
      end
    y=x(1);  % sauvegarde de y(k)

    % Mise à jour du vecteur d'observation
    phi=[u1; u2; y1; y2];
   
   % mise à jour des paramètres du régulateur
    epsilon=(y+c*y1)-phi'*Theta;
    Theta = Theta + P*phi*epsilon/( 1 + phi'*P*phi );   %estimateur des MC  (Tâche 16)
    P     = P - (P*phi*phi'*P)/( 1 + phi'*P*phi );      %estimateur des MC  (Tâche 16)
    
    %sauvergarde des variables 
    THETA=[THETA Theta]; %sauvrgarde des paramètres estimés
    THETAv=[THETAv Thetav];   %sauvrgarde des vrais paramètres

     % Loi de commande adaptative à MR  
    beta0=Theta(1);
    if (abs(Theta(1))<0.001)
        beta0=(0.001)*sign(Theta(1));
    end

   if mode_BO==1
    u=Consigne(k);   % Tâche 16 : test de l'estimateur en boucle ouverte
   else
    u = -Theta(2)*u1 - Theta(3)*y - Theta(4)*y1 + Yref(k+1) + c*Yref(k);  % loi de commande  (Tâche 17)
    u = u/beta0;   % loi de commande suite
   end
 
   e=[y-Yref(k) y-Consigne(k)]';  % calcul des erreurs de poursuite
  X =[X x]; % sauvegarde du vecteur d'état
  Y=[Y y];  % sauvegarde du signal de sortie
  U=[U u];  % sauvegarde du signal de commande
  E=[E e];  % sauvegarde des erreurs de poursuite
  Eps=[Eps epsilon];  % sauvegarde des erreurs de poursuite
 
  % Changement de la valeur de la masse du système (sans toucher au régulateur)
       if k==floor(n1)
        m0=m0_fin; 
        % Calcul des nouveaux vrais paramètres du régulateur en vue de
        % comparaison avec les estimés
        a1c=K/(m0*(r^2)); a0c=g/r; b0c=1/(m0*(r^2));
        G=tf(b0c,[1 a1c a0c]);
        Ge=c2d(G,Te,"zoh");
        B=Ge.num{1};  A=Ge.den{1};
        b0=B(2);  b1=B(3);
        a0=A(1); a1=A(2); a2=A(3);
         s0=c-a1;  s1=-a2;
         Thetav=[b0; b1; s0; s1];   % mise à jour des vrais paramètres
 
         P=1000*eye(4);  % Reset de P
    end

 end 

% %% Affichage des résultats
 figure (1)
  plot(t,Yref,'b-',t,Consigne,'r--');
  title('Signal de consigne et signal de référence');
  xlabel('Temps (s)');
  ylabel('Angle (rad)');
  legend('référence y^*','consigne');
  grid on;

 figure (2)
  plot(t,Y,'b-', t,Yref,'r-',t, Consigne,'g--');
  title('Position pendule -signal référence yref');
  xlabel('Temps (s)');
  ylabel('Angle (rad)');
  legend('sortie y','référence y^*','consigne');
  grid on;

   figure (3)
   plot(t,U);
   title('Signal de commande (CAMR)');
  xlabel('Temps (s)');
  ylabel('commande (N*m)');
  legend('commande u');
  grid on;

   figure (4)
   plot(t,E(1,:),'r-',t,E(2,:),'g-');
   title('Erreurs de poursuite (CAMR)');
  xlabel('Temps (s)');
  ylabel('erreurs (rad)');
  legend('y - y^*','y - consigne');
  grid on;

   figure (5)
   plot(t,Eps,'b-');
   title('Erreur d estimation epsilon');
  xlabel('Temps (s)');
  ylabel('erreur d estimation');
  grid on;

   figure (6)
   plot(t,THETA(1,:),'r-',t,THETAv(1,:),'r--',t,THETA(2,:),'b-',t,THETAv(2,:),'b--');
   title('Paramètres b0 et b1 et leurs estimés');
  xlabel('Temps (s)');
  ylabel('paramètre b0 et b1');
  legend('b0 estimé','b0 vrai','b1 estimé','b1 vrai');
  grid on;

   figure (7)
   plot(t,THETA(3,:),'r-',t,THETAv(3,:),'r--',t,THETA(4,:),'b-',t,THETAv(4,:),'b--');
   title('Paramètres s0 et s1 et leurs estimés');
  xlabel('Temps (s)');
  ylabel('paramètres s0 et s1');
  legend('s0 estimé','s0 vrai','s1 estimé','s1 vrai');
  grid on;

%% Sauvegarde automatique des figures
dossier = fullfile(fileparts(mfilename('fullpath')),'figures');
if ~exist(dossier,'dir')
    mkdir(dossier);
end

figs = findobj('Type','figure');
for kf = 1:numel(figs)
    nomfig = sprintf('%s_fig%d.png', prefixe, figs(kf).Number);
    exportgraphics(figs(kf), fullfile(dossier,nomfig), 'Resolution', 300);
end
disp(['Figures enregistrees dans : ' dossier]);