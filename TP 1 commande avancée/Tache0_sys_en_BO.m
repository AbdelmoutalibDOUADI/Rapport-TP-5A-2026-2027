% Script de simulation du système pendulaire en boucle ouverte (Tâche 0)
clear all; close all; clc;

%% Paramètres du système
m0 = 0.25;      % masse (kg)
r = 1;          % longueur du bras (m)
k = 0.1;        % coefficient de frottement visqueux
g = 10;         % gravité (m/s²)
J = m0*(r^2);   % moment d'inertie

% Période d'échantillonnage 
Te = 0.05;  % Pas pour RK2

% Pas de calcul à utiliser dans la méthode de Runge-Kutta 
Tc = 0.005;
m=Te/Tc;

%% Paramètres de simulation
t = 0:Te:3000*Te;   % grille de temps de simulation (s)
n=length(t);   %  durée de simulation avec Te comme unité de temps  

%% Définition de l'équation différentielle
f = @(x, u, m0) [x(2); (u - k*x(2) - m0*g*r*sin(x(1)))/(m0*r^2)];

%% Amplitudes des échelons de commande à tester
Amplitudes=[0.1 0.5 1.0 1.5];   % (N.m)

for j=1:length(Amplitudes)

    %% signal de commande : échelon d'amplitude A0
    A0=Amplitudes(j);
    consigne= @(t) A0 * (t >= 0); % Échelon de commande A0 (N.m)

    %% Initialisations vecteur d'état, commande
    x=[0 0]';      %  initialisation du vecteur d'état [position, vitesse]'
    X=x;           % X sauvegarde des valeurs de x(k) à différents instants k
    u = 0;   % u commande à l'instant présent k
    u1=0;    % u1 commande à l'instant présent k-1
    U=u;     % Sauvegarde des valeurs de la commande u(k) 
    y=x(1);  % sortie (angle theta) 
    Y=y;       % sauvegarde des sorties y(k)

    %% Début de la boucle de simulation (unité temps horloge = Te)
       %% Acquisition de la mesure de sortie
       % Résolution de l'équation différentielle du pendule 
       % par la méthode de Runge-Kutta d'ordre 2 (méthode de Heun)
       % pas de calcul=Tc
     for kk=1:n-1
         for i = 1:m
            x1=x+Tc*f(x,u1, m0);
            x = x+(Tc/2)*( f(x1,u, m0) + f(x,u1, m0) );
         end

       % Lecture de sortie
        y=x(1);

       u1=u;  % sauvegarde de u(k) pour le futur
       u=consigne(t(kk));  % commande en BO : échelon

       % Sauvgardes
        Y=[Y y];
        X =[X x];    % sauvegarde
        U=[U u];     % sauvegarde
     end 

    %% Affichage des résultats (une figure par amplitude)
     figure (j)
     plot(t,Y,'b-', t,U,'r-');
     title(sprintf('Réponse en BO à un échelon de %.1f N.m', A0));
     xlabel('Temps (s)');
     ylabel('Angle (°)');
     legend('sortie \theta','échelon de commande');
     xlim([0 40]);
     grid on;

end