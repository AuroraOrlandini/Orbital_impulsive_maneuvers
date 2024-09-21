% PROVA FINALE
% FONTANA, GERARDINI, ORLANDINI

% NOTA:le funzioni di verifica "verifica_orbita" sono tutte commentate onde
% evitare di appesantire ogni volta la compilazione del file.
% Eliminando il commento viene avviata la verifica grafica. Per ottenere la
% verifica dei parametri (senza plot) è necessario cambiare l'ultimo valore
% di input da 1 a 0

% dati globali
mu=398600;

% orbita iniziale (1)

rx=-7203.4669;
ry=-4889.1562;
rz=-503.0398;

vx=1.8390;
vy=-3.1230;
vz=-5.4180;

dati1=[rx; ry; rz; vx; vy; vz];

[a1,e1,i1,Omega1,o1,theta1]=car2kep(dati1,mu);
kepE1=[a1,e1,i1,Omega1,o1,theta1];
h1_vett=cross(dati1(1:3),dati1(4:6));
e1_vett=(cross(dati1(4:6),h1_vett))/mu-dati1(1:3)/norm(dati1(1:3));

% orbita finale (2)

a2=10810.0000;
e2=0.2329;
i2=1.5080;
Omega2=1.5980;
o2=0.4748;
theta2=1.0100;
kepE2=[a2,e2,i2,Omega2,o2,theta2];

dati2=kep2car(a2,e2,i2,Omega2,o2,theta2,mu);
h2_vett=cross(dati2(1:3),dati2(4:6));
e2_vett=(cross(dati2(4:6),h2_vett))/mu-dati2(1:3)/norm(dati2(1:3));

% plottaggio orbite

Terra3d;
[X1,Y1,Z1]=plotOrbit(kepE1,mu);
plot3(X1,Y1,Z1);
h = plot3(nan,nan,nan,'or');
step_animation = 10;
grid on
view([1 1 1]);
hold on
[X2,Y2,Z2]=plotOrbit(kepE2,mu);
plot3(X2,Y2,Z2);
hold on
plot3(dati2(1),dati2(2),dati2(3),'o'); % punto di inizio
hold on
plot3(dati1(1),dati1(2),dati1(3),'o') % punto finale


% dati generali

p1=a1*(1-e1^2);
va1=sqrt(mu/p1)*(1-e1);
vp1=sqrt(mu/p1)*(1+e1);
ra1=p1/(1-e1);
rp1=p1/(1+e1);

p2=a2*(1-e2^2);
va2=sqrt(mu/p2)*(1-e2);
vp2=sqrt(mu/p2)*(1+e2);
ra2=p2/(1-e2);
rp2=p2/(1+e2);

%% Trasferimento 1: STANDARD
% cambio piano orbitale, anomalia di pericentro, forma (velocità)

% cambio piano
[Deltav_1_1,thetaint1,u1,omegatr1]=cambio_piano(i1,i2,Omega1,Omega2,o1,p1,e1,a1,theta1,0); 

% sistemo anomalia pericentro
[Deltav_1_2,thetamanovra1_1,thetamanovra2_1]=manovra_secante(omegatr1,o2,p1,e1,thetaint1); 

% cambio forma
[Deltav_1_3,dT_1_4,~,~,Deltav_1_3a,Deltav_1_3b]=changeOrbitShape(a1,e1,a2,e2,0); 
aT_1=(rp1+ra2)/2;
eT_1=(-rp1+ra2)/(rp1+ra2);


% verifica

rv_1_1=kep2car(a1,e1,i1,Omega1,o1,thetaint1,mu);
rv_1_1=rv_1_1';
rv_1_2=kep2car(a1,e1,i2,Omega2,omegatr1,thetaint1,mu);
rv_1_2=rv_1_2';
dv_1_1=rv_1_2(4:6)-rv_1_1(4:6);
rv_1_3=kep2car(a1,e1,i2,Omega2,omegatr1,thetamanovra1_1,mu);
rv_1_3=rv_1_3';
rv_1_4=kep2car(a1,e1,i2,Omega2,o2,thetamanovra2_1,mu);
rv_1_4=rv_1_4';
dv_1_2=rv_1_4(4:6)-rv_1_3(4:6);
rv_1_5=kep2car(a1,e1,i2,Omega2,o2,0,mu);
rv_1_5=rv_1_5';
rv_1_6=kep2car(aT_1,eT_1,i2,Omega2,o2,0,mu);
rv_1_6=rv_1_6';
dv_1_3=rv_1_6(4:6)-rv_1_5(4:6);
rv_1_7=kep2car(aT_1,eT_1,i2,Omega2,o2,pi,mu);
rv_1_7=rv_1_7';
rv_1_8=kep2car(a2,e2,i2,Omega2,o2,pi,mu);
rv_1_8=rv_1_8';
dv_1_4=rv_1_8(4:6)-rv_1_7(4:6);

DV_1=[dv_1_1,dv_1_2,dv_1_3,dv_1_4];
THETA_1=[theta1,thetaint1;thetaint1,thetamanovra1_1;thetamanovra2_1,0;0,pi;pi,theta2];
% kep1_ver=verifica_orbita(dati1,DV_1,THETA_1,kepE1,kepE2,1);

Deltavtot_1=norm(dv_1_1)+norm(dv_1_2)+norm(dv_1_3)+norm(dv_1_4);

% calcolo tempi

dT_1_1=calcolo_tempi(theta1,thetaint1,a1,e1); % tempo da immissione a cambio piano
dT_1_2=calcolo_tempi(thetaint1,thetamanovra1_1,a1,e1); % da cambio piano a cambio anomalia
dT_1_3=calcolo_tempi(thetamanovra2_1,0,a1,e1); % da cambio anomalia a pericentro 
% da thetamanovra1 dT_1_4
dT_1_5=calcolo_tempi(pi,theta2,a2,e2); % da apocentro2 a theta2

Deltat_tot_1=dT_1_1+dT_1_2+dT_1_3+dT_1_4+dT_1_5; % delta t tot

%% Trasferimento 2: STANDARD-VARIANTE
% cambio piano orbitale, anomalia pericentro, forma (tempo)

% cambio piano 
[Deltav_2_1,thetaint2,u2,omegatr2,dT_2_1]=cambio_piano(i1,i2,Omega1,Omega2,o1,p1,e1,a1,theta1,1);

% sistemo anomalia pericentro
[Deltav_2_2,thetamanovra1_2,thetamanovra2_2]=manovra_secante(omegatr2,o2,p1,e1,thetaint2);

% cambio forma
[Deltav_2_3,dT_2_4,~,~,Deltav_2_3a,Deltav_2_3b]=changeOrbitShape(a1,e1,a2,e2,1);

aT_2=(rp2+ra1)/2;
eT_2=(-rp2+ra1)/(rp2+ra1);


% verifica

rv_2_1=kep2car(a1,e1,i1,Omega1,o1,thetaint2,mu);
rv_2_1=rv_2_1';
rv_2_2=kep2car(a1,e1,i2,Omega2,omegatr2,thetaint2,mu);
rv_2_2=rv_2_2';
dv_2_1=rv_2_2(4:6)-rv_2_1(4:6);
rv_2_3=kep2car(a1,e1,i2,Omega2,omegatr2,thetamanovra1_2,mu);
rv_2_3=rv_2_3';
rv_2_4=kep2car(a1,e1,i2,Omega2,o2,thetamanovra2_2,mu);
rv_2_4=rv_2_4';
dv_2_2=rv_2_4(4:6)-rv_2_3(4:6);
rv_2_5=kep2car(a1,e1,i2,Omega2,o2,pi,mu);
rv_2_5=rv_2_5';
rv_2_6=kep2car(aT_2,eT_2,i2,Omega2,o2,pi,mu);
rv_2_6=rv_2_6';
dv_2_3=rv_2_6(4:6)-rv_2_5(4:6);
rv_2_7=kep2car(aT_2,eT_2,i2,Omega2,o2,0,mu);
rv_2_7=rv_2_7';
rv_2_8=kep2car(a2,e2,i2,Omega2,o2,0,mu);
rv_2_8=rv_2_8';
dv_2_4=rv_2_8(4:6)-rv_2_7(4:6);

DV_2=[dv_2_1,dv_2_2,dv_2_3,dv_2_4];
THETA_2=[theta1,thetaint2;thetaint2,thetamanovra1_2;thetamanovra2_2,pi;pi,0;0,theta2];
% kep2_ver=verifica_orbita(dati1,DV_2,THETA_2,kepE1,kepE2,1);

Deltavtot_2=norm(dv_2_1)+norm(dv_2_2)+norm(dv_2_3)+norm(dv_2_4);

% calcolo tempi

% dT_2_1 da immissione a cambio piano
dT_2_2=calcolo_tempi(thetaint2,thetamanovra1_2,a1,e1); % da thetaint2 a thetamanovra1_2 dT_2_2
dT_2_3=calcolo_tempi(thetamanovra2_2,pi,a1,e1); % da cambio anomalia a apocentro
% da thetamanovra2 dT_2_4
dT_2_5=calcolo_tempi(0,theta2,a2,e2);% da pericentro2 a theta2

Deltat_tot_2=dT_2_1+dT_2_2+dT_2_3+dT_2_4+dT_2_5;

%% Trasferimento 3: forma-piano da pericentro
% forma, cambio di piano orbitale, anomalia pericentro (velocità)

% cambio forma 
[Deltav_3_1,dT_3_2,thetaf3,~,Deltav_3_1a,Deltav_3_1b]=changeOrbitShape(a1,e1,a2,e2,0);

% cambio piano
[Deltav_3_2,thetaint3,u1_3,omegaint3,dT_3_3]=cambio_piano(i1,i2,Omega1,Omega2,o1,p2,e2,a2,thetaf3,0);

% cambio anomalia pericentro
[Deltav_3_3, thetamanovra1_3,thetamanovra2_3]=manovra_secante(omegaint3,o2,p2,e2,thetaint3);

aT_3=(rp1+ra2)/2;
eT_3=(ra2-rp1)/(ra2+rp1);


% verifica

rv_3_1=kep2car(a1,e1,i1,Omega1,o1,0,mu);
rv_3_1=rv_3_1';
rv_3_2=kep2car(aT_3,eT_3,i1,Omega1,o1,0,mu);
rv_3_2=rv_3_2';
dv_3_1=rv_3_2(4:6)-rv_3_1(4:6);
rv_3_3=kep2car(aT_3,eT_3,i1,Omega1,o1,pi,mu);
rv_3_3=rv_3_3';
rv_3_4=kep2car(a2,e2,i1,Omega1,o1,pi,mu);
rv_3_4=rv_3_4';
dv_3_2=rv_3_4(4:6)-rv_3_3(4:6);
rv_3_5=kep2car(a2,e2,i1,Omega1,o1,thetaint3,mu);
rv_3_5=rv_3_5';
rv_3_6=kep2car(a2,e2,i2,Omega2,omegaint3,thetaint3,mu);
rv_3_6=rv_3_6';
dv_3_3=rv_3_6(4:6)-rv_3_5(4:6);
rv_3_7=kep2car(a2,e2,i2,Omega2,omegaint3,thetamanovra1_3,mu);
rv_3_7=rv_3_7';
rv_3_8=kep2car(a2,e2,i2,Omega2,o2,thetamanovra2_3,mu);
rv_3_8=rv_3_8';
dv_3_4=rv_3_8(4:6)-rv_3_7(4:6);

DV_3=[dv_3_1,dv_3_2,dv_3_3,dv_3_4];
THETA_3=[theta1,0;0,pi;pi,thetaint3;thetaint3,thetamanovra1_3;thetamanovra2_3,theta2];
% kep3_ver=verifica_orbita(dati1,DV_3,THETA_3,kepE1,kepE2,1);

Deltavtot_3=norm(dv_3_1)+norm(dv_3_2)+norm(dv_3_3)+norm(dv_3_4);

% tempi

dT_3_1=calcolo_tempi(theta1,0,a1,e1);% da immissione a pericentro 1
% da pericentro 1 ad apocentro 2 dT_3_2
% da apocentro 2 a thetaint3 dT_3_3
dT_3_4=calcolo_tempi(thetaint3,thetamanovra1_3,a2,e2);% da thetaint3 a thetamanovra_1
dT_3_5=calcolo_tempi(thetamanovra2_3,theta2,a2,e2);% da thetamanovra_2 a theta2

Deltat_tot_3=dT_3_1+dT_3_2+dT_3_3+dT_3_4+dT_3_5;

%% Trasferimento 4: forma-piano da apocentro
% forma, cambio di piano orbitale, anomalia pericentro (tempi)

% cambio forma 
[Deltav_4_1,dT_4_2,thetaf4,kepet4_1t,Deltav_4_1a,Deltav_4_1b]=changeOrbitShape(a1,e1,a2,e2,1);
pt_4=kepet4_1t(1)*(1-kepet4_1t(2)^2);

% cambio piano
[Deltav_4_2,thetaint4,u4_1,omegaint4,dT_4_3]=cambio_piano(i1,i2,Omega1,Omega2,o1,p2,e2,a2,thetaf4,1);

% cambio anomalia pericentro
[Deltav_4_3, thetamanovra1_4,thetamanovra2_4]=manovra_secante(omegaint4,o2,p2,e2,thetaint4);

% verifica

rv_4_1=kep2car(a1,e1,i1,Omega1,o1,pi,mu);
rv_4_2=kep2car(kepet4_1t(1),kepet4_1t(2),i1,Omega1,o1,pi,mu);
dv_4_1=rv_4_2(4:6)-rv_4_1(4:6);
dv_4_1=dv_4_1';
rv_4_3=kep2car(kepet4_1t(1),kepet4_1t(2),i1,Omega1,o1,0,mu);
rv_4_4=kep2car(a2,e2,i1,Omega1,o1,0,mu);
dv_4_2=rv_4_4(4:6)-rv_4_3(4:6);
dv_4_2=dv_4_2';
rv_4_5=kep2car(a2,e2,i1,Omega1,o1,thetaint4,mu);
rv_4_6=kep2car(a2,e2,i2,Omega2,omegaint4,thetaint4,mu);
dv_4_3=rv_4_6(4:6)-rv_4_5(4:6);
dv_4_3=dv_4_3';
rv_4_7=kep2car(a2,e2,i2,Omega2,omegaint4,thetamanovra1_4,mu);
rv_4_8=kep2car(a2,e2,i2,Omega2,o2,thetamanovra2_4,mu);
dv_4_4=rv_4_8(4:6)-rv_4_7(4:6);
dv_4_4=dv_4_4';

DV_4=[dv_4_1,dv_4_2,dv_4_3,dv_4_4];
THETA_4=[theta1,pi; pi,0; 0,thetaint4; thetaint4,thetamanovra1_4; thetamanovra2_4,theta2];
% kep4_ver=verifica_orbita(dati1,DV_4,THETA_4,kepE1,kepE2,1);

Deltavtot_4=norm(dv_4_1)+norm(dv_4_2)+norm(dv_4_3)+norm(dv_4_4);

% tempi

dT_4_1=calcolo_tempi(theta1,pi,a1,e1);% da immissione ad apocentro 1
% da apocentro 1 a pericentro 2 dT_4_2
% da pericentro 2 a thetaint1 dT_1_3
dT_4_4=calcolo_tempi(thetaint4,thetamanovra1_4,a2,e2);% da thetaint4 a thetamanovra1_4
dT_4_5=calcolo_tempi(thetamanovra2_4,theta2,a2,e2);% da thetamanovra2_4 a theta2

Deltat_tot_4=dT_4_1+dT_4_2+dT_4_3+dT_4_4+dT_4_5;

%% Trasferimento 5: dopo cambio piano manovra diretta all'apocentro finale
% Il primo impulso effettua lo stesso cambio piano visto per i casi
% precedenti (minimizzante il tempo impiegato). 
% Il secondo deltav immette il satellite su un'orbita intermedia che ha lo
% stesso apocentro dell'orbita finale e interseca l'orbita iniziale ad
% un'anomalia vera specifica che permette di effettuare contemporaneamente
% anche il cambio di anomalia al pericentro. 
% L'ultimo impulso è dato dalla differenza tra le velocità vettoriali
% all'apocentro delle orbite finale ed intermedia e permette di raggiungere
% l'orbita finale.

% considero lo stesso cambio di piano della manovra standard
[Deltav_1_1,thetaint1,u1,omegatr1,dT_1_2]=cambio_piano(i1,i2,Omega1,Omega2,o1,p1,e1,a1,theta1,0);

deltao_5=abs(o2-omegatr1);
thetaman_5=deltao_5;

% caratterizzo orbita di trasferimento 
rat_5=ra2;
rpt_5=p1/(1+e1*cos(2*pi-deltao_5));
et_5=(rat_5-rpt_5)/(rat_5+rpt_5);
at_5=(rat_5+rpt_5)/2;
pt_5=at_5*(1-et_5^2);
kepel_t_5=[at_5,et_5,i2,Omega2,o2,2*pi-deltao_5];


% verifica
rv_5_1=kep2car(a1,e1,i1,Omega1,o1,thetaint2,mu); % iniziale
rv_5_1=rv_5_1';
rv_5_2=kep2car(a1,e1,i2,Omega2,omegatr1,thetaint2,mu); % dopo cambio piano
rv_5_2=rv_5_2';
dv_5_1=rv_5_2(4:6)-rv_5_1(4:6);
rv_5_3=kep2car(a1,e1,i2,Omega2,omegatr1,2*pi-thetaman_5,mu); % dopo cambio piano al punto di manovra
rv_5_3=rv_5_3';
rv_5_4=kep2car(at_5,et_5,i2,Omega2,o2,0,mu); % trasferimento al punto di manovra
rv_5_4=rv_5_4';
dv_5_2=rv_5_4(4:6)-rv_5_3(4:6);
rv_5_5=kep2car(at_5,et_5,i2,Omega2,o2,pi,mu);
rv_5_5=rv_5_5';
rv_5_6=kep2car(a2,e2,i2,Omega2,o2,pi,mu);
rv_5_6=rv_5_6';
dv_5_3=rv_5_6(4:6)-rv_5_5(4:6);

DV_5=[dv_5_1,dv_5_2,dv_5_3];
THETA_5=[theta1,thetaint2;thetaint2,2*pi-thetaman_5; 0,pi; pi,theta2];

% kep5_ver=verifica_orbita(dati1,DV_5,THETA_5,kepE1,kepE2,1);
Deltavtot_5=norm(dv_5_1)+norm(dv_5_2)+norm(dv_5_3);
Deltav_5_1=norm(dv_5_1);
Deltav_5_2=norm(dv_5_2);
Deltav_5_3=norm(dv_5_3);

% minimizzando tempi:
dT_5_1=dT_2_1;
dT_5_2=calcolo_tempi(thetaint2,2*pi-thetaman_5,a1,e1);
dT_5_3=pi*sqrt(at_5^3/mu);
dT_5_4=calcolo_tempi(pi,theta2,a2,e2);

Deltat_tot_5=dT_5_4+dT_5_3+dT_5_2+dT_5_1;

%% Trasferimento 6: dopo cambio piano manovra diretta al pericentro dell'orbita finale
% Il primo impulso effettua lo stesso cambio piano visto per i casi
% precedenti (minimizzante il tempo impiegato). 
% Il secondo deltav immette il satellite su un'orbita intermedia che ha lo
% stesso pericentro dell'orbita finale e interseca l'orbita iniziale ad
% un'anomalia vera specifica che permette di effettuare contemporaneamente
% anche il cambio di anomalia al pericentro. 
% L'ultimo impulso è dato dalla differenza tra le velocità vettoriali
% al pericentro delle orbite finale ed intermedia e permette di raggiungere
% l'orbita finale.

% considero lo stesso cambio di piano della manovra standard
[Deltav_1_1,thetaint1,u1,omegatr1,dT_1_2]=cambio_piano(i1,i2,Omega1,Omega2,o1,p1,e1,a1,theta1,0);

deltao_6=abs(o2-omegatr1);
thetaman_6=pi-(deltao_6);

% caratterizzo orbita di trasferimento
rpt_6=rp2;
rat_6=p1/(1+e1*cos(thetaman_6));
at_6=(rat_6+rpt_6)/2;
et_6=(rat_6-rpt_6)/(rat_6+rpt_6);
pt_6=at_6*(1-et_6^2);

rv_6_1=kep2car(a1,e1,i1,Omega1,o1,thetaint2,mu);
rv_6_1=rv_6_1';
rv_6_2=kep2car(a1,e1,i2,Omega2,omegatr1,thetaint2,mu);
rv_6_2=rv_6_2';
dv_6_1=rv_6_2(4:6)-rv_6_1(4:6);
rv_6_3=kep2car(a1,e1,i2,Omega2,omegatr1,thetaman_6,mu);
rv_6_3=rv_6_3';
rv_6_4=kep2car(at_6,et_6,i2,Omega2,o2,pi,mu);
rv_6_4=rv_6_4';
dv_6_2=rv_6_4(4:6)-rv_6_3(4:6);
rv_6_5=kep2car(at_6,et_6,i2,Omega2,o2,0,mu);
rv_6_5=rv_6_5';
rv_6_6=kep2car(a2,e2,i2,Omega2,o2,0,mu);
rv_6_6=rv_6_6';
dv_6_3=rv_6_6(4:6)-rv_6_5(4:6);

DV_6=[dv_6_1,dv_6_2,dv_6_3];
THETA_6=[theta1,thetaint2;thetaint2,thetaman_6;pi,0;0,theta2];
%kep6_ver=verifica_orbita(dati1,DV_6,THETA_6,kepE1,kepE2,1);

Deltavtot_6=norm(dv_6_3)+norm(dv_6_2)+norm(dv_6_1);
Deltav_6_1=norm(dv_6_1);
Deltav_6_2=norm(dv_6_2);
Deltav_6_3=norm(dv_6_3);

% calcolo tempo
dT_6_1=dT_2_1;
dT_6_2=calcolo_tempi(thetaint2,thetaman_6,a1,e1);
dT_6_3=pi*sqrt(at_6^3/mu);
dT_6_4=calcolo_tempi(0,theta2,a2,e2);

Deltat_tot_6=dT_6_4+dT_6_3+dT_6_2+dT_6_1;

%% Trasferimento 7: secante-tangente al pericentro finale
% In questa manovra usiamo due singoli impulsi per arrivare all'orbita
% finale. La prima variazione di velocità consente di immettere il
% satellite su un'orbita intermedia che ha la stessa anomalia al pericentro
% dell'orbita finale ed interseca l'orbita iniziale ad un'anomalia vera
% specifica che consente di effettuare contemporaneamente anche un cambio
% piano. Il secondo impulso è dato dalla differenza delle due velocità
% vettoriali al pericentro (finale ed intermedia) e immette il satellite
% sull'orbita finale corretta, al suo pericentro.

% cambio piano 
[dv7,thetaman_7,u1,otrasf7,dT_7_1,u2]=cambio_piano(i1,i2,Omega1,Omega2,o1,p1,e1,a1,theta1,1);

deltao_7=abs(o2-o1);
theta_tra_7_int=u2-o2+pi;

% caratterizzo orbita di trasferimento
rpt_7=rp2;
r1_int=p1/(1+e1*cos(thetaman_7));
pt_7=r1_int*(1-cos(theta_tra_7_int))*(1-r1_int/rp2*cos(theta_tra_7_int))^(-1);
et_7=pt_7/rpt_7-1;
at_7=pt_7/(1-et_7^2);
rat_7=pt_7/(1-et_7);

rv_7_1=kep2car(a1,e1,i1,Omega1,o1,thetaman_7,mu);% velocità al punto di manovra sull'orbita 1
rv_7_1=rv_7_1';
rv_7_2=kep2car(at_7,et_7,i2,Omega2,o2,theta_tra_7_int,mu);% velocità al punto di manovra sull'orbita di trasferimento
rv_7_2=rv_7_2';
dv_7_1=rv_7_2(4:6)-rv_7_1(4:6);
rv_7_3=kep2car(at_7,et_7,i2,Omega2,o2,0,mu);
rv_7_4=kep2car(a2,e2,i2,Omega2,o2,0,mu);
rv_7_4=rv_7_4';
rv_7_3=rv_7_3';
dv_7_2=rv_7_4(4:6)-rv_7_3(4:6);

THETA_7=[theta1, thetaman_7; theta_tra_7_int,0; 0,theta2]; 
DV_7=[dv_7_1, dv_7_2];
% kep_7_ver=verifica_orbita(dati1,DV_7,THETA_7,kepE1,kepE2,1);

Deltavtot_7=norm(dv_7_1)+norm(dv_7_2);
Deltav_7_1=norm(dv_7_1);
Deltav_7_2=norm(dv_7_2);


% tempi
% dT_7_1
dT_7_2=calcolo_tempi(theta_tra_7_int,0,at_7,et_7);
dT_7_3=calcolo_tempi(0,theta2,a2,e2);
Deltat_tot_7=dT_7_3+dT_7_2+dT_7_1;

%% Trasferimento 8: da nodo ascendento orbita iniziale a punto finale
% In questa manovra usiamo due singoli impulsi per arrivare all'orbita
% finale. Il primo delta-v trasferisce il satellite dal nodo ascendente 
% dell'orbita iniziale ad un'anomalia vera corrispondente sull'orbita
% intermedia di trasferimento. La seconda variazione di velocità
% trasferisce il satellite direttamente al punto di anomalia vera finale
% sull'orbita finale, dove si impone l'intersezione con l'orbita di
% trasferimento. 
% L'orbita di trasferimento è scelta in modo da 
% intersecare l'orbita iniziale al suo nodo ascendente e l'orbita finale
% al punto di anomalia vera finale. La scelta dell'orbita è fatta in modo
% da minimizzare il tempo speso totale al variare dell'anomalia al
% pericentro tra 0° e 360°.

rv1_8=kep2car(a1,e1,i1,Omega1,o1,2*pi-o1,mu); % punto corrispondente all'asse dei nodi dell'orbita 1 (nodo ascendente)
R1N_8=rv1_8(1:3);
R2F_8=dati2(1:3); % punto di arrivo
R_8=norm(R2F_8)/norm(R1N_8);

phi_8=acos(dot(R1N_8,R2F_8)/(norm(R1N_8)*norm(R2F_8))); % angolo tra i raggi del punto R1N e RF2 sull'orbita di trasferimento

beta_8=asin(sin(Omega2-Omega1)/sin(phi_8)*sin(i2)); 

gamma_8=asin(sin(o2+theta2)*sin(beta_8)/sin(Omega2-Omega1)); % inclinazione orbita di trasferimento

Omegat_8=Omega1; % ascensione retta del nodo ascendente dell'orbita di trasferimento pari a quella dell'orbita iniziale

ot_8=0:deg2rad(1):2*pi; % anomalia del pericentro

% ottimizzazione: per ogni valore di anomalia del pericentro calcolo le
% possibili orbite ellittiche passanti per i due punti scelti

for  i=1:length(ot_8)

thetaC1_8=2*pi-ot_8(i); % punto di partenza 

et_8=(1-R_8)/(R_8*cos(thetaC1_8+phi_8)-cos(thetaC1_8));

pt_8(i)=norm(R1N_8)*(1+et_8*cos(thetaC1_8));

at_8=pt_8(i)/(1-(et_8)^2);

kepEt_8(i,:)=[at_8,et_8,gamma_8,Omegat_8,ot_8(i),thetaC1_8];

rv_8_1=kep2car(a1,e1,i1,Omega1,o1,2*pi-o1,mu);
rv_8_2=kep2car(at_8,et_8,gamma_8,Omegat_8,ot_8(i),thetaC1_8,mu);
dv_8_1(:,i)=rv_8_2(4:6)-rv_8_1(4:6);
rv_8_3=kep2car(at_8,et_8,gamma_8,Omegat_8,ot_8(i),thetaC1_8+phi_8,mu);
rv_8_4=kep2car(a2,e2,i2,Omega2,o2,theta2,mu);
dv_8_2(:,i)=rv_8_4(4:6)-rv_8_3(4:6);

Deltavtot_8(i)=norm(dv_8_2(:,i))+norm(dv_8_1(:,i));

dT_8_1=calcolo_tempi(theta1,2*pi-o1,a1,e1);
dT_8_2=calcolo_tempi(thetaC1_8,phi_8+thetaC1_8,at_8,et_8);

mat_Deltat_tot_8(i)=dT_8_2+dT_8_1;

end

Valori=[kepEt_8,mat_Deltat_tot_8',Deltavtot_8',dv_8_1',dv_8_2',pt_8'];

j=1;

for i=1:length(ot_8)
    if Valori(i,2)>0 && Valori(i,2)<1 % controllo sull'eccentricità
        rp8(i)=Valori(i,end)/(1+Valori(i,2));
        if rp8(i)>=6471 % nel suo punto più vicino al fuoco l'orbita deve essere a una distanza sufficiente dalla terra
            Valori_ottimali(j,:)=Valori(i,:);
            j=j+1;
        end
    end
end

[Deltat_tot_8,temp]=min(Valori_ottimali(:,7));
Deltav_8 = Valori_ottimali(temp,8);


DV_8=[Valori_ottimali(temp,9:11)',Valori_ottimali(temp,12:14)'];
THETA_8=[theta1,2*pi-o1;Valori_ottimali(temp,6),phi_8+Valori_ottimali(temp,6);theta2,theta2];
% kep8_ver=verifica_orbita(dati1,DV_8,THETA_8,kepE1,kepE2,1);