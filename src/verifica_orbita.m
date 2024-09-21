function [kep,traj]=verifica_orbita(rv_in,dv,theta,kepE1,kepE2,par)
% ++Descrizione++
%
% Verifica la correttezza dei Dv calcolati.
% Date le coordinate e le componenti di velocità iniziali, la funzione
% calcola di volta in volta l'orbita su cui si trova il satellite 
% sommando alla velocità del punto caratterizzato dall'anomalia vera 
% specificata il delta V in forma vettoriale.
%
% ++Input++
%
% rv_in[km][km/s]: vettore con posizione e velocità del punto iniziale 
% 
% dv[km/s]: matrice 3xn avente come righe le componenti cartesiane x,y,z degli impulsi 
%   dati per effettuare le manovre. Ogni colonna equivale a un impulso diverso
% 
% theta[rad]: matrice nx2 con angoli in cui si effettuano le manovre: theta(1,1)
%   angolo di partenza, theta(1,2) angolo di arrivo dove si effettua la
%   manovra, theta(2,1) angolo di partenza per la nuova orbita (e spesso è
%   uguale a quello precedente) eccetera fino a theta(n,2) che è l'angolo
%   finale sull'orbita finale
%
% par=1 per plottare, par=0 per verificare solo i parametri
%
% kepE1[km] [rad]  parametri orbitali dell'orbita iniziale 
% kepE2 [km] [rad] parametri orbitali dell'orbita  finale
% 
% ++output++
%
% kep: vettore con i parametri orbitali dell'orbita finale (da usare come 
%   controllo) 
% traj: matrice di tutta la traiettoria percorsa dal satellite
%   con componenti raggio e velocità in vettore colonna, le colonne sono i tempi
% 
mu=398600;

if par==1
    tol=0.25*deg2rad(1);
else
    tol=1e-6;
end

l=size(theta);
[kep(1),kep(2),kep(3),kep(4),kep(5),kep(6)]=car2kep(rv_in,mu);
pos=plotOrbit_tiniz(kep,mu,theta(1,2),kep(6),tol);
traj=pos;

for i=2:l(1)
    v_new=traj(4:6,end)+dv(:,i-1);
    rv = [traj(1:3,end);v_new];
    [kep(1),kep(2),kep(3),kep(4),kep(5),kep(6)]=car2kep(rv,mu);
    pos = plotOrbit_tiniz(kep,mu,theta(i,2),theta(i,1),tol);
    traj = [traj,pos];
end

[af,ef,i_f,Omegaf,of,thetaf]=car2kep(traj(:,end),mu);
kep=[af,ef,i_f,Omegaf,of,thetaf];

if par==1

figure, hold on
hold on
[X1,Y1,Z1]=plotOrbit(kepE1,mu);
plot3(X1,Y1,Z1);
hold on
[X2,Y2,Z2]=plotOrbit(kepE2,mu);
plot3(X2,Y2,Z2);
hold on
Terra3d
for j=1:length(traj)
    plot3(traj(1,j),traj(2,j),traj(3,j),'.','Color','r')
    drawnow
end

end



