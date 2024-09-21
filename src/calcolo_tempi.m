function Deltat = calcolo_tempi(thetain,thetaf,a,e)
% ++descrizione++ 
% Questa funzione calcola il tempo necessario per spostarsi tra due punti
% (ad anomalie vere note) sulla stessa orbita.
%
% ++input++ 
% thetain[rad]: anomalia vera della posizione iniziale sull'orbita corrente
%
% thetaf[rad]: anomalia vera della posizione finale sull'orbita corrente
%
% a[km]: semiasse maggiore dell'orbita
%
% e[-]: eccentricità dell'orbita
%
% ++output++ 
% deltat[s]: tempo che occore per spostarsi dalla posizione a thetain a
% quella a thetaf

mu = 398600;

dt1 = theta_T(thetain,a,e,mu);
dt2 = theta_T(thetaf,a,e,mu);

periodo = 2 * pi * sqrt(a^3/mu);

if dt2 < dt1
    Deltat = periodo + dt2 - dt1;
else
    Deltat = dt2 - dt1;
end

