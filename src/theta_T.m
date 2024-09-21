function deltaT=theta_T(theta,a,e,mu)
% ++Descrizione++
% 
% Dati semiassemaggiore, eccentricità di un orbita e un'anomalia restituisce 
% il tempo necessario per raggiungere tale anomalia dal perigeo
%
% ++input++
% theta [rad]: anomalia di interesse
% a[km]: semiassemaggiore
% e[-]: eccentricità
% mu[km^3/s^2]: costante planetaria
% 
% ++output++
% deltaT[s]: tempo necessario per raggiungere l'anomalia theta
%

E = 2*atan(sqrt((1-e)/(1+e))*tan(theta/2));

if E<0
    E = E+2*pi;
end

deltaT = sqrt(a^3/mu)*(E-e*sin(E));