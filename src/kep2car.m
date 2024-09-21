function [rv] = kep2car(a,e,i,omega,o,theta,mu)
% ++descrizione++  
% Questa funzione, presi in input i parametri orbitali di un'orbita calcola 
% e resistuisce i parametri cartesiani (vettore posizione e velocità)
% riferiti all'anomalia vera "theta" in input
%
% ++input++  
% -a[km]: semiasse maggiore dell'orbita
% -e[ - ]: eccentricità dell'orbita
% -i[rad]: inclinazione dell'orbita nel piano tridimensionale
% -omega[rad]: ascensione retta del nodo ascendente
% -o[rad]: anomalia del pericentro 
% -theta[rad]: anomalia vera relativa alla posizione data in input
% -mu[km^3/s^2]: costante planetaria 
%
% ++output++ 
% -rv: vettore colonna contenente nelle prime tre posizioni le
%  componenti x,y,z della posizione[km] e nelle seguenti tre posizioni le
%  componenti x,y,z della velocità[km/s]


p = a * (1 - e^2);

r = p/(1 + e * cos(theta));
rr(1,1) = r * cos(theta);
rr(2,1) = r * sin(theta);
rr(3,1) = 0;

vv(1,1) =  - sqrt(mu/p) * sin(theta);
vv(2,1) = sqrt(mu/p) * (e + cos(theta));
vv(3,1) = 0;

R3O = [cos(omega) sin(omega) 0;  -sin(omega) cos(omega) 0; 0 0 1];
r1i = [1 0 0; 0 cos(i) sin(i); 0  -sin(i) cos(i)];
r3o = [cos(o) sin(o) 0;  -sin(o) cos(o) 0; 0 0 1];

mr = r3o * r1i * R3O;
MR = mr';

rv(1,1:3) = MR * rr;
rv(1,4:6) = MR * vv;

