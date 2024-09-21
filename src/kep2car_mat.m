function [rv] = kep2car_mat(a,e,i,omega,o,theta,mu)
% ++descrizione++  
% Questa funzione, presi in input i parametri orbitali di un'orbita calcola 
% e resistuisce i parametri cartesiani (vettore posizione e velocità)
% riferiti ad un range di anomalie vere "theta" in input
%
% ++input++  
% -a[km]: semiasse maggiore dell'orbita
% -e[ - ]: eccentricità dell'orbita
% -i[rad]: inclinazione dell'orbita nel piano tridimensionale
% -omega[rad]: ascensione retta del nodo ascendente
% -o[rad]: anomalia del pericentro 
% -theta[rad]: vettore di anomalie vere raggiunte dallo spacraft in tempi
%  diversi
% -mu[km^3/s^2]: costante planetaria 
%
% ++output++ 
% -rv: matrice 6xn. Una singola colonna contiene nelle prime tre posizioni le
%  componenti x,y,z della posizione[km] e nelle tre seguenti le
%  componenti x,y,z della velocità[km/s]. I calcoli sono ripetuti ad ogni
%  valore del vettore theta (n: length(theta)) e riportati lungo le righe


p = a * (1 - e^2);

for j = 1:length(theta)

r = p/(1 + e * cos(theta(j)));
rr(1,j) = r * cos(theta(j));
rr(2,j) = r * sin(theta(j));
rr(3,j) = 0;

vv(1,j) =  - sqrt(mu/p) * sin(theta(j));
vv(2,j) = sqrt(mu/p) * (e + cos(theta(j)));
vv(3,j) = 0;

end

R3O = [cos(omega) sin(omega) 0;  -sin(omega) cos(omega) 0; 0 0 1];
r1i = [1 0 0; 0 cos(i) sin(i); 0  -sin(i) cos(i)];
r3o = [cos(o) sin(o) 0;  -sin(o) cos(o) 0; 0 0 1];

mr = r3o * r1i * R3O;
MR = mr';

for j = 1:length(theta)

rv(1:3,j) = MR * rr(:,j);
rv(4:6,j) = MR * vv(:,j);

end