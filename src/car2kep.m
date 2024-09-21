function [a,e,i,omega,o,theta] = car2kep(s,mu)
% ++descrizione++ 
% Questa funzione, presi in input la posizione e la velocità vettoriali (nel
% sistema di riferimento cartesiano), calcola e resistuisce i parametri
% orbitali
%
% ++input++ 
% -s: vettore (orizzontale o verticale) contenente nelle prime tre posizioni le
%  componenti x,y,z della posizione[km] e nelle seguenti tre posizioni le
%  componenti x,y,z della velocità[km/s]
% -mu[km^3/s^2]: costante planetaria
%
% ++output++ 
% -a[km]: semiasse maggiore dell'orbita
% -e[-]: eccentricità dell'orbita
% -i[rad]: inclinazione dell'orbita nel piano tridimensionale
% -omega[rad]: ascensione retta del nodo ascendente
% -o[rad]: anomalia del pericentro 
% -theta[rad]: anomalia vera relativa alla posizione data in input


h = cross(s(1:3),s(4:6));
i = acos(h(3)/norm(h));
e1 = cross(s(4:6),h)/mu - s(1:3)/norm(s(1:3));
e = norm(e1);
a =  - mu/2 * ((norm(s(4:6)))^2/2 - mu/norm(s(1:3)))^( - 1);
z = [0,0,1];
x = [1,0,0];

if i == 0
    omega = 0;
    o = 0;
else
    N = cross(z,h)/norm(cross(z,h));
    if N(2) >= 0
        omega = acos(dot(x,N)/norm(N));
    else
        omega = 2 * pi - acos(dot(x,N)/norm(N));
    end
    if e1(3) >= 0
        o = acos(dot(e1,N)/e * norm(N));
    else
        o = 2 * pi - acos(dot(e1,N)/e * norm(N));
    end
end

v_r = dot(s(4:6),s(1:3)/norm(s(1:3)));
if v_r >= 0
    theta = acos(dot(e1,s(1:3))/(e * norm(s(1:3))));
else
    theta = 2 * pi - acos(dot(e1,s(1:3))/(e * norm(s(1:3))));
end

