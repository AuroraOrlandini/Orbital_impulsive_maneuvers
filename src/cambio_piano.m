function [Deltav,thetaint,u1,omegatrasf,dt,u2] = cambio_piano(i1,i2,Omega1,Omega2,o1,p1,e1,a1,thetadovesono,k)
% ++descrizione++
% Questa funzione, dati i parametri relativi ad orbita finale ed iniziale,
% restituisce i risultati relativi ad una generica manovra di cambio di
% piano, che modifica quindi inclinazione i, ascensione retta del nodo
% ascendente Omega e anomalia del pericentro o. Grazie al parametro k
% questa funzione esegue i calcoli in modo da minimizzare la variazione di
% velocità o di tempo, in base alla richiesta in input.
%
% ++input++
% -i1[rad]: inclinazione dell'orbita precedente alla manovra
% -i2[rad]: inclinazione dell'orbita dopo la manovra
% -Omega1[rad]: ascensione retta del nodo ascendente dell'orbita precedente
%  alla manovra
% -Omega2[rad]: ascensione retta del nodo ascendente dell'orbita dopo
%  la manovra
% -o1[rad]: anomalia del pericentro dell'orbita precedente alla manovra 
% -p1[km]: semilato retto dell'orbita precedente alla manovra
% -e1[-]: eccentricità dell'orbita precedente alla manovra
% -a1[km]: semiasse maggiore dell'orbita precedente alla manovra
% -thetadovesono[rad]: anomalia vera della posizione corrente sull'orbita,
% prima di iniziare la manovra
% -k[-]: parametro per la scelta della manovra, k = 0 convenienza in velocità, k = 1 convenienza in tempo
% 
% ++output++
% -Deltav[km/s]: variazione di velocità necessaria per effettuare questa
%  manovra
% -thetaint[rad]: anomalia vera relativa alla posizione in cui si effettua la manovra di cambio piano
% -u1[rad]: lato del triangolo sferico relativo all'orbita precedente alla
%  manovra
% -omegatrasf[rad]: anomalia al pericentro che caratterizza la nuova orbita
%  a seguito della manovra
% -dt[s]: tempo necessario per spostarmi dal thetadovesono al thetaint
% -u2[rad]: lato del triangolo sferico relativo all'orbita ottenuta a
%  seguito della manovra


mu = 398600;

deltaOmega = Omega2 - Omega1;
deltai = i2 - i1;

alpha = acos(cos(i1) * cos(i2)+sin(i1) * sin(i2) * cos(deltaOmega));

u1 = asin(sin(deltaOmega)/sin(alpha) * sin(i2));
u2 = asin(sin(deltaOmega)/sin(alpha) * sin(i1));


if deltaOmega > 0
    if deltai > 0 || deltai == 0
        thetaint = u1 - o1;
        omegatrasf = u2 - thetaint;
    else
        thetaint = 2 * pi - u1 - o1;
        omegatrasf = 2 * pi - u2 - thetaint;
    end
else
    if deltai > 0 || deltai == 0
        thetaint = 2 * pi - u1 - o1;
        omegatrasf = 2 * pi - u2 - thetaint;
    else
        thetaint = u1 - o1;
        omegatrasf = u2 - thetaint;
    end
end

v1 = sqrt(mu/p1) * (1+e1 * cos(thetaint));
v2 = sqrt(mu/p1) * (1+e1 * cos(thetaint + pi));
dt1 = calcolo_tempi(thetadovesono,thetaint,a1,e1);
dt2 = calcolo_tempi(thetadovesono,thetaint + pi,a1,e1);

if k == 0
    if v2 < v1
        thetaint = thetaint+pi;
        v = v2;
        dt = dt2;
    else
        v = v1;
        dt = dt1;
    end
else
    if dt2 < dt1
        thetaint = thetaint + pi;
        v = v2;
        dt = dt2;
    else
        v = v1;
        dt = dt1;
    end
end


Deltav = 2 * v * sin(alpha/2);


