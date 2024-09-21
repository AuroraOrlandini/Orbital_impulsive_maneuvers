function [Deltavtot,Deltat,thetaf,kepet,Deltav1,Deltav2] = changeOrbitShape(a1,e1,a2,e2,k)
% ++descrizione++
% Questa funzione, dati i parametri relativi ad orbita finale ed iniziale,
% restituisce i risultati relativi ad una generica manovra di cambio di
% forma,scelta adeguatamente seguento la conformità dei parametri. Grazie
% al parametro k questa funzione esegue i calcoli in modo da minimizzare
% la variazione di velocità o di tempo, in base alla scelta in input.
%
% ++input++
% -a1[km]: semiasse maggiore dell'orbita di partenza della manovra
% -e1[-]: eccentricità dell'orbita di partenza della manovra
% -a2[km]: semiasse maggiore dell'orbita su cui voglio trasferirmi tramite
%  questa manovra
% -e2[-]: eccentricità dell'orbita su cui voglio trasferirmi tramite
%  questa manovra
% -k[-]: parametro per la scelta della manovra, k = 0 convenienza in
%  velocità (minimizza il deltav), k = 1 convenienza in tempo (minimizza il
%  deltat)
%
% ++output++
% -Deltav1[km/s]: variazione di velocità necessaria per trasferirsi
%  dall'orbita 1 iniziale al pericentro o apocentro dell'orbita di
%  trasferimento
% -Deltav2[km/s]: variazione di velocità necessaria per trasferirsi
%  dall'orbita di trasferimento bitangente al pericentro o apocentro
%  dell'orbita 2 finale
% -Deltavtot[km/s]: variazione di velocità necessaria per effettuare questa
%  manovra, inclusiva di entrambi gli impulsi che comprende
% -thetaint[rad]: anomalia vera relativa alla posizione in cui si effettua 
%  la manovra di cambio piano
% -Deltat[s]: tempo trascorso sull'orbita di trasferimento (corrispondente
%  a metà del suo periodo) 
% -thetaf[rad]: anomalia vera della posizione (pericentro o apocentro)
%  raggiunta a seguito della manovra di cambio di forma scelta
% -kepet[a e]: vettore riga contenente in ordine semiasse maggiore
%  dell'orbita ed eccentricità dell'orbita di trasferimento 

mu = 398600;

p1 = a1 * (1 - e1^2);
p2 = a2 * (1 - e2^2);

r1_p = p1/(1 + e1);
r1_a = p1/(1 - e1);
r2_p = p2/(1 + e2);
r2_a = p2/(1 - e2);

v1_p = sqrt(mu/p1) * (1 + e1);
v1_a = sqrt(mu/p1) * (1 - e1);
v2_p = sqrt(mu/p2) * (1 + e2);
v2_a = sqrt(mu/p2) * (1 - e2);

[Deltava,kepeta,Dva] = trasf_bitangente(r1_p,v1_p,r2_a,v2_a);
[Deltavb,kepetb,Dvb] = trasf_bitangente(r1_a,v1_a,r2_p,v2_p);

if k == 0
    if Deltava<Deltavb
        fprintf('\n pericentro - apocentro più conveniente \n')
        Deltav1 = Dva(1);
        Deltav2 = Dva(2);
        Deltavtot = Deltava;
        thetaf = pi;
        kepet = kepeta;
    else
        fprintf('\n apocentro - pericento più conveniente \n')
        Deltav1 = Dvb(1);
        Deltav2 = Dvb(2);
        Deltavtot = Deltavb;
        thetaf = 0;
        kepet = kepetb;
    end
    Deltat = pi * sqrt((kepet(1)^3/mu));
end

if k == 1
    Deltata = pi * sqrt((kepeta(1)^3/mu));
    Deltatb = pi * sqrt((kepetb(1)^3/mu));
    if Deltata<Deltatb
        fprintf('\n pericentro - apocentro più conveniente \n')
        Deltav1 = Dva(1);
        Deltav2 = Dva(2);
        Deltavtot = Deltava;
        thetaf = pi;
        Deltat = Deltata;
        kepet = kepeta;
    else
        fprintf('\n apocentro - pericento più conveniente \n')
        Deltav1 = Dvb(1);
        Deltav2 = Dvb(2);
        Deltavtot = Deltavb;
        thetaf = 0;
        Deltat = Deltatb;
        kepet = kepetb;
    end
end





