function [Deltav,kepEt,Dv]=trasf_bitangente(R1,v1,R2,v2)
% ++Descrizione++
% Dati due coppie di vettore R,V restituisce i Dv (totali e parziali) 
% necessari per completare un trasferimento bitangente dal punto indicato 
% da {R1,v1} al punto indicato da {R2,v2}. La funzione riconosce quale dei 
% due raggi è maggiore e assegna in tal modo i valori a perigeo e apogeo per
% l'orbita di trasferimento.
%
% ++input++
%
% R1 [km}: vettore posizione del punto iniziale 
% v1 [km/s]: vettore velocità del punto iniziale
% R2 [km}: vettore posizione del punto finale 
% v2 [km/s]: vettore velocità del punto finale
%
% ++output++
%
% Deltav[km/s]: Deltav complessivo necessario per completare la manovra
% kepEt [km]: parametri orbitali dell'orbita di trasferimento {a,e} 
% Dv [km/s]: vettore (1 x 2) in cui sono presenti i valori dei due Dv

mu=398600;

aT = (R1+R2)/2;

if R1>R2
    eT = (R1-R2)/(R1+R2);
else
    eT = (R2-R1)/(R1+R2);
end

pT = aT*(1-eT^2);

vaT = sqrt(mu/pT)*(1-eT);
vpT = sqrt(mu/pT)*(1+eT);

if R1 > R2
    Deltav1 = vaT-v1;
    Deltav2 = v2-vpT;
else
    Deltav1 = vpT-v1;
    Deltav2 = v2-vaT;
end

Deltav = Deltav2+Deltav1;

kepEt = [aT,eT];
Dv = [Deltav1; Deltav2];