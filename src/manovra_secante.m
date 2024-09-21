function [Deltav, thetamanovra_1,thetamanovra_2] = manovra_secante(theta1,theta2,p1,e1,thetaacuisono)
% ++descrizione++
% Questa funzione, dati i parametri necessari relativi ad orbita finale ed
% iniziale, restituisce i risultati relativi ad una generica manovra di 
% cambio di anomalia del pericentro.
%
% ++input++
% -theta1[rad]: anomalia del pericentro dell'orbita di partenza della manovra
% -theta2[rad]: anomalia del pericentro dell'orbita finale desiderata
% -p1[km]: semilato retto dell'orbita di partenza
% -e1[-]: eccentricità dell'orbita di partenza
% -thetaincuisono[rad]: anomalia vera della posizione corrente, precedente
%  alla manovra
%
% ++output++
% -Deltav[km/s]: variazione di velocità necessaria per compiere la manovra
% -thetamanovra_1[rad]: anomalia vera del punto di manovra calcolata
%  rispetto ai parametri dell'orbita di partenza
% -thetamanovra_2[rad]: anomalia vera del punto di manovra calcolata
%  rispetto ai parametri dell'orbita di arrivo


mu = 398600;

delta_o = theta2 - theta1;

if delta_o < 0
    if thetaacuisono<pi - abs(delta_o/2)
        thetamanovra_2 = pi + abs(delta_o/2);
        thetamanovra_1 = pi - abs(delta_o/2);
    else
        thetamanovra_1 = 2 * pi - abs(delta_o/2);
        thetamanovra_2 = abs(delta_o/2);
    end
else
    if thetaacuisono<pi - abs(delta_o/2)
        thetamanovra_1 = pi + abs(delta_o/2);
        thetamanovra_2 = pi - abs(delta_o/2);
    else
        thetamanovra_2 = 2 * pi - abs(delta_o/2);
        thetamanovra_1 = abs(delta_o/2);
    end
end


Deltav = 2 * sqrt(mu/p1) * e1 * sin(thetamanovra_2);

