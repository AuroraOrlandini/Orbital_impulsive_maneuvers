function [theta,E,n,M]=T_theta(deltaT,a,e,mu)

% ricavo theta(t) passando per E anomalia eccentrica
% theta restituito in radianti

p=a*(1-e^2);
h=sqrt(mu*p);
n=sqrt(mu/a^3); % velocità angolare media
M=n*deltaT; % anomalia media

if e==0
    theta=deltaT*mu^2/h^3;
    return
end

if e==1
    phi=@(tet) mu^2/h^2*deltaT-1/2*tan(tet/2)*(1+1/3*(tan(tet/2))^2);
    theta=ptofis(0,phi);
    return
end

% f=@(E) 1/(1-e^2)^3*(E-e*sin(E))-mu^2/h^3*deltaT;
% df=@(E) 1/(1-e^2)^3*(1-e*cos(E));

f=@(E) E-e*sin(E)-M;
df=@(E) 1-e*cos(E); 
E=newton(0,f,df);

theta=2*atan(((1+e)/(1-e)^(-1/2))*tan(E/2));
