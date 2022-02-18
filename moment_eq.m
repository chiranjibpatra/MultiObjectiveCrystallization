%%
%Hemalatha, K., & Rani, K. Y. (2017). Multiobjective optimization of unseeded and seeded batch cooling crystallization processes. Industrial & Engineering Chemistry Research, 56(20), 6012-6021.
%%
function y = moment_eq(t,x,T,tspan)

u0 = x(1);
u1 = x(2);
u2 = x(3);
u3 = x(4);
C = x(5);

T = interp1(tspan,T,t);
Csat = -2.0282 + 0.36592*(T-273) - 0.025618*(T-273)^2 + 9.7964e-4*(T-273)^3 - 2.1062e-5*(T-273)^4 + 2.4309e-7*(T-273)^5 - 1.18e-9*(T-273)^6;

G = 1.926e-2*exp(-25800/(8.314*T))*Csat*(C/Csat-1);
B = 6.9e22*exp(-76700/(8.314*T))*exp(-0.16/(log(C/Csat))^2);
rhoc = 1.4;
kv = 1;

y =[B
    G*u0
    G*u1*2
    G*u2*3
    -3*rhoc*G*u2*kv];
end