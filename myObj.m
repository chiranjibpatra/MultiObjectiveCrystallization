%%
%Hemalatha, K., & Rani, K. Y. (2017). Multiobjective optimization of unseeded and seeded batch cooling crystallization processes. Industrial & Engineering Chemistry Research, 56(20), 6012-6021.
%%
function obj = myObj(dT)

tf = 600;               
step = 14;             
tspan = linspace(0,tf,step+1)';
T = 273 + 48 + cumsum([0 dT]);
[~,x] = ode45(@(t,x) moment_eq(t,x,T,tspan),tspan,[0 0 0 0 0.56]);

cv = sqrt(x(end,3).*x(end,1)./x(end,2).^2-1);  
nms = (-x(end,2)/x(end,1));                       

obj = [nms*1e6;cv*1e2];                         
end