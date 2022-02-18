%%
%Hemalatha, K., & Rani, K. Y. (2017). Multiobjective optimization of unseeded and seeded batch cooling crystallization processes. Industrial & Engineering Chemistry Research, 56(20), 6012-6021.
%%
%% Setup
tf = 600;               
step = 14;                      
Pop_Size = 200;          
Aeq = ones(1,step);
Beq = 18-48;

dT_init = repmat(-(48-18)/step,Pop_Size,step);
figure
plot(linspace(0,tf,step+1),48 + [zeros(size(dT_init,1),1) cumsum(dT_init,2)])
title("Input Temperature Trajectory")
xlabel('Time (min)'); ylabel(['Temperature (',char(176),'C)'])
%% Optimisation part
options = optimoptions('gamultiobj','InitialPopulationMatrix',dT_init,'PopulationSize',Pop_Size,'PlotFcn','gaplotpareto','FunctionTolerance',1e5);
[dT_Opt,fval] = gamultiobj(@myObj,step,[],[],Aeq,Beq,repmat(-1.35*tf/step,step,1),zeros(step,1),[],options);
%% Describing results
results = [fval,dT_Opt];
sorted_results = sortrows(results);
sortedfval = sorted_results(:,1:2);
figure(1)
plot(-sortedfval(:,1),sortedfval(:,2),'.',"LineWidth",2), xlabel('NMS ($\mu m$)',"Interpreter","latex"); ylabel('% CV')

tspan = linspace(0,tf,step+1);
figure
plot(tspan,48+cumsum([0 dT_Opt(find(fval(:,1) == min(fval(:,1)),1),:)]),'--o'),hold on
plot(tspan,48+cumsum([0 dT_Opt(find(fval(:,2) == min(fval(:,2)),1),:)]),'-.o'),
title("Temperature Trajectory for min %CV and Max NMS")
xlabel('Time (min)'); ylabel(['Temperature (',char(176),'C)'])
legend("For Min. CV","For Max. NMS")

p = 32;
figure(1)
hold on, plot(-sortedfval(p,1), sortedfval(p,2), 'O',"LineWidth",2)
T_Opt = 273+48+cumsum([0 dT_Opt(p,:)]);
figure
plot(tspan,T_Opt - 273,'-.x');title("Temperature Trajectory for optimal solution")
xlabel('Time (min)'); ylabel(['Temperature (',char(176),'C)'])

[t,x] = ode45(@(t,x) moment_eq(t,x,T_Opt,tspan),[0,600],[ 0 0 0 0 0.56]);
Cc = x(:,5);
figure
subplot 321, plot(t,x(:,1)); xlabel('Time (min)'); ylabel('Zeroth Moment')
subplot 322, plot(t,x(:,2)); xlabel('Time (min)'); ylabel('First Moment')
subplot 323, plot(t,x(:,3)); xlabel('Time (min)'); ylabel('Second Moment')
subplot 324, plot(t,x(:,4)); xlabel('Time (min)'); ylabel('third Moment')
subplot 325, plot(t,x(:,5)); xlabel('Time (min)'); ylabel('Concentration (g/ml)')

