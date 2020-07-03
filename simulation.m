state0 = [0;0;0];
%statef = [rand(1)*1000;rand(1)*1000;rand(1)*pi];
statef = [980;676;2.87];
t0 = 0;
tf = 1e7;
sim_dt = 0.001;
state = state0;
t = t0;
Tau_h = 90;
a = [10;0;0];
b = [0.35;0;0];
parameters = [a;b];
figure;
plot(state0(1),state0(2),'-kx','MarkerSize',15, 'LineWidth',2); hold on;

while t<tf
    [parameters,state,del_statef_p] = modelPredictiveControl(state,statef,parameters,t,Tau_h);
    t = t + sim_dt;
    plot(state(1),state(2),'-kx'); pause(0.01);
    if del_statef_p < 1
        statef
        state
        parameters
        break;
    end
end