epsilon = 1e-7;
state0 = [0;0;0];
statef = [rand(1)*1000;rand(1)*1000;rand(1)*pi];
%statef = [234;189;1.75];
dist = norm(statef(1:2)-state0(1:2));
a = [10;0;0];
b = [0.35;0;0];
L = 100;
dt = 0.01;
N = 50;
NN = 100;
vel = dist/(N*dt);
[parameters,states] = modelPredictiveControl(N,NN,dt,L,a,b,state0,statef,vel);

function [parameters,states] = modelPredictiveControl(N,NN,dt,L,a,b,state0,statef,vel)
    x = zeros(1,NN+1);
    y = zeros(1,NN+1);
    theta = zeros(1,NN+1);
    parameters = [a;b];
    setGlobalState([x;y;theta]);
    states = getGlobalState;
    %options = optimoptions('fmincon','Display','iter');
    for k = 1:NN-N+1
        [parameters,fval] = fmincon(@(parameters) trajectoryGen(parameters,state0,statef,N,L,dt,k),parameters,[],[],[],[],...
            [-inf;-inf;-inf;-inf;-inf;-inf],[inf;inf;inf;inf;inf;inf],...
            @(parameters) nonLinearConstraint(parameters,states,statef,N,NN,dt,vel,k));
        del_statef_p=trajectoryGen(parameters,state0,statef,N,L,dt,k);
    end
    fval
    states = getGlobalState;
    step = 2;
    hold on;
    plot(states(1,1),states(2,1),'-xk','MarkerSize',15,'LineWidth',2);
    xlabel("X coordinate");
    ylabel("Y coordinate");
    hold on;
    plot(statef(1),statef(2),'-xr','MarkerSize',15,'LineWidth',2) 
    hold on;
    %Sf = polyval(parameters(1:3,1),(N)*dt);
    %quiver(statef(1),statef(2),Sf.*cos(statef(3)),Sf.*sin(statef(3)),...
    % 'Color','black', 'LineWidth',1.5,'MaxHeadSize', 1.2,'AutoScaleFactor',1);
    %hold on;
    curve2 = animatedline('LineWidth',3);
    grid on;
    for i=1:length(states)
        addpoints(curve2,states(1,i),states(2,i));
        drawnow
    end
    hold on;
    %S = polyval(parameters(1:3,1),(1:step:N)*dt);
    %quiver(states(1,2:step:end),states(2,2:step:end),S.*cos(states(3,2:step:end)),S.*sin(states(3,2:step:end)),...
    % 'Color','red', 'LineWidth',0.5,'MaxHeadSize', 0.2,'AutoScaleFactor',1);
    %hold on;
    %Sf = polyval(parameters(1:3,1),(N)*dt);
    %quiver(states(1,end),states(2,end),Sf.*cos(states(3,end)),Sf.*sin(states(3,end)),...
    % 'Color','blue', 'LineWidth',1.5,'MaxHeadSize', 1.2,'AutoScaleFactor',1);
    %hold on;
    plot(states(1,end),states(2,end),'-xb','MarkerSize',15,'LineWidth',2);
    statef
    states(:,end)
end

function del_statef_p = trajectoryGen(parameters, state0, statef, N, L, dt, k)
    states = getGlobalState;
    states(:,1)=state0;
    for j = 1+k:N+k 
    %     inputs(1,j-1) = polyval(parameters(1:3,1),(j-1)*dt);
    %     inputs(2,j-1) = polyval(parameters(4:6,1),(j-1)*dt);
        states(1,j) = states(1,j-1) + polyval(parameters(1:3,1),(j-1)*dt).*cos(states(3,j-1))*dt;
        states(2,j) = states(2,j-1) + polyval(parameters(1:3,1),(j-1)*dt).*sin(states(3,j-1))*dt;
        states(3,j) = states(3,j-1) + polyval(parameters(1:3,1),(j-1)*dt).*tan(pi_to_pi(polyval(parameters(4:6,1),(j-1)*dt)))'*dt/L;
        states(3,j) = pi_to_pi(states(3,j));
    end
    setGlobalState(states);
    del_statef_p = immse([states(1:2,N+k);states(3,N+k)*2] , [statef(1:2);statef(3)*2]);
end

function [c,ceq] = nonLinearConstraint(parameters, states,statef, N, NN, dt, vel, k)
    c = zeros(2*(NN),1);
    for j = k:N+k-1
        c(j) =  abs(polyval(parameters(1:3,1),(j)*dt)) - vel;
    end
    for j = (length(c)/2)+k:(length(c)/2)+N-1+k
        c(j) =  abs(polyval(parameters(4:6,1),(j-N)*dt)) - pi/4;
    end
    c(N+k-1) = abs(polyval(parameters(1:3,1),(N+k-1)*dt));
    c(end/2 + N+k-1) = abs(polyval(parameters(4:6,1),(N+k-1)*dt));
    ceq=[];
end

function angle = pi_to_pi(angle)
    angle = mod((angle+pi),(2*pi))-pi;
end

function setGlobalState(val)
    global state
    state = val;
end

function r = getGlobalState
    global state
    r = state;
end