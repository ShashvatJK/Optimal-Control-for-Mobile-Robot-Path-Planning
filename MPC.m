epsilon = 1e-7;
state0 = [0;0;0];
statef = [300;400;0.55];
dist = norm(statef(1:2)-state0(1:2));
a = [10;0;0];
b = [0.35;0;0];
L = 10;
dt = 0.001;
N = 100;
vel = dist/(N*dt);
[parameters,states] = modelPredictiveControl(N,dt,L,a,b,state0,statef,vel);

function [parameters,states] = modelPredictiveControl(N,dt,L,a,b,state0,statef,vel)
    x = zeros(1,N+1);
    y = zeros(1,N+1);
    theta = zeros(1,N+1);
    parameters = [a;b];
    setGlobalState([x;y;theta]);
    %c = nonLinearConstraint
    states = getGlobalState;
    parameters = fmincon(@(parameters) trajectoryGen(parameters,state0,statef,N,L,dt),parameters,[],[],[],[],...
        [-inf;-inf;-inf;-inf;-inf;-inf],[inf;inf;inf;inf;inf;inf],...
        @(parameters) nonLinearConstraint(parameters,states,statef,N,dt,vel));
    del_statef_p=trajectoryGen(parameters,state0,statef,N,L,dt);
    states = getGlobalState;
    step = 2;
    hold on;
    plot(states(1,1),states(2,1),'-xk','MarkerSize',15,'LineWidth',2);
    xlabel("X coordinate");
    ylabel("Y coordinate");
    hold on;
    plot(statef(1),statef(2),'-xk','MarkerSize',15,'LineWidth',2) 
    hold on;
    hold on;
    Sf = polyval(parameters(1:3,1),(N)*dt);
    quiver(statef(1),statef(2),Sf.*cos(statef(3)),Sf.*sin(statef(3)),...
     'Color','black', 'LineWidth',1.5,'MaxHeadSize', 1.2,'AutoScaleFactor',1);
    hold on;
    curve2 = animatedline('LineWidth',3);
    %set(gca,'Xlim',[-500 500],'Ylim',[-500 500]);
    grid on;
    for i=1:length(states)
        %plot(states(1,:),states(2,:),'Color','black', 'LineWidth',2); grid on;
        addpoints(curve2,states(1,i),states(2,i));
        drawnow
    end
    hold on;
    S = polyval(parameters(1:3,1),(1:step:N)*dt);
    quiver(states(1,2:step:end),states(2,2:step:end),S.*cos(states(3,2:step:end)),S.*sin(states(3,2:step:end)),...
     'Color','red', 'LineWidth',0.5,'MaxHeadSize', 0.2,'AutoScaleFactor',1);
    hold on;
    Sf = polyval(parameters(1:3,1),(N)*dt);
    quiver(states(1,end),states(2,end),Sf.*cos(states(3,end)),Sf.*sin(states(3,end)),...
     'Color','blue', 'LineWidth',1.5,'MaxHeadSize', 1.2,'AutoScaleFactor',1);
    hold on;
    plot(states(1,end),states(2,end),'-xb','MarkerSize',15,'LineWidth',2);
    statef
    states(:,end)
    (del_statef_p)
    %Us = [1 dt*k (dt^2)*(k^2)]*parameters(1:3,:);
    %Uphi = [1 dt*k (dt^2)*(k^2)]*parameters(4:6,:);

end

function del_statef_p = trajectoryGen(parameters, state0, statef, N, L, dt)
    states = getGlobalState;
    states(:,1)=state0;
    for j = 2:N+1 
    %     inputs(1,j-1) = polyval(parameters(1:3,1),(j-1)*dt);
    %     inputs(2,j-1) = polyval(parameters(4:6,1),(j-1)*dt);
        states(1,j) = states(1,j-1) + polyval(parameters(1:3,1),(j-1)*dt).*cos(states(3,j-1))*dt;
        states(2,j) = states(2,j-1) + polyval(parameters(1:3,1),(j-1)*dt).*sin(states(3,j-1))*dt;
        states(3,j) = states(3,j-1) + polyval(parameters(1:3,1),(j-1)*dt).*tan(pi_to_pi(polyval(parameters(4:6,1),(j-1)*dt)))'*dt/L;
        states(3,j) = pi_to_pi(states(3,j));
    end
    setGlobalState(states);
    del_statef_p = immse([states(1:2,N+1);states(3,N+1)], [statef(1:2);statef(3)]);
end

function [c,ceq] = nonLinearConstraint(parameters, states,statef, N, dt, vel)
    c = zeros(2*N+3,1);
    %ceq = zeros(2*N+3,1);
    for j = 1:N 
        c(j) =  abs(polyval(parameters(1:3,1),(j)*dt)) - 3*vel;
    end
    for j = N+1:2*N
        c(j) =  abs(polyval(parameters(4:6,1),(j-N)*dt)) - pi/6;
    end
    c(N) = abs(polyval(parameters(1:3,1),N*dt));
    c(2*N) = abs(polyval(parameters(4:6,1),N*dt));
    ceq(1)= immse(states(1,end),statef(1));
    ceq(2)= immse(states(2,end),statef(2));
    ceq(3)= immse(states(3,end),statef(3));
    %ceq(4:5) = [abs(polyval(parameters(1:3,1),N*dt));abs(polyval(parameters(4:6,1),N*dt))];
end
% inputs(1,1) = a(1,1);
% inputs(2,1) = b(1,1);
% syms t

% for i = 2:N
%     
%     inputs(1,i) = subs([1 t t^2]*parameters(1:3,i),(i-1)*dt); %speed
%     inputs(2,i) = subs([1 t t^2]*parameters(4:6,i),(i-1)*dt); %phi
% end

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
