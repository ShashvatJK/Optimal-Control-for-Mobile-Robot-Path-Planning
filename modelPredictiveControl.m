function [parameters,state,del_statef_p] = modelPredictiveControl(state,statef,parameters,t,Tau_h)

    L = 100;
    dt = 0.01;
    dist = norm(statef(1:2)-state(1:2));
    vel = dist/(Tau_h*dt);
    x = zeros(1,Tau_h+1);
    y = zeros(1,Tau_h+1);
    theta = zeros(1,Tau_h+1);
    setGlobalState([x;y;theta]);
    parameters = fmincon(@(parameters) trajectoryGen(parameters,state,statef,Tau_h,L,t,dt),parameters,[],[],[],[],...
        [-inf;-inf;-inf;-inf;-inf;-inf],[inf;inf;inf;inf;inf;inf],...
        @(parameters) nonLinearConstraint(parameters,Tau_h,dt,vel));
    del_statef_p=trajectoryGen(parameters,state,statef,Tau_h,L,t,dt)
    states = getGlobalState;
    state = states(:,Tau_h+1);
    
    hold on;
    xlabel("X coordinate");
    ylabel("Y coordinate");
    hold on;
    plot(statef(1),statef(2),'-xr','MarkerSize',15,'LineWidth',2) 
    hold on;
    curve2 = animatedline('LineWidth',1.5);
    grid on;
    for i=1:length(states)
        %plot(states(1,:),states(2,:),'Color','black', 'LineWidth',2); grid on;
        addpoints(curve2,states(1,i),states(2,i));
        drawnow
    end
    hold on;
end

function del_statef_p = trajectoryGen(parameters,state,statef,Tau_h,L,t,dt)
    states = getGlobalState;
    states(:,1)=state;
    for j = 2:Tau_h+1 
        states(1,j) = states(1,j-1) + polyval(parameters(1:3,1),t+(j-1)*dt).*cos(states(3,j-1))*dt;
        states(2,j) = states(2,j-1) + polyval(parameters(1:3,1),t+(j-1)*dt).*sin(states(3,j-1))*dt;
        states(3,j) = states(3,j-1) + polyval(parameters(1:3,1),t+(j-1)*dt).*tan(pi_to_pi(polyval(parameters(4:6,1),t+(j-1)*dt)))'*dt/L;
        states(3,j) = pi_to_pi(states(3,j));
    end
    setGlobalState(states);
    del_statef_p = immse([states(1:2,Tau_h+1);states(3,Tau_h+1)*1] , [statef(1:2);statef(3)*1]);
end

function [c,ceq] = nonLinearConstraint(parameters, Tau_h, dt, vel)
    c = zeros(2*Tau_h,1);
    for j = 1:Tau_h 
        c(j) =  abs(polyval(parameters(1:3,1),(j)*dt)) - 2*vel;
    end
    for j = Tau_h+1:2*Tau_h
        c(j) =  abs(polyval(parameters(4:6,1),(j-Tau_h)*dt)) - pi/6;
    end
    c(Tau_h) = abs(polyval(parameters(1:3,1),Tau_h*dt));
    c(2*Tau_h) = abs(polyval(parameters(4:6,1),Tau_h*dt));
    ceq= [];
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