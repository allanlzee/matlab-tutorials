[t, x, y] = EMsolver([0 -1;1 0], 1, 0, 10, 1000); 
x_exact = cos(t);
y_exact = sin(t);

disp(size(t))
disp(size(x))

[tIEM, xIEM, yIEM] = IEMsolver([0 -1;1 0], 1, 0, 10, 1000); 

subplot(2, 1, 1)
hold on
plot(t, x_exact, "b--")
plot(t, x)
plot(t, xIEM)
title("x(t)")
xlabel("Time")
ylabel("Value")
hold off

subplot(2, 1, 2)
hold on
plot(t, y_exact, "r--")
plot(t, y)
plot(t, yIEM)
title("y(t)")
xlabel("Time")
ylabel("Value")
hold off


function [t,x,y] = EMsolver(A,x_0,y_0,T,N)
    dt = T/N;
    t = 0:dt:T;

    SOL = NaN(2,length(t));
    % Initial values
    SOL(1,1) = x_0;
    SOL(2,1) = y_0;

    for(k = 2:length(t))
        SOL(:,k) = SOL(:,k-1) + dt*A*SOL(:,k-1);
    end

    x = SOL(1,:);
    y = SOL(2,:);
end


function [t,x,y] = IEMsolver(A,x_0,y_0,T,N)
    dt = T/N;
    t = 0:dt:T;

    SOL = NaN(2,length(t));
    % initial values
    SOL(1,1) = x_0;
    SOL(2,1) = y_0;

    for(k = 2:length(t))
        % Zn+1 = Zn,IEM + dt / 2 * (Azn + Azn+1,EM)
        % Last Zn       
        Zn = SOL(:,k-1); 
        Zn_1 = Zn + dt * (A * Zn); 
        SOL(:,k) = Zn + (dt/2) * (A * Zn + A * Zn_1); 
    end

    x = SOL(1,:);
    y = SOL(2,:);
end
