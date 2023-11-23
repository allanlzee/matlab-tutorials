% Frequency of Building 
% Second Order Differential Equations 
% Euler's Method to Solve IVP 

% Measurements in mm
x1_0 = 100; 
x2_0 = 50; 

% Solve for Horizontal Displacements of the First 
% and Second Storeys (x1(t), x2(t))

% First derivatives at 0. 
x1_diff_0 = 0; 
x2_diff_0 = 0;

% Constants related to the second order differentials. 
k1 = 4.66; 
k2 = 4.66; 
m1 = 0.0917; 
m2 = 0.0765; 
% Equations 
% m1x1'' = -k1x1 - k2(x1 - x2) 
% m2x2'' = -k2(x2 - x1)

% A Matrix 
A = [0 1 0 0; (-k1-k2)/m1 0 k2/m1 0; 0 0 0 1; k2/m2 0 -k2/m2 0]; 

T = 10; 

[t, x1, x_diff, x2, y_diff] = IEMsolver(A, x1_0, x2_0, x1_diff_0, x2_diff_0, T, 1000); 

disp(t)
disp(x)
disp(y)

subplot(2, 1, 1)
plot(t, x1) 
title("x1")
xlabel("mm")
ylabel("kN/mm")

subplot(2, 1, 2)
plot(t, x2)
title("x2")
xlabel("mm")
ylabel("kN/mm")

function [t,x,x_diff,y,y_diff] = IEMsolver(A,x_0,y_0, x_diff_0, y_diff_0, T,N)
    dt = T/N;
    t = 0:dt:T;

    SOL = NaN(4,length(t));
    % initial values
    SOL(1,1) = x_0;
    SOL(3,1) = y_0;
    SOL(2,1) = x_diff_0; 
    SOL(4,1) = y_diff_0;

    for(k = 2:length(t))
        % Zn+1 = Zn,IEM + dt / 2 * (Azn + Azn+1,EM)
        % Last Zn       
        Zn = SOL(:,k-1); 
        Zn_1 = Zn + dt * (A * Zn); 
        SOL(:,k) = Zn + (dt/2) * (A * Zn + A * Zn_1); 
    end

    x = SOL(1,:);
    x_diff = SOL(2,:);
    y = SOL(3,:);
    y_diff = SOL(4,:);
end