% Arithmetic Series 
% Sn = n(a1 + an / 2)

final_area = (1 / (1 - 1/4)) + zeros(20);

n = 1:20;

% Vector for each point of arithmetic series 
arithmetic = 1:20; 
area = zeros(20);

for i = 0:19 
    arithmetic(i + 1) = (2 ^ i) * (1 / 8) ^ (i); 
end 

disp(arithmetic)

for i = 1:20
    sum = 0; 
    for arith = 1:i
        sum = sum + arithmetic(arith);
    end 

    area(i) = sum;
end

ylabel("Area")
xlabel("Term Number")
title("Arithmetic Sequence Area Approximation")
plot(n, final_area, "-")
hold on 

plot(n, area)
legend("Final Area", "Area")
hold off
