load lawdata.mat

x = gpa;
y = lsat;

x_plot = linspace(2.0, 4);

% Linear Best Line of Fit
% y = a + bx 
A_linear = ones(length(gpa), 2); 

% Put y in the first column.
for i = 1:length(gpa)
    A_linear(i, 2) = gpa(i);
end

% Solve for A*. This will be AT * A. This will result in a 
% square matrice. 
A_star = transpose(A_linear) * A_linear;
B_star = transpose(A_linear) * y; 

x_linear = A_star\B_star;
y_linear = x_linear(1) + x_linear(2) * x_plot;

error_linear = (norm((y - A_linear * (A_star\B_star))).^2) / 15;
disp(error_linear)

% ------------------------------------------------%

% Quadratic Best Line of Fit 
% y = a + bx + cx^2 
A_quadratic = ones(length(gpa), 3); 

% Put y values in the first column and y^2 in the second column.
for col = 2:3
    for i = 1:length(gpa) 
        A_quadratic(i, col) = gpa(i).^(col - 1);
    end
end

A_star_quadratic = transpose(A_quadratic) * A_quadratic; 
B_star_quadratic = transpose(A_quadratic) * y; 

x_quadratic = A_star_quadratic\B_star_quadratic; 
y_quadratic = x_quadratic(1) + x_quadratic(2) * x_plot + x_quadratic(3) * (x_plot.^2);

error_quadratic = (norm((y - A_quadratic* (A_star_quadratic\B_star_quadratic))).^2) / 15; 
disp(norm(error_quadratic))
% ------------------------------------------------%

% Cubic Best Line of Fit 
% y = a + bx + cx^2 + dx^3
A_cubic = ones(length(gpa), 4);

% Put y in first column, y^2 in second column, and y^3 in third column. 
for col = 2:4 
    for i = 1:length(gpa)
        A_cubic(i, col) = gpa(i).^(col - 1); 
    end
end 

A_star_cubic = transpose(A_cubic) * A_cubic; 
B_star_cubic = transpose(A_cubic) * y; 

x_cubic = A_star_cubic\B_star_cubic; 
y_cubic = x_cubic(1) + x_cubic(2) * x_plot + x_cubic(3) * (x_plot.^2) + x_cubic(4) * (x_plot.^3); 

% Error Vectors 
% Calculate associated error vector 
% e = b - AxLS 
% Known: xLS = 
error_cubic = (norm((y - A_cubic * (A_star_cubic\B_star_cubic))).^2) / 15;
disp(norm(error_cubic))

% Least Squares Problem
plot(x, y, 'o')

hold on 

plot(x_plot, y_linear, Color="red")

plot(x_plot, y_quadratic, Color="blue")

plot(x_plot, y_cubic, '--', Color="green")

hold off