%% The interior point method is used to solve constrained nonlinear optimization problems
% The system of equations is converted to residual form
residuals = @(x) biological_interactions(x)
    
% Define the objective function (sum of squares of residuals)
objective = @(x) sum(residuals(x).^2)
y = randi([0, 3], 1, 14);
for i=1:100
    x0 = randi([0, 10], 1, 14);   
    % There are no linear or nonlinear constraints, so A, b, Aeq, beq, lb, ub are all empty
    A = -eye(14);
    b = zeros(14, 1);
%     A = [];
%     b = [];
    Aeq = [];
    beq = [];
    lb = []; % If there is a lower bound for a variable, it is defined here
    ub = []; % If there is an upper bound on a variable, it is defined here
    % fmincon solve
        options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
        [x, fval] = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, [], options);
        x0 = round(x * 10^8) / 10^8;
        % displaying results
        disp('Solution:');
        disp(x0);
        disp('Objective function value at the solution:');
        disp(fval);

        % Check if the residual is small enough
        residual_norm = norm(residuals(x));
        disp(['Residual norm:', num2str(residual_norm)]);

        % If the residual is small enough, the solution to the system can be considered found
        if residual_norm < 1e-6
            disp('The solution is considered acceptable.');
            y = [y;x0];
        else
            disp('The solution may not be accurate enough.');
        end
end
tabulate(y(:,1))
%% GRN without BRD
function F = biological_interactions(x) %the biological interactions between genes
a=1.1;
b=2;
k=1;
n=3;
s=0.5;
F(1)=a*x(4)^n/(s^n+x(4)^n)+a*x(5)^n/(s^n+x(5)^n)+a*x(8)^n/(s^n+x(8)^n)-k*x(1);
F(2)=a*x(9)^n/(s^n+x(9)^n)-k*x(2);
F(3)=a*x(5)^n/(s^n+x(5)^n)-k*x(3);
F(4)=a*x(1)^n/(s^n+x(1)^n)-k*x(4);
F(5)=a*x(1)^n/(s^n+x(1)^n)+a*x(2)^n/(s^n+x(2)^n)+a*x(3)^n/(s^n+x(3)^n)+a*x(13)^n/(s^n+x(13)^n)+b*s^n/(s^n+x(9)^n)+b*s^n/(s^n+x(10)^n)+b*s^n/(s^n+x(6)^n)+b*s^n/(s^n+x(14)^n)-k*x(5);
F(6)=a*x(1)^n/(s^n+x(1)^n)+a*x(12)^n/(s^n+x(12)^n)+b*s^n/(s^n+x(5)^n)+b*s^n/(s^n+x(10)^n)-k*x(6);
F(7)=a*x(5)^n/(s^n+x(5)^n)-k*x(7);
F(8)=a*x(2)^n/(s^n+x(2)^n)+a*x(14)^n/(s^n+x(14)^n)+b*s^n/(s^n+x(9)^n)-k*x(8);
F(9)=a*x(2)^n/(s^n+x(2)^n)+a*x(3)^n/(s^n+x(3)^n)+a*x(4)^n/(s^n+x(4)^n)+a*x(6)^n/(s^n+x(6)^n)+a*x(1)^n/(s^n+x(1)^n)-k*x(9);
F(10)=a*x(1)^n/(s^n+x(1)^n)+a*x(2)^n/(s^n+x(2)^n)+a*x(6)^n/(s^n+x(6)^n)+b*s^n/(s^n+x(10)^n)-k*x(10);
F(11)=b*s^n/(s^n+x(11)^n)+a*x(2)^n/(s^n+x(2)^n)+a*x(4)^n/(s^n+x(4)^n)+a*x(6)^n/(s^n+x(6)^n)-k*x(11);
F(12)=b*s^n/(s^n+x(5)^n)+b*s^n/(s^n+x(6)^n)+b*s^n/(s^n+x(1)^n)-k*x(12);
F(13)=a*x(5)^n/(s^n+x(5)^n)+a*x(12)^n/(s^n+x(12)^n)-k*x(13);
F(14)=a*x(14)^n/(s^n+x(14)^n)-k*x(14);
end


