%% 内点法求解约束非线性最优化问题
% 将方程组转换为残差形式
residuals = @(x) biological_interactions(x)
    
% 定义目标函数（残差的平方和）
objective = @(x) sum(residuals(x).^2)
y = randi([0, 3], 1, 14);
for i=1:100
    x0 = randi([0, 10], 1, 14);   
    % 没有线性或非线性约束，因此 A, b, Aeq, beq, lb, ub 均为空
    A = -eye(14);
    b = zeros(14, 1);
%     A = [];
%     b = [];
    Aeq = [];
    beq = [];
    lb = []; % 如果有变量下界，则在这里定义
    ub = []; % 如果有变量上界，则在这里定义
    % 使用 fmincon 求解
        options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
        [x, fval] = fmincon(objective, x0, A, b, Aeq, beq, lb, ub, [], options);
        x0 = round(x * 10^8) / 10^8;
        % 显示结果
        disp('Solution:');
        disp(x0);
        disp('Objective function value at the solution:');
        disp(fval);

        % 检查残差是否足够小
        residual_norm = norm(residuals(x));
        disp(['Residual norm:', num2str(residual_norm)]);

        % 如果残差足够小，可以认为找到了方程组的解
        if residual_norm < 1e-6
            disp('The solution is considered acceptable.');
            y = [y;x0];
        else
            disp('The solution may not be accurate enough.');
        end
end
tabulate(y(:,1))
%% 网络
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


