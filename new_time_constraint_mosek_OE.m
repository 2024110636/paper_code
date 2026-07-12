%% 使用 MOSEK 求解器（CVX接口）
% 单位说明：功率 W，时间 s，频率 MHz，速率 Mbit/s

clear all;  
cvx_clear;  
cvx_solver mosek;
clear all;
cvx_clear;
cvx_solver mosek;
cvx_precision low
cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', 300);
cvx_solver_settings('MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-3);
cvx_expert true;

%% 参数设置
M = 10         % 光纤数量
N = 7;         % 设备数量
L0 = 5;         % 基础距离 (km)
r = 1.5;        % 每个设备增加距离 (km)
alpha_loss = 0.2;   % 光纤损耗 (dB/km)
P_total = 25;       % 总供电功率 (W)
T_max = 1.4;        % 最大允许时延 (s)    
oe_eta = 0.3; % 新增光电转化效率

% 设备功率模型
P_static = 0.1;    % 静态功率 (W)
C_S = 0.2;            % 采样功率系数 (W/MHz)
C_CPU = 0.001;       % CPU功率系数 (W/GHz^3)
s = 2;              % 每个样本大小 (bit)
T0 = 1;             % 采样时间 (s)
n_flop = 100;       % 每个样本浮点操作数
R = 500;            % 光纤速率 (Mbit/s)

% CPU & 采样频率 (单位：Hz)
f_CPU_min = 500; f_CPU_max = 1000;   % MHz
f_s_min = 1;      f_s_max = 2;        % MHz

% 光纤功率约束
P_fiber_min = 0; P_fiber_max = 3;     % W

% 光纤损耗系数 alpha_f(j)
L = @(j) L0 + r * (2 * j - 1);
alpha_f = arrayfun(@(j) 10^(-alpha_loss * L(j)/10), 1:N);

% 发射功率 P_tx(j)
P_tx = arrayfun(@(j) 0.2 * 10^((2 + 0.2 * L(j)) / 10) / 1000, 1:N);

% Dinkelbach 参数
   
theta = 0;
epsilon = 1e-4;    
max_iter = 50;    
iter = 0;          

%% Dinkelbach 迭代
while iter < max_iter  
    iter = iter + 1;  
    
    cvx_begin
        
        cvx_solver mosek;
        % 决策变量
        variable x(M, N) binary
        variable P_fiber(M) nonnegative
        variable f_CPU(N) nonnegative
        variable f_s(N) nonnegative

        % 辅助变量
        variable t(N) nonnegative          % 1 / f_CPU
        variable f(N) nonnegative          % f_s * t
        variable w(M, N) nonnegative       % P_fiber 分配 x * P
        variable z(M, N) nonnegative       % w * f
        variable u(M, N) nonnegative       % w * f_s

        % A_x 与 B_x（目标函数）
        A_x = s * T0 * sum(f_s);  
        
        B_x = 0;
        for j = 1:N
            a_j = alpha_f(j);
            B_x = B_x + a_j * T0 * sum(w(:,j)) ...
                + a_j * s * T0  * n_flop * sum(z(:,j)) ...
                + a_j * (s * T0 / R) * sum(u(:,j));
        end

        maximize(A_x - theta * B_x)

        subject to
            %光纤分配约束
            for i = 1:M
                sum(x(i, :)) <= 1;
            end
             
            % 总功率限制
             sum(sum(w)) <= P_total;
             
            % 时延约束（用 f 近似表示 f_s / f_CPU）
            for j = 1:N
                T0 +  s * T0 * n_flop * f(j)  + ( s * T0 * f_s(j) ) / R <= T_max;

            end 

            % 每个设备功率限制
            for j = 1:N
                %P_static + C_S * f_s(j) + C_CPU / (100^3) * pow_p(t(j), -3) + P_tx(j) <= alpha_f(j) * sum(w(:, j));
                P_static + C_S * f_s(j) + C_CPU * pow_pos(f_CPU(j)/100, 3) + P_tx(j) <= oe_eta * alpha_f(j) * sum(w(:, j));
            end

            % 变量范围约束
            f_CPU_min <= f_CPU <= f_CPU_max;
            f_s_min <= f_s <= f_s_max;
            P_fiber_min <= P_fiber <= P_fiber_max;

            %1 / f_CPU_max <= t <= 1 / f_CPU_min;

            f_min =  f_s_min / f_CPU_max;
            f_max =  f_s_max / f_CPU_min;
            f_min <= f <= f_max;

            
            for j = 1:N


                f(j) <= f_s(j)/(f_CPU_min);
                % f(j) <= (f_s_max/f_CPU_min) * (f_CPU_max - f_CPU(j))/(f_CPU_max - f_CPU_min);
                % {2 * f(j), f_CPU(j) - f_s(j), f_CPU(j) + f_s(j)} == rotated_lorentz(3);
                % f(j) >= f_s(j)/(f_CPU_max);
                f(j) >= (f_s_max/f_CPU_max) * (f_CPU_max - f_CPU(j))/(f_CPU_max - f_CPU_min); % 避免下界过于保守
            end
            

            % Big-M + McCormick for w(i,j), z, u
            for i = 1:M
                for j = 1:N

                    % w(i,j) = x(i,j) * P_fiber(i)
                    w(i,j) >= P_fiber(i) - P_fiber_max * (1 - x(i,j));
                    w(i,j) <= P_fiber(i);
                    w(i,j) <= P_fiber_max * x(i,j);

                    % 放松 McCormick for t(i,j) = P_fiber(i) * f(j)

                    % t(i,j) >= P_fiber_min * f(j) + P_fiber(i) * f_min - P_fiber_min * f_min;
                    % t(i,j) >= P_fiber_max * f(j) + P_fiber(i) * f_max - P_fiber_max * f_max;
                    % t(i,j) <= P_fiber_max * f(j) + P_fiber(i) * f_min - P_fiber_max * f_min;
                    % t(i,j) <= P_fiber_min * f(j) + P_fiber(i) * f_max - P_fiber_min * f_max;


                    % z(i,j) = w(i,j) * f(j)
                    z(i,j) >= 0 * f(j) + w(i,j) * f_min - 0 * f_min;
                    z(i,j) >= P_fiber_max * f(j) + w(i,j) * f_max - P_fiber_max * f_max;
                    z(i,j) <= 0 * f(j) + w(i,j) * f_max - 0 * f_max;
                    z(i,j) <= P_fiber_max * f(j) + w(i,j) * f_min - P_fiber_max * f_min;
                    % z(i, j) >= w(i, j) * f_s_min / f_CPU_max;  
                    % z(i, j) <= w(i, j) * f_s_max / f_CPU_min;  

                   
                    % % 放宽 Big-M 范围
                    % t_ij_L = 0.001;
                    % t_ij_U = 0.008;
                    % z(i,j) <= t(i,j) - t_ij_L * (1 - x(i,j));
                    % z(i,j) >= t(i,j) - t_ij_U * (1 - x(i,j));
                    % z(i,j) <= t_ij_U * x(i,j);
                    % z(i,j) >= t_ij_L * x(i,j);

                    % McCormick: u(i,j) = w(i,j) * f_s(j)
                    u(i,j) >= P_fiber_min * f_s(j) + f_s_min * w(i,j) - P_fiber_min * f_s_min;
                    u(i,j) >= P_fiber_max * f_s(j) + f_s_max * w(i,j) - P_fiber_max * f_s_max;
                    u(i,j) <= P_fiber_min * f_s(j) + f_s_max * w(i,j) - P_fiber_min * f_s_max;
                    u(i,j) <= P_fiber_max * f_s(j) + f_s_min * w(i,j) - P_fiber_max * f_s_min;

                    % u(i,j) >= f_s_min * w(i,j);
                    % % u(i,j) >= P_fiber_max * f_s(j) + f_s_max * w(i,j)  - P_fiber_max * f_s_max;
                    % % u(i,j) <= P_fiber_max * f_s(j) + f_s_min * w(i,j)  - P_fiber_max * f_s_min;
                    % u(i,j) <= f_s_max * w(i,j);
                    % u(i,j) >= f_s(j) * P_fiber_min;
                    % u(i,j) <= f_s(j) * P_fiber_max;

                    
                end
            end
    cvx_end

    % Dinkelbach 更新
    if abs(A_x - theta * B_x) < epsilon || abs(B_x) < 1e-6
        break;
    else
        theta = A_x / max(B_x, 1e-6);
    end
end

%% 输出结果
disp('Dinkelbach 迭代收敛:');  
disp(['最终 theta: ', num2str(theta)]);  
disp(['迭代次数: ', num2str(iter)]);  
disp('光纤分配矩阵 x:');  
disp(x);  
disp('光纤功率 P_fiber:');  
disp(P_fiber);  
disp('CPU 频率 f_CPU:');  
disp(f_CPU);  
disp('采样频率 f_s:');  
disp(f_s);  
disp('最优目标值:');  
disp(cvx_optval);  
disp('最终的分子 A_x:');
disp(A_x);
disp('最终的分母 B_x:');
disp(B_x);
disp('最终 theta:');
disp(theta);
disp('手动计算的最优目标值:');
disp(A_x - theta * B_x);

disp('CVX 计算的最优目标值:');
disp(cvx_optval);

%% 计算误差
f_CPU_approx = f_s ./ f;  % 必须用点除
relative_error = (f_CPU - f_CPU_approx) ./ f_CPU;

disp('设备频率估算结果：');
disp(table(f_CPU', f_CPU_approx', relative_error', ...
    'VariableNames', {'理论值\n','估算值\n','相对误差\n'}));

% 检查每个设备的时延约束是否满足
disp('--- 时延约束检查 ---');
for j = 1:N
    delay = T0 + f(j) * s * T0 * n_flop + (f_s(j) * s * T0)/R;
    fprintf('设备 %d 时延: %.6f (限制 %.2f) --> %s\n', j, delay, T_max, ...
        ternary(delay <= T_max, '满足', '不满足'));
end

% 检查每个设备的功率约束是否满足
disp('--- 功率约束检查 ---');
for j = 1:N
    % power_lhs = P_static + C_S * f_s(j) + C_CPU / (100^3) * pow_p(t(j), -3) + P_tx(j);
    power_lhs = P_static + C_S * f_s(j) + C_CPU * (f_CPU(j)/100)^3 + P_tx(j);
    power_rhs = sum(full(w(:, j)) .* alpha_f(j)) * oe_eta;
    fprintf('设备 %d 功率: LHS = %.6f, RHS = %.6f --> %s\n', j, power_lhs, power_rhs, ...
        ternary(power_lhs <= power_rhs, '满足', '不满足'));
end

% 辅助函数
function out = ternary(condition, trueStr, falseStr)
    if condition
        out = trueStr;
    else
        out = falseStr;
    end
end