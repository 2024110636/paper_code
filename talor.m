%% SCP + Dinkelbach
scp_max_iter = 10;    % SCP最大迭代
scp_tol = 1e-3;       % SCP收敛阈值
scp_iter = 0;
converged = false;

%% 参数设置
M = 8         % 光纤数量
N = 5;         % 设备数量
L0 = 5;         % 基础距离 (km)
r = 1.5;        % 每个设备增加距离 (km)
alpha_loss = 0.2;   % 光纤损耗 (dB/km)
P_total = 14;       % 总供电功率 (W)
T_max = 1.3;        % 最大允许时延 (s)

% 设备功率模型
P_static = 0.1;    % 静态功率 (W)
C_S = 0.5;            % 采样功率系数 (W/MHz)
C_CPU = 0.001;       % CPU功率系数 (W/GHz^3)
s = 2;              % 每个样本大小 (bit)
T0 = 1;             % 采样时间 (s)
n_flop = 100;       % 每个样本浮点操作数
R = 500;            % 光纤速率 (Mbit/s)

% CPU & 采样频率 (单位：Hz)
f_CPU_min = 500; f_CPU_max = 1000;   % MHz
f_s_min = 1;      f_s_max = 2;        % MHz

% 光纤功率约束
P_fiber_min = 0; P_fiber_max = 2;     % W

% 光纤损耗系数 alpha_f(j)
L = @(j) L0 + r * (2 * j - 1);
alpha_f = arrayfun(@(j) 10^(-alpha_loss * L(j)/10), 1:N);

% 发射功率 P_tx(j)
P_tx = arrayfun(@(j) 0.2 * 10^((2 + 0.2 * L(j)) / 10) / 1000, 1:N);

% 初始化泰勒展开的参考点
f_s_bar = (f_s_min + f_s_max)/2 * ones(N,1);
f_CPU_bar = (f_CPU_min + f_CPU_max)/2 * ones(N,1);

while ~converged && scp_iter < scp_max_iter
    scp_iter = scp_iter + 1;

    fprintf('\n=== SCP 外层迭代 %d ===\n', scp_iter);

    theta = 0;
    epsilon = 5e-4;
    dinkelbach_iter = 0;
    max_inner_iter = 50;

    while dinkelbach_iter < max_inner_iter
        dinkelbach_iter = dinkelbach_iter + 1;

        cvx_begin
            cvx_solver mosek
            cvx_precision low
            cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', 300);
            cvx_solver_settings('MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-3);
            cvx_expert true;

            % 决策变量
            variable x(M,N) binary
            variable P_fiber(M) nonnegative
            variable f_CPU(N) nonnegative
            variable f_s(N) nonnegative

            variable f(N) nonnegative
            variable w(M,N) nonnegative
            variable z(M,N) nonnegative
            variable u(M,N) nonnegative

            % 目标
            A_x = s * T0 * sum(f_s);

            B_x = 0;
            for j=1:N
                a_j = alpha_f(j);
                B_x = B_x + a_j * T0 * sum(w(:,j)) ...
                + a_j * s * T0  * n_flop * sum(z(:,j)) ...
                + a_j * (s * T0 / R) * sum(u(:,j));
            end

            maximize(A_x - theta*B_x)

            subject to
                % 光纤分配
                for i=1:M
                    sum(x(i,:)) <= 1;
                end

                % 总功率限制
                sum(sum(w)) <= P_total;

                % 时延
                for j=1:N
                    T0 + s*T0*n_flop*f(j) + (s*T0*f_s(j))/R <= T_max;
                end

                % 功率
                for j=1:N
                    P_static + C_S * f_s(j) + C_CPU * pow_pos(f_CPU(j)/100, 3) + P_tx(j) <= alpha_f(j) * sum(w(:, j));
                end

                % 范围
                f_CPU_min <= f_CPU <= f_CPU_max;
                f_s_min <= f_s <= f_s_max;
                P_fiber_min <= P_fiber <= P_fiber_max;

                % 泰勒一阶线性化
                for j=1:N
                    f(j) == ( f_s_bar(j)/f_CPU_bar(j) ) ...
                          + (1/f_CPU_bar(j)) * ( f_s(j) - f_s_bar(j) ) ...
                          - ( f_s_bar(j)/(f_CPU_bar(j)^2) ) * ( f_CPU(j) - f_CPU_bar(j) );
                end

                % Big-M 和 McCormick
                for i=1:M
                    for j=1:N
                        w(i,j) >= P_fiber(i) - P_fiber_max*(1-x(i,j));
                        w(i,j) <= P_fiber(i);
                        w(i,j) <= P_fiber_max*x(i,j);

                        z(i,j) >= 0*f(j) + w(i,j)*f_min - 0*f_min;
                        z(i,j) >= P_fiber_max*f(j) + w(i,j)*f_max - P_fiber_max*f_max;
                        z(i,j) <= 0*f(j) + w(i,j)*f_max - 0*f_max;
                        z(i,j) <= P_fiber_max*f(j) + w(i,j)*f_min - P_fiber_max*f_min;

                        u(i,j) >= P_fiber_min*f_s(j) + f_s_min*w(i,j) - P_fiber_min*f_s_min;
                        u(i,j) >= P_fiber_max*f_s(j) + f_s_max*w(i,j) - P_fiber_max*f_s_max;
                        u(i,j) <= P_fiber_min*f_s(j) + f_s_max*w(i,j) - P_fiber_min*f_s_max;
                        u(i,j) <= P_fiber_max*f_s(j) + f_s_min*w(i,j) - P_fiber_max*f_s_min;
                    end
                end

        cvx_end

        % Dinkelbach 更新
        if abs(A_x - theta*B_x) < epsilon || abs(B_x) < 1e-6
            break
        else
            theta = A_x / max(B_x, 1e-6);
        end

    end % Dinkelbach

    % SCP收敛性检查
    diff_s = norm(f_s - f_s_bar) / (norm(f_s_bar) + 1e-6);
    diff_CPU = norm(f_CPU - f_CPU_bar) / (norm(f_CPU_bar) + 1e-6);

    fprintf('SCP误差：采样 %.6f，CPU %.6f\n', diff_s, diff_CPU);

    if diff_s < scp_tol && diff_CPU < scp_tol
        converged = true;
    else
        % 更新线性化点
        f_s_bar = f_s;
        f_CPU_bar = f_CPU;
    end
end

%% 结果
disp('=== SCP + Dinkelbach 迭代收敛 ===');
disp(['总SCP迭代: ', num2str(scp_iter)]);
disp(['最后theta: ', num2str(theta)]);
disp('最终光纤分配 x:'); disp(x);
disp('最终采样频率 f_s:'); disp(f_s);
disp('最终CPU频率 f_CPU:'); disp(f_CPU);
