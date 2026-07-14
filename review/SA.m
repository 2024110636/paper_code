clear;clc;close all;

% ========== 参数设置（与standard_GA.m一致） ==========
M = 10;         % 光纤数量
N = 7;          % 设备数量
L0 = 5;         % 基础距离 (km)
r = 1.5;        % 距离增量 (km/设备)
alpha_loss = 0.2;   % 光纤损耗 (dB/km)
P_total = 25;       % 总功率 (W)
T_max = 1.4;        % 最大时延 (s)
P_static = 0.1;     % 静态功耗 (W)
C_S = 0.2;          % 采样功耗系数 (W/MHz)
C_CPU = 0.001;      % CPU功耗系数 (W/GHz^3)
s = 2;              % 样本大小 (bit)
T0 = 1;             % 采样时间 (s)
n_flop = 100;       % 浮点操作数/样本
R = 500;            % 传输速率 (Mbit/s)
f_CPU_min = 500;    % CPU最小频率 (MHz)
f_CPU_max = 1000;   % CPU最大频率 (MHz)
f_s_min = 1;        % 采样最小频率 (MHz)
f_s_max = 2;        % 采样最大频率 (MHz)
P_fiber_min = 0;    % 光纤最小功率 (W)
P_fiber_max = 3;    % 光纤最大功率 (W)

rng(444) % 固定随机种子

% ========== SA算法参数 ==========
T_init = 1e3;       % 初始温度
T_min = 1e-4;       % 终止温度
alpha = 0.999;      % 冷却系数（几何降温）
L_inner = 50;       % 每个温度下的内循环迭代次数
penalty_coeff = 1e3; % 约束违反惩罚系数
output_interval = 100; % 每隔多少代输出一次

% ========== 数据记录 ==========
best_fitness_history = [];
current_fitness_history = [];
temperature_history = [];
iter_count = 0;

% ========== 初始化当前解 ==========
% x_ij：每根光纤50%概率分配给一个随机设备
x_ij = zeros(M, N);
for i = 1:M
    if rand() < 0.5
        j = randi(N);
        x_ij(i, j) = 1;
    end
end

% 连续变量均匀初始化
P_fiber = P_fiber_min + (P_fiber_max - P_fiber_min)*rand(M, 1);
f_s = f_s_min + (f_s_max - f_s_min)*rand(N, 1);
f_CPU = f_CPU_min + (f_CPU_max - f_CPU_min)*rand(N, 1);

% 计算初始适应度
current_fitness = evaluateFitness(x_ij, P_fiber, f_s, f_CPU, ...
    M, N, L0, r, alpha_loss, P_static, C_S, C_CPU, s, T0, n_flop, ...
    R, P_total, T_max, penalty_coeff);

% 记录全局最优
best_x_ij = x_ij;
best_P_fiber = P_fiber;
best_f_s = f_s;
best_f_CPU = f_CPU;
best_fitness = current_fitness;

% ========== SA主循环 ==========
T = T_init;

while T > T_min
    for k = 1:L_inner
        iter_count = iter_count + 1;

        % 生成邻域解（在当前解基础上随机扰动一个维度）
        new_x_ij = x_ij;
        new_P_fiber = P_fiber;
        new_f_s = f_s;
        new_f_CPU = f_CPU;

        % 随机选择扰动类型（4类变量等概率）
        perturb_type = randi(4);

        switch perturb_type
            case 1  % 扰动x_ij：随机选一根光纤重新分配
                i = randi(M);
                j = randi(N);
                new_x_ij(i,:) = 0;
                new_x_ij(i,j) = 1;

            case 2  % 扰动P_fiber：高斯扰动一根光纤
                i = randi(M);
                delta = 0.4 * (P_fiber_max - P_fiber_min) * randn();
                new_P_fiber(i) = P_fiber(i) + delta;
                new_P_fiber(i) = max(P_fiber_min, min(P_fiber_max, new_P_fiber(i)));

            case 3  % 扰动f_s：高斯扰动一个设备
                j = randi(N);
                delta = 0.4 * (f_s_max - f_s_min) * randn();
                new_f_s(j) = f_s(j) + delta;
                new_f_s(j) = max(f_s_min, min(f_s_max, new_f_s(j)));

            case 4  % 扰动f_CPU：高斯扰动一个设备
                j = randi(N);
                delta = 0.4 * (f_CPU_max - f_CPU_min) * randn();
                new_f_CPU(j) = f_CPU(j) + delta;
                new_f_CPU(j) = max(f_CPU_min, min(f_CPU_max, new_f_CPU(j)));
        end

        % 计算新解适应度
        new_fitness = evaluateFitness(new_x_ij, new_P_fiber, new_f_s, new_f_CPU, ...
            M, N, L0, r, alpha_loss, P_static, C_S, C_CPU, s, T0, n_flop, ...
            R, P_total, T_max, penalty_coeff);

        % Metropolis接受准则
        delta_f = new_fitness - current_fitness;
        if delta_f > 0 || rand() < exp(delta_f / T)
            % 接受新解
            x_ij = new_x_ij;
            P_fiber = new_P_fiber;
            f_s = new_f_s;
            f_CPU = new_f_CPU;
            current_fitness = new_fitness;

            % 更新全局最优
            if current_fitness > best_fitness
                best_fitness = current_fitness;
                best_x_ij = x_ij;
                best_P_fiber = P_fiber;
                best_f_s = f_s;
                best_f_CPU = f_CPU;
            end
        end

        % 记录
        best_fitness_history(iter_count) = best_fitness;
        current_fitness_history(iter_count) = current_fitness;
        temperature_history(iter_count) = T;
    end

    % 降温
    T = T * alpha;

    % 进度输出
    if mod(iter_count, output_interval * L_inner) == 0
        fprintf('T=%.4e | Iter=%d | Best=%.4f | Current=%.4f\n', ...
            T, iter_count, best_fitness, current_fitness);
    end
end

fprintf('\nSA终止：最终温度=%.4e，总迭代=%d\n', T, iter_count);

% ========== 结果输出 ==========
disp('========== Optimal Solution ==========');

fprintf('\n[光纤-设备连接矩阵 x_ij]:\n');
disp(best_x_ij);

fprintf('\n[光纤发射功率 P_fiber_i (W)]:\n');
for i = 1:M
    fprintf('Fiber %2d: %.4f\n', i, best_P_fiber(i));
end

fprintf('\n[设备采样频率 f_s_j (MHz)]:\n');
for j = 1:N
    fprintf('Device %d: %.2f\n', j, best_f_s(j));
end

fprintf('\n[CPU频率 f_CPU_j (MHz)]:\n');
for j = 1:N
    fprintf('Device %d: %.2f\n', j, best_f_CPU(j));
end

fprintf('\n[性能指标]\n');
fprintf('最佳适应度: %.4f\n', best_fitness);
fprintf('总消耗功率: %.2f W (限额: %d W)\n', ...
    sum(sum(best_x_ij .* best_P_fiber)), P_total);

% 每台设备的详细信息
for j = 1:N
    Lj = L0 + (2*j-1)*r;
    fiber_power = sum(best_x_ij(:,j) .* best_P_fiber .* 10.^(-alpha_loss*Lj/10));
    electrical_power = fiber_power * 0.3;
    P_tx_j = 10^((0.2 .* Lj + 5)/10) / 1000;
    required = P_static + C_S*best_f_s(j) + ...
        C_CPU*(best_f_CPU(j)/100)^3 + P_tx_j;
    delay_j = T0 + (best_f_s(j)*s*T0*n_flop)/best_f_CPU(j) + ...
        (best_f_s(j)*s*T0)/R;
    fprintf('设备 %d: 距离=%.1fkm, 接收光功率=%.4fW, 电功率=%.4fW, 需求功率=%.4fW, 时延=%.4fs\n', ...
        j, Lj, fiber_power, electrical_power, required, delay_j);
end

% ========== 收敛曲线可视化 ==========
figure;

% 左图：全局适应度曲线
subplot(1,2,1);
plot(1:iter_count, best_fitness_history, 'b-', 'LineWidth', 1.5); hold on;
plot(1:iter_count, current_fitness_history, 'r-', 'LineWidth', 0.5);
title('SA Fitness Evolution');
xlabel('Iteration'); ylabel('Fitness');
legend('Best','Current','Location','southeast');
grid on; set(gca, 'FontSize', 11);

% 右图：温度曲线
subplot(1,2,2);
semilogy(1:iter_count, temperature_history, 'k-', 'LineWidth', 1);
title('Temperature Schedule');
xlabel('Iteration'); ylabel('Temperature');
grid on; set(gca, 'FontSize', 11);

% ========== 适应度计算函数 ==========
function f = evaluateFitness(x_ij, P_fiber, f_s, f_CPU, ...
    M, N, L0, r, alpha, P_static, C_S, C_CPU, s, T0, n_flop, ...
    R, P_total, T_max, penalty)

    % 目标函数分子
    numerator = sum(f_s * s * T0);

    % 分母及约束违反量
    denominator = 0;
    v1 = 0; v2 = 0; v3 = 0;

    for j = 1:N
        Lj = L0 + (2*j-1)*r;

        % 设备j接收到的光纤功率（含衰减）
        fiber_power = sum(x_ij(:,j) .* P_fiber .* 10.^(-alpha*Lj/10));

        % 光电转换后的电功率
        eta_oe = 0.3;
        electrical_power = fiber_power * eta_oe;

        % 距离相关的发射功率 P_tx(j)
        P_tx_j = 0.2 * 10^((2 + 0.2 * Lj) / 10) / 1000;

        % 设备j的功耗需求
        required = P_static + C_S*f_s(j) + C_CPU*(f_CPU(j)/100)^3 + P_tx_j;

        % 约束1：光电转换后的供给功率 >= 设备需求功率
        if electrical_power < required
            v1 = v1 + (required - electrical_power);
        end

        % 分母累加项（使用电功率）
        term_j = fiber_power*(T0 + (f_s(j)*s*T0*n_flop)/f_CPU(j)...
            + (f_s(j)*s*T0)/R);
        denominator = denominator + term_j;

        % 约束3：时延 <= T_max
        delay = T0 + (f_s(j)*s*T0*n_flop)/f_CPU(j) + (f_s(j)*s*T0)/R;
        if delay > T_max
            v3 = v3 + (delay - T_max);
        end
    end

    % 约束2：总功率 <= P_total
    total_power = sum(sum(x_ij .* P_fiber));
    if total_power > P_total
        v2 = total_power - P_total;
    end

    % 计算适应度（能效 = 数据量 / 能耗）
    if denominator == 0
        f = -penalty * (v1 + v2 + v3);
    else
        f = numerator / denominator - penalty * (v1 + v2 + v3);
    end
end
