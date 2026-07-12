clear;clc;close all;

% ========== 参数定义 ==========
M = 10;         % 光纤数量
N = 7;          % 设备数量
L0 = 5;         % 基础距离 (km)
r = 1.5;        % 距离增量 (km/设备)
alpha_loss = 0.2;   % 光纤损耗 (dB/km)
P_total = 15;       % 总功率 (W)
T_max = 1.4;        % 最大时延 (s)
P_static = 0.1;     % 静态功耗 (W)
C_S = 0.5;          % 采样功耗系数 (W/MHz)
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
P_fiber_max = 2;    % 光纤最大功率 (W)
P_tx = 1e-4;         % 假设传输功耗 (W)
output_interval = 1;

rng(444) % 设定随机数种子

% ========== 遗传算法参数 ==========
pop_size = 5000;      % 种群大小，即有5000个解的组合。在这些解之间进行选择，搜索出最优或近似最优的解组合
max_gen =1000;      % 最大迭代次数
crossover_prob = 0.8; % 交叉概率
mutation_prob = 0.2; % 变异概率
penalty_coeff = 1e5; % 约束违反惩罚系数

% ========== 迭代数据记录 ==========
best_fitness_history = zeros(max_gen, 1);
avg_fitness_history = zeros(max_gen, 1);

% ========== 初始化种群 ==========
% 初始化种群结构体，每个个体包含x_ij、P_fiber、f_s、f_CPU四个字段
% 结构体数组population初始为空，后通过循环向内填入pop_size个体
population = struct('x_ij', {}, 'P_fiber', {}, 'f_s', {}, 'f_CPU', {});

% 遍历种群中每个个体
for ind = 1:pop_size
    % 初始化每个个体
    
    % 初始化x_ij（保证每个光纤最多连接一个设备，但不保证每个设备至少被连接）
    x_ij = zeros(M, N);
    for i = 1:M
        if rand() < 0.5  % 50%概率连接设备
            j = randi(N); %在N个设备中随机选一个
            x_ij(i, j) = 1; %第i行第j列设为1
        end
    end
    
    % 初始化其他变量
    % 随机初始化光纤功率，在[P_fiber_min, P_fiber_max]范围内
    P_fiber = P_fiber_min + (P_fiber_max - P_fiber_min)*rand(M, 1); % 随机生成一个M*1的列向量
    
    % 随机初始化采样频率，在[f_s_min, f_s_max]范围内
    f_s = f_s_min + (f_s_max - f_s_min)*rand(N, 1);
    
    % 随机初始化CPU频率，在[f_CPU_min, f_CPU_max]范围内
    f_CPU = f_CPU_min + (f_CPU_max - f_CPU_min)*rand(N, 1);
    
    % 构造结构体个体并加入种群 作为population(ind)中的第ind个个体
    population(ind).x_ij = x_ij;
    population(ind).P_fiber = P_fiber;
    population(ind).f_s = f_s;
    population(ind).f_CPU = f_CPU;
end

% ========== 主循环 ==========
best_fitness = -inf;
best_solution = [];
for gen = 1:max_gen
    % 迭代主循环
    
    % 计算适应度
    fitness = zeros(pop_size, 1); % 存放每个个体的适应度，初始化全为0 
    for ind = 1:pop_size
        sol = population(ind); % 取出第ind个个体sol（结构体，包含该个体的四个变量）
        
        % 调用适应度计算函数，对个体进行评估
        [f, v1, v2, v3] = calculateFitness(sol, M, N, L0, r, alpha_loss,...
            P_static, C_S, C_CPU, s, T0, n_flop, R, P_total, T_max,...
            P_tx, penalty_coeff);
        % 计算适应度值和约束违反量
        fitness(ind) = f - penalty_coeff*(v1 + v2 + v3); % 最终适应值=原始目标函数-违反约束惩罚
        % 适应度函数，考虑约束惩罚。如果个体违反某个约束，适应度显著降低，在选择阶段就被淘汰
    end
    
    % 更新最优解
    [max_fit, idx] = max(fitness);
    if max_fit > best_fitness
        best_fitness = max_fit;
        best_solution = population(idx);
    end
        
    % 记录迭代数据
    best_fitness_history(gen) = max_fit;
    avg_fitness_history(gen) = mean(fitness);

    % 记录迭代数据
    best_fitness_history(gen) = max_fit;
    avg_fitness_history(gen) = mean(fitness);

    % === 收敛判断模块 ===
    if gen > 50
        delta = abs(diff(best_fitness_history(gen-49:gen)));
        std_recent = std(best_fitness_history(gen-49:gen));
        avg_gap = mean(best_fitness_history(gen-49:gen) - avg_fitness_history(gen-49:gen));

        if std_recent < 1e-4 && mean(delta(end-10:end)) < 1e-4 && avg_gap < 1e-3
            fprintf('提前终止：第 %d 代后收敛稳定（std = %.2e, gap = %.2e）\n', ...
                gen, std_recent, avg_gap);
            break;
        end
    end

    
    % ========== 迭代进度输出 ==========
    if mod(gen, output_interval) == 0 || gen == 1 %每代结果输出
        fprintf('Generation %4d:  Best Fitness=%.4f  |  Avg Fitness=%.4f\n',...
                gen, best_fitness_history(gen), avg_fitness_history(gen));
    end
    
    % 选择（锦标赛选择）
    % 创建一个空结构体数组存放新的下一代种群
    new_pop = struct('x_ij', {}, 'P_fiber', {}, 'f_s', {}, 'f_CPU', {});
    for i = 1:pop_size %共生成pop_size个新个体
        % 随机选择两个个体进行比较
        candidates = randperm(pop_size, 2); %随机选出2个不同个体作为候选者
        if fitness(candidates(1)) > fitness(candidates(2)) %比较两者的适应度，选择适应度更高的个体作为优胜者
            new_pop(i) = population(candidates(1)); % 选出的更优个体复制到下一代种群
        else
            new_pop(i) = population(candidates(2)); 
        end
    end
    
    % 交叉操作
    for i = 1:2:pop_size-1 %取出两个相邻个体进行交叉
        if rand() < crossover_prob %通过概率选择是否交叉，防止交叉太频繁
            % 选择父代
            p1 = new_pop(i);
            p2 = new_pop(i+1);
            
            % 单点交叉
            cross_point = randi(M*N + M + 2*N - 1); % 每个个体按变量拼接成一个一维长向量，即x_ij,P_fiber,f_s,f_cpu合并
            % 生成随机交叉点
            
            % 调用自定义的crossover函数，执行交叉并生成子代
            c1 = crossover(p1, p2, cross_point, M, N);
            c2 = crossover(p2, p1, cross_point, M, N);
            new_pop(i) = c1;
            new_pop(i+1) = c2;
        end
    end
    
    % 变异操作 防止陷入局部最优
    for i = 1:pop_size %对种群内每个个体依次进行变异
        if rand() < mutation_prob %以设定的变异概率为依据，决定是否对该个体进行变异
            mutated = mutate(new_pop(i), M, N, P_fiber_min, P_fiber_max,...
                f_s_min, f_s_max, f_CPU_min, f_CPU_max); %对第i个个体进行变异操作，变异的范围由变量的上下限控制
            new_pop(i) = mutated; % 变异后的个体替换原个体
        end
    end
    
    population = new_pop;
end

% ========== 结果输出 ==========
disp('========== Optimal Solution ==========');
% 1. 显示连接矩阵
fprintf('\n[光纤-设备连接矩阵 x_ij]:\n');
disp(best_solution.x_ij);


% 2. 显示光纤功率
fprintf('\n[光纤发射功率 P_fiber_i (W)]:\n');
for i = 1:M
    fprintf('Fiber %2d: %.4f\n', i, best_solution.P_fiber(i));
end

% 3. 显示设备采样频率
fprintf('\n[设备采样频率 f_s_j (MHz)]:\n');
for j = 1:N
    fprintf('Device %d: %.2f\n', j, best_solution.f_s(j));
end

% 4. 显示CPU频率
fprintf('\n[CPU频率 f_CPU_j (MHz)]:\n');
for j = 1:N
    fprintf('Device %d: %.2f\n', j, best_solution.f_CPU(j));
end

% 5. 性能指标
fprintf('\n[性能指标]\n');
fprintf('最佳适应度: %.4f\n', best_fitness);
fprintf('总消耗功率: %.2f W (限额: %d W)\n', ...
    sum(sum(best_solution.x_ij .* best_solution.P_fiber)), P_total);

% ========== 绘制迭代曲线 ==========
% ==== 可视化：完整曲线 vs 收敛放大视角（仅上下容差线） ==== 
figure;

% ==== 子图1：全局趋势 ====
subplot(1,2,1);
plot(1:max_gen, best_fitness_history, 'b-', 'LineWidth', 1.5); hold on;
plot(1:max_gen, avg_fitness_history, 'r--', 'LineWidth', 1.5);
title('Global Fitness Evolution');
xlabel('Generation'); ylabel('Fitness');
legend('Best','Average','Location','southeast');
grid on; set(gca, 'FontSize', 11);


% ==== 子图2：收敛区域放大（收敛线 + 容差线 + 百分比标注） ====
subplot(1,2,2);
plot(1:max_gen, best_fitness_history, 'b-', 'LineWidth', 1.5); hold on;

% 1. 提取最后150代最常出现的值作为收敛值
tail_vals = best_fitness_history(end-149:end);
[unique_vals, ~, idx_map] = unique(round(tail_vals, 6)); % 避免浮点误差
counts = accumarray(idx_map, 1);
[~, max_idx] = max(counts);
converged_val = unique_vals(max_idx);

% 2. 计算容差和百分比差
gen_687_val = best_fitness_history(687);
tolerance = abs(gen_687_val - converged_val);
percent_diff = 100 * tolerance / converged_val;

% 3. 画收敛线和容差线
% yline(converged_val, 'g-', 'LineWidth', 1.5, 'DisplayName', 'Converged');
% yline(converged_val + tolerance, 'r--', 'LineWidth', 1, 'DisplayName', '+Tolerance');
% yline(converged_val - tolerance, 'r--', 'LineWidth', 1, 'DisplayName', '-Tolerance');

% 4. 标注百分比差距（收敛线和上容差线之间）
% x_txt = max_gen - 100;  % 标注横坐标位置
% y_txt = converged_val + tolerance/2;  % 标注纵坐标位置（中间）
% label_str = sprintf('%.2f%%', percent_diff);
% text(x_txt, y_txt, label_str, 'Color', 'm', 'FontSize', 11, 'FontWeight', 'bold');

% 5. 视觉调整
ylim([converged_val - 1.5*tolerance, converged_val + 1.5*tolerance]);
title('Zoomed View: Convergence Band with Percent Gap');
xlabel('Generation'); ylabel('Fitness');
legend('Best','Converged','+/- Tolerance','Location','southeast');
grid on; set(gca, 'FontSize', 11);

% figure;
% hold on;
% plot(1:max_gen, best_fitness_history, 'b-', 'LineWidth', 1.5);
% plot(1:max_gen, avg_fitness_history, 'r--', 'LineWidth', 1.5);
% xlabel('Generation');
% ylabel('Fitness');
% title('Evolution of Fitness');
% legend('Best Fitness', 'Average Fitness', 'Location', 'southeast');
% grid on;
% set(gca, 'FontSize', 12);

% ========== 适应度计算函数 ==========
function [f, v1, v2, v3] = calculateFitness(sol, M, N, L0, r, alpha,...
    P_static, C_S, C_CPU, s, T0, n_flop, R, P_total, T_max,...
    P_tx,penalty)
    x_ij = sol.x_ij;
    P_fiber = sol.P_fiber;
    f_s = sol.f_s;
    f_CPU = sol.f_CPU;
    
    % 计算目标函数分子
    numerator = sum(f_s * s * T0);
    
    % 初始化分母和约束违反量
    denominator = 0;
    v1 = 0; % 约束1违反
    v2 = 0; % 约束2违反
    v3 = 0; % 约束3违反
    
    % 对每个设备计算
    for j = 1:N
        Lj = L0 + (2*j-1)*r;
        % 计算第j个设备到光纤的距离
        
        fiber_power = sum(x_ij(:,j) .* P_fiber .* 10.^(-alpha*Lj/10));
        
       
        % 计算光纤功率（考虑损耗）
        % right1 = P_static + C_S*f_s(j) + C_CPU*(f_CPU(j)/1e3)^3 + P_tx;
        right1 = P_static + C_S*f_s(j) + C_CPU*(f_CPU(j)/100)^3 + P_tx;
        % 计算设备功率需求
        
        % 约束1检查（光纤功率 >= 设备功率需求）
        if fiber_power < right1
            v1 = v1 + (right1 - fiber_power);
        end
        
        % 分母计算
        term_j = fiber_power*(T0 + (f_s(j)*s*T0*n_flop)/f_CPU(j)...
            + (f_s(j)*s*T0)/R);
        denominator = denominator + term_j;
        
        % 约束3检查（时延 <= T_max）
        delay = T0 + (f_s(j)*s*T0*n_flop)/f_CPU(j) + (f_s(j)*s*T0)/R;
        if delay > T_max
            v3 = v3 + (delay - T_max);
        end
    end
    
    % 约束2检查（总功率 <= P_total）
    total_power = sum(sum(x_ij .* P_fiber));
    if total_power > P_total
        v2 = total_power - P_total;
    end
    
    % 计算最终适应度
    if denominator == 0
        f = 0;
    else
        f = numerator / denominator;
    end
end


% ========== 交叉操作 ==========
function child = crossover(p1, p2, point, M, N)
    % 将解结构转换为向量
    vec1 = [p1.x_ij(:); p1.P_fiber; p1.f_s; p1.f_CPU];
    vec2 = [p2.x_ij(:); p2.P_fiber; p2.f_s; p2.f_CPU];
    
    % 执行单点交叉
    child_vec = vec1;
    child_vec(point+1:end) = vec2(point+1:end);
    
    % 转换回结构体，从child_vec的不同区段中提取，分别恢复决策变量
    child.x_ij = reshape(child_vec(1:M*N), M, N);
    child.P_fiber = child_vec(M*N+1:M*N+M);
    child.f_s = child_vec(M*N+M+1:M*N+M+N);
    child.f_CPU = child_vec(M*N+M+N+1:end);
    
    % 修复x_ij约束（每光纤最多连接一个设备）
    for i = 1:M
        idx = find(child.x_ij(i,:) == 1);
        if numel(idx) > 1
            % 如果有多个连接，违反"每个光纤最多连接一个设备"，随机保留一个
            child.x_ij(i,:) = 0;
            child.x_ij(i, idx(randi(numel(idx)))) = 1;
        end
    end
end

% ========== 变异操作 ==========
function mutated = mutate(sol, M, N, Pf_min, Pf_max, fs_min, fs_max, fCPU_min, fCPU_max)
    mutated = sol; % 复制原始个体
    
    % 变异x_ij（随机翻转连接）
    for i = 1:M
        if rand() < 0.1
            j = randi(N);
            mutated.x_ij(i,:) = 0;
            mutated.x_ij(i,j) = 1;
        end
    end
    
    % 变异P_fiber（随机选择一个光纤并重新生成其功率）
    mut_idx = randi(M);
    mutated.P_fiber(mut_idx) = Pf_min + (Pf_max - Pf_min)*rand();
    
    % 变异f_s（随机选择一个设备并重新生成其采样频率）
    mut_idx = randi(N);
    mutated.f_s(mut_idx) = fs_min + (fs_max - fs_min)*rand();
    
    % 变异f_CPU（随机选择一个设备并重新生成其CPU频率）
    mut_idx = randi(N);
    mutated.f_CPU(mut_idx) = fCPU_min + (fCPU_max - fCPU_min)*rand();
end