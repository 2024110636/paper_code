clear; clc; close all;

% ========== 参数定义 ==========
M = 10; 
N = 8;
L0 = 5; 
r = 1.5; 
alpha_loss = 0.2;
P_total = 15; 
T_max = 1.4;
P_static = 0.1; 
C_S = 0.5; 
C_CPU = 0.001;
s = 2; 
T0 = 1; 
n_flop = 100; 
R = 500;
f_CPU = 600; % 固定值 MHz
f_s_min = 1; f_s_max = 2;
P_fiber_min = 0; P_fiber_max = 2;
output_interval = 1;


% ========== 遗传算法参数 ==========
pop_size = 5000; 
max_gen = 1000;
crossover_prob = 0.8; 
mutation_prob = 0.3;
penalty_coeff = 1e5;
rng(444)  % 随机种子

% ========== 数据记录 ==========
best_fitness_history = zeros(max_gen, 1);
avg_fitness_history = zeros(max_gen, 1);


%% 初始化种群
population = struct('x_ij', {}, 'P_fiber', {}, 'f_s', {}); % 不用优化fCPU
for ind = 1:pop_size
    x_ij = zeros(M, N); 
    used_fibers = [];

    for j = 1:6
        available = setdiff(1:M, used_fibers);
        i = available(randi(length(available)));
        x_ij(i, j) = 1; used_fibers(end+1) = i;
    end
     
    for j = 7:8
        available = setdiff(1:M, used_fibers);
        chosen = randsample(available, 2);
        for k = 1:2
            x_ij(chosen(k), j) = 1;
            used_fibers(end+1) = chosen(k);
        end
    end

    % for j = 3:4
    %     available = setdiff(1:M, used_fibers);
    %     chosen = randsample(available, 3);
    %     for k = 1:3
    %         x_ij(chosen(k), j) = 1;
    %         used_fibers(end+1) = chosen(k);
    %     end
    % end

    
    % 初始化其他变量
    P_fiber = P_fiber_min + (P_fiber_max - P_fiber_min)*rand(M, 1);
    f_s = f_s_min + (f_s_max - f_s_min)*rand(N, 1);
    
    population(ind).x_ij = x_ij;
    population(ind).P_fiber = P_fiber;
    population(ind).f_s = f_s;
end

%% ========== 遗传主循环 ==========
best_fitness = -inf; 
best_solution = [];

for gen = 1:max_gen

    % 计算适应度
    fitness = zeros(pop_size, 1);

    for ind = 1:pop_size
        % 计算适应度值和约束违反量
        sol = population(ind);
        [f, v1, v2, v3] = calculateFitness(sol, M, N, L0, r, alpha_loss,...
            P_static, C_S, C_CPU, s, T0, n_flop, R, P_total, T_max, penalty_coeff, f_CPU);
        
        % 适应度函数，考虑约束惩罚
        fitness(ind) = f - penalty_coeff*(v1 + v2 + v3);
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

    
    % ========== 迭代进度输出 ==========
    if mod(gen, output_interval) == 0 || gen == 1
        fprintf('Generation %4d:  Best Fitness=%.4f  |  Avg Fitness=%.4f\n',...
                gen, best_fitness_history(gen), avg_fitness_history(gen));
    end
    
    
    % 选择
    new_pop = struct('x_ij', {}, 'P_fiber', {}, 'f_s', {});  % 不用选择fCPU
    
    for i = 1:pop_size
        candidates = randperm(pop_size, 2);
        if fitness(candidates(1)) > fitness(candidates(2))
            new_pop(i) = population(candidates(1));
        else
            new_pop(i) = population(candidates(2));
        end
    end

    % 交叉
    for i = 1:2:pop_size-1
        
        if rand() < crossover_prob
            % 选择父代
            c1 = crossover(new_pop(i), new_pop(i+1), randi(M*N + M + N - 1), M, N);
            c2 = crossover(new_pop(i+1), new_pop(i), randi(M*N + M + N - 1), M, N);
            new_pop(i) = c1; 
            new_pop(i+1) = c2;
        end
    end

    % 变异
    for i = 1:pop_size
        if rand() < mutation_prob
            new_pop(i) = mutate(new_pop(i), M, N, P_fiber_min, P_fiber_max, f_s_min, f_s_max);
        end
    end

    population = new_pop;

end

%% ========== 可视化输出 ==========

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
fprintf('\n[CPU频率 f_CPU (MHz)]:\n');
for j = 1:N
    fprintf('Device %d: %.2f\n', j, f_CPU);
end

% 5. 性能指标
fprintf('\n[性能指标]\n');
fprintf('最佳适应度: %.4f\n', best_fitness);
fprintf('总消耗功率: %.2f W (限额: %d W)\n', ...
    sum(sum(best_solution.x_ij .* best_solution.P_fiber)), P_total);

% 计算每个设备的能耗与时延
device_energy = zeros(N, 1);
device_delay  = zeros(N, 1);

for j = 1:N
    Lj = L0 + (2*j-1)*r;
    fiber_power = sum(best_solution.x_ij(:,j) .* best_solution.P_fiber .* 10.^(-alpha_loss*Lj/10));
    P_tx_j = 0.2 * 10^((2 + 0.2 * Lj) / 10) / 1000;
    required = P_static + C_S*best_solution.f_s(j) + ...
        C_CPU*(f_CPU/100)^3 + P_tx_j;

    device_energy(j) = required; % 记录总能耗
    
    device_delay(j) = T0 + (best_solution.f_s(j)*s*T0*n_flop)/f_CPU + ...
        (best_solution.f_s(j)*s*T0)/R; % 记录时延
end

disp("---每台设备的时延---");
for j = 1:N
    fprintf('设备 %d 时延: %.6f s (限制 %.2f s)\n', j, device_delay(j), T_max);
   
end

disp("---每台设备的能耗---");
for j = 1:N
    fprintf('设备 %d 时延: %.6f W \n', j, device_energy(j));
end


%% ========== 主图（完整适应度曲线） ==========

% ==== 可视化：完整曲线 vs 收敛放大视角（仅上下容差线） ====
% ========== 主图 ==========
figure;

% 主图（只显示 fitness > 0）
plot(1:max_gen, best_fitness_history, 'b-', 'LineWidth', 1.5); hold on;
% plot(1:max_gen, avg_fitness_history, 'r--', 'LineWidth', 1.5);
xlabel('Generation'); ylabel('Fitness');
title('Fitness Evolution');
legend('Best Fitness', 'Average', 'Location', 'southeast');
grid on;
set(gca, 'FontSize', 12);
ylim([0, max(best_fitness_history) * 1.05]);  % 只显示正数部分

% ========== 子图（550~1000代局部） ==========
sub_start = 550;
sub_end = max_gen;
sub_idx = sub_start:sub_end;
sub_vals = best_fitness_history(sub_idx);

% 计算局部最大值 & 出现频率最多的值（收敛值）
[unique_vals, ~, idx_map] = unique(round(sub_vals, 6));
counts = accumarray(idx_map, 1);
[~, idx_max] = max(counts);
val_converged = unique_vals(idx_max);  % 收敛值
val_peak = max(sub_vals);              % 局部最大值
gap_ratio = (val_peak - val_converged) / val_converged;  % 百分差

% 容差范围（±gap_ratio%）
eps_tol = gap_ratio * val_converged;

% 插入子图
inset_pos = [0.5, 0.4, 0.35, 0.30];  % 可微调位置
inset_axes = axes('Position', inset_pos);
box on; grid on; hold on;

plot(sub_idx, sub_vals, 'k-', 'LineWidth', 1.2);
yline(val_converged, 'g-.', 'LineWidth', 1.2);
yline(val_converged + eps_tol, 'r--', 'LineWidth', 1);
% yline(val_converged - eps_tol, 'r--', 'LineWidth', 1);

% 添加注释
text(sub_start+5, val_converged + eps_tol + 0.001, ...
    sprintf('%.2f%% deviation', gap_ratio*100), ...
    'FontSize', 12, 'Color', 'm');

xlabel('Generation', 'FontSize', 16);
ylabel('Fitness', 'FontSize', 16);

ylim([val_converged - 0.01, val_peak + 0.01]);  % 自动适配
set(gca, 'FontSize', 12);




%% ========== 函数：适应度函数 ==========
function [f, v1, v2, v3] = calculateFitness(sol, M, N, L0, r, alpha, P_static, C_S, C_CPU, s, T0, n_flop, R, P_total, T_max, penalty, f_CPU)
x_ij = sol.x_ij; 
P_fiber = sol.P_fiber; 
f_s = sol.f_s;

numerator = sum(f_s * s * T0); 
denominator = 0; 
v1 = 0; 
v2 = 0; 
v3 = 0;


for j = 1:N
    Lj = L0 + (2*j-1)*r;

    fiber_power = sum(x_ij(:,j) .* P_fiber .* 10.^(-alpha*Lj/10));
    
    P_tx_j = 0.2 * 10^((2 + 0.2 * Lj) / 10) / 1000;
    
    required = P_static + C_S*f_s(j) + C_CPU*(f_CPU/100)^3 + P_tx_j;
    
    if fiber_power < required
        v1 = v1 + (required - fiber_power);
    end
    
    delay = T0 + (f_s(j)*s*T0*n_flop)/f_CPU + (f_s(j)*s*T0)/R;
    
    if delay > T_max
        v3 = v3 + (delay - T_max);
    end
    
    denominator = denominator + fiber_power * delay;
end

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


% ========== 交叉函数 ==========
function child = crossover(p1, p2, point, M, N)

vec1 = [p1.x_ij(:); p1.P_fiber; p1.f_s];
vec2 = [p2.x_ij(:); p2.P_fiber; p2.f_s];

child_vec = vec1; 
child_vec(point+1:end) = vec2(point+1:end);

child.x_ij = reshape(child_vec(1:M*N), M, N);
child.P_fiber = child_vec(M*N+1:M*N+M);
child.f_s = child_vec(M*N+M+1:end);

% 修复光纤最多连一个设备，每设备连接数约束
for i = 1:M
    idx = find(child.x_ij(i,:) == 1);

    if numel(idx) > 1
        child.x_ij(i,:) = 0;
        child.x_ij(i, idx(randi(numel(idx)))) = 1;
    end
end

for j = 1:N
    % t = 1;
    t = (j <= 2)*1 + (j > 2)*2;
    % t = (j <= 2)*2 + (j > 2)*3;

    idx = find(child.x_ij(:, j) == 1);

    if numel(idx) ~= t
        child.x_ij(:, j) = 0;
        available = find(sum(child.x_ij, 2) == 0);

        if length(available) < t
            available = 1:M;
        end

        chosen = randsample(available, t);
        child.x_ij(chosen, j) = 1;

    end
end
end


% ========== 变异函数 ==========
function mutated = mutate(sol, M, N, Pf_min, Pf_max, fs_min, fs_max)
mutated = sol;

i = randi(M);
mutated.P_fiber(i) = Pf_min + (Pf_max - Pf_min)*rand();

j = randi(N);
mutated.f_s(j) = fs_min + (fs_max - fs_min)*rand();

for i = 1:M
    idx = find(mutated.x_ij(i,:) == 1);

    if numel(idx) > 1
        mutated.x_ij(i,:) = 0;
        mutated.x_ij(i, idx(randi(numel(idx)))) = 1;
    end

end

for j = 1:N
    % t = 1;
    t = (j <= 2)*1 + (j > 2)*2;
    % t = (j <= 2)*2 + (j > 2)*3;

    idx = find(mutated.x_ij(:, j) == 1);
    
    if numel(idx) ~= t
        mutated.x_ij(:, j) = 0;
        available = find(sum(mutated.x_ij, 2) == 0);

        if length(available) < t
            available = 1:M;
        end

        chosen = randsample(available, t);
        mutated.x_ij(chosen, j) = 1;

    end
end
end