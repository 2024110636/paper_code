%% 启发式算法对比实验：Standard GA vs GA-SFDR vs SA
% 三种部署场景：近端分布、远端分布、两端分布
% N=5 设备，M=10 光纤
% 统一参数：P_total=25, P_fiber_max=3, C_S=0.2, η_oe=0.3

clear; clc; close all;

%% ========== 公共参数 ==========
M = 10;
N = 5;
alpha_loss = 0.2;
P_total = 15;
T_max = 1.4;
P_static = 0.1;
C_S = 0.2;
C_CPU = 0.001;
s = 2;
T0 = 1;
n_flop = 100;
R = 500;
f_CPU_min = 500; f_CPU_max = 1000;
f_s_min = 1;     f_s_max = 2;
P_fiber_min = 0; P_fiber_max = 2;
eta_oe = 0.3;
penalty_coeff = 1e3;

%% ========== 启发式算法参数（精简，缩短运行时间） ==========
pop_size = 2000;
max_gen = 500;
crossover_prob = 0.8;
mutation_prob = 0.2;

% SA参数
T_init = 1e3;
T_min = 1e-4;
alpha_cool = 0.999;
L_inner = 50;

%% ========== 三种部署场景 ==========
rng(42);  L1 = sort(5 + 10 * rand(1, N));   % 近端: [5, 15]
rng(43);  L2 = sort(15 + 15 * rand(1, N));  % 远端: [15, 30]
rng(44);  L3 = sort(5 + 25 * rand(1, N));   % 两端: [5, 30]

scenarios = {L1, L2, L3};
scenario_names = {'Near-end', 'Far-end', 'Both-ends'};

%% ========== 结果存储 ==========
EE_results = zeros(3, 3);  % 行=场景, 列=算法(StdGA, GA-SFDR, SA)

%% ========== 主实验循环 ==========
for sc = 1:3
    L = scenarios{sc};
    alpha_f = 10.^(-alpha_loss * L / 10);
    P_tx_arr = 10.^((5 + 0.2 * L) / 10) / 1000;

    fprintf('\n========================================\n');
    fprintf('Scenario %d: %s\n', sc, scenario_names{sc});
    fprintf('L = ['); fprintf('%.1f ', L); fprintf('] km\n');
    fprintf('========================================\n');

    %% ---------- 算法1: Standard GA ----------
    rng(444);
    [ee_ga, ~] = runStandardGA(M, N, alpha_f, P_tx_arr, ...
        P_total, T_max, P_static, C_S, C_CPU, s, T0, n_flop, R, ...
        f_CPU_min, f_CPU_max, f_s_min, f_s_max, P_fiber_min, P_fiber_max, ...
        eta_oe, penalty_coeff, pop_size, max_gen, crossover_prob, mutation_prob);
    fprintf('  Standard GA: EE = %.4f\n', ee_ga);
    EE_results(sc, 1) = ee_ga;

    %% ---------- 算法2: GA-SFDR ----------
    rng(444);
    % 根据距离自适应分配：近端设备1根光纤，远端设备2根
    [~, sort_idx] = sort(L);
    n_near = 3;  % 前3近设备各1根
    n_far = N - n_near;  % 后2远设备各2根
    fiber_per_device = zeros(1, N);
    fiber_per_device(sort_idx(1:n_near)) = 1;
    fiber_per_device(sort_idx(n_near+1:end)) = 2;

    [ee_sfdr, ~] = runGASFDR(M, N, alpha_f, P_tx_arr, fiber_per_device, ...
        P_total, T_max, P_static, C_S, C_CPU, s, T0, n_flop, R, ...
        f_CPU_min, f_CPU_max, f_s_min, f_s_max, P_fiber_min, P_fiber_max, ...
        eta_oe, penalty_coeff, pop_size, max_gen, crossover_prob, mutation_prob);
    fprintf('  GA-SFDR:     EE = %.4f\n', ee_sfdr);
    EE_results(sc, 2) = ee_sfdr;

    %% ---------- 算法3: SA ----------
    rng(444);
    [ee_sa, ~] = runSA(M, N, alpha_f, P_tx_arr, ...
        P_total, T_max, P_static, C_S, C_CPU, s, T0, n_flop, R, ...
        f_CPU_min, f_CPU_max, f_s_min, f_s_max, P_fiber_min, P_fiber_max, ...
        eta_oe, penalty_coeff, T_init, T_min, alpha_cool, L_inner);
    fprintf('  SA:          EE = %.4f\n', ee_sa);
    EE_results(sc, 3) = ee_sa;
end

%% ========== 汇总结果 ==========
fprintf('\n\n========================================\n');
fprintf('       Heuristic Algorithm Comparison\n');
fprintf('========================================\n');
fprintf('%-12s  %10s  %10s  %10s\n', 'Scenario', 'StdGA', 'GA-SFDR', 'SA');
fprintf('----------------------------------------\n');
for sc = 1:3
    fprintf('%-12s  %10.4f  %10.4f  %10.4f\n', ...
        scenario_names{sc}, EE_results(sc,1), EE_results(sc,2), EE_results(sc,3));
end
fprintf('----------------------------------------\n');
fprintf('%-12s  %10.4f  %10.4f  %10.4f\n', 'Average', ...
    mean(EE_results(:,1)), mean(EE_results(:,2)), mean(EE_results(:,3)));

%% ========== 可视化 ==========
figure('Name', 'Heuristic Algorithm Comparison', 'Position', [200, 200, 700, 450]);

bar_data = EE_results;
b = bar(bar_data);
b(1).FaceColor = [0.56 0.80 0.56]; b(1).EdgeColor = 'k'; b(1).LineWidth = 0.8;
b(2).FaceColor = [0.96 0.69 0.52]; b(2).EdgeColor = 'k'; b(2).LineWidth = 0.8;
b(3).FaceColor = [0.61 0.76 0.91]; b(3).EdgeColor = 'k'; b(3).LineWidth = 0.8;

set(gca, 'XTickLabel', scenario_names, 'FontSize', 12);
xlabel('Deployment Scenarios');
ylabel('Energy Efficiency (Mbit/J)');
title('Heuristic Algorithm Comparison under Non-uniform Deployment');
legend('Standard GA', 'GA-SFDR', 'SA', 'Location', 'best');

hold on;
for i = 1:3
    for j = 1:3
        if bar_data(i,j) > 0
            text(j + 0.22*(i-2), bar_data(i,j) + max(bar_data(:))*0.02, ...
                sprintf('%.4f', bar_data(i,j)), ...
                'HorizontalAlignment', 'center', 'FontSize', 9);
        end
    end
end
hold off;
grid off;

%% ========================================================================
%%                        函 数 定 义
%% ========================================================================

%% ========== Standard GA ==========
function [best_ee, best_sol] = runStandardGA(M, N, alpha_f, P_tx_arr, ...
    P_total, T_max, P_static, C_S, C_CPU, s, T0, n_flop, R, ...
    f_CPU_min, f_CPU_max, f_s_min, f_s_max, P_fiber_min, P_fiber_max, ...
    eta_oe, penalty_coeff, pop_size, max_gen, crossover_prob, mutation_prob)

    % 初始化种群
    pop = struct('x_ij', {}, 'P_fiber', {}, 'f_s', {}, 'f_CPU', {});
    for ind = 1:pop_size
        x_ij = zeros(M, N);
        for i = 1:M
            if rand() < 0.5
                j = randi(N);
                x_ij(i, j) = 1;
            end
        end
        pop(ind).x_ij = x_ij;
        pop(ind).P_fiber = P_fiber_min + (P_fiber_max - P_fiber_min)*rand(M, 1);
        pop(ind).f_s = f_s_min + (f_s_max - f_s_min)*rand(N, 1);
        pop(ind).f_CPU = f_CPU_min + (f_CPU_max - f_CPU_min)*rand(N, 1);
    end

    best_fitness = -inf; best_sol = [];

    for gen = 1:max_gen
        fitness = zeros(pop_size, 1);
        for ind = 1:pop_size
            sol = pop(ind);
            [f, pen] = calcFitness(sol.x_ij, sol.P_fiber, sol.f_s, sol.f_CPU, ...
                M, N, alpha_f, P_tx_arr, P_total, T_max, P_static, C_S, C_CPU, ...
                s, T0, n_flop, R, eta_oe, penalty_coeff);
            fitness(ind) = f - pen;
        end
        [max_fit, idx] = max(fitness);
        if max_fit > best_fitness
            best_fitness = max_fit;
            best_sol = pop(idx);
        end

        % 选择
        new_pop = struct('x_ij', {}, 'P_fiber', {}, 'f_s', {}, 'f_CPU', {});
        for i = 1:pop_size
            c = randperm(pop_size, 2);
            if fitness(c(1)) > fitness(c(2))
                new_pop(i) = pop(c(1));
            else
                new_pop(i) = pop(c(2));
            end
        end
        % 交叉
        for i = 1:2:pop_size-1
            if rand() < crossover_prob
                pt = randi(M*N + M + 2*N - 1);
                c1 = doCrossover(new_pop(i), new_pop(i+1), pt, M, N);
                c2 = doCrossover(new_pop(i+1), new_pop(i), pt, M, N);
                new_pop(i) = c1; new_pop(i+1) = c2;
            end
        end
        % 变异
        for i = 1:pop_size
            if rand() < mutation_prob
                new_pop(i) = doMutate(new_pop(i), M, N, P_fiber_min, P_fiber_max, ...
                    f_s_min, f_s_max, f_CPU_min, f_CPU_max);
            end
        end
        pop = new_pop;
    end

    % 可行性检查
    feasible = true;
    for j = 1:N
        fp = sum(best_sol.x_ij(:,j) .* best_sol.P_fiber .* alpha_f(j));
        ep = fp * eta_oe;
        req = P_static + C_S*best_sol.f_s(j) + C_CPU*(best_sol.f_CPU(j)/100)^3 + P_tx_arr(j);
        if ep < req || fp < 1e-10
            feasible = false;
            break;
        end
    end
    if ~feasible
        best_ee = 0;
    else
        [best_ee, ~] = calcFitness(best_sol.x_ij, best_sol.P_fiber, best_sol.f_s, best_sol.f_CPU, ...
            M, N, alpha_f, P_tx_arr, P_total, T_max, P_static, C_S, C_CPU, ...
            s, T0, n_flop, R, eta_oe, 0);
    end
end

%% ========== GA-SFDR ==========
function [best_ee, best_sol] = runGASFDR(M, N, alpha_f, P_tx_arr, fiber_per_device, ...
    P_total, T_max, P_static, C_S, C_CPU, s, T0, n_flop, R, ...
    f_CPU_min, f_CPU_max, f_s_min, f_s_max, P_fiber_min, P_fiber_max, ...
    eta_oe, penalty_coeff, pop_size, max_gen, crossover_prob, mutation_prob)

    % 初始化种群（结构化分配）
    pop = struct('x_ij', {}, 'P_fiber', {}, 'f_s', {}, 'f_CPU', {});
    for ind = 1:pop_size
        x_ij = zeros(M, N);
        used = [];
        for j = 1:N
            t = fiber_per_device(j);
            avail = setdiff(1:M, used);
            if length(avail) < t
                avail = 1:M;
            end
            chosen = randsample(avail, t);
            for k = 1:t
                x_ij(chosen(k), j) = 1;
                used(end+1) = chosen(k);
            end
        end
        pop(ind).x_ij = x_ij;
        pop(ind).P_fiber = P_fiber_min + (P_fiber_max - P_fiber_min)*rand(M, 1);
        pop(ind).f_s = f_s_min + (f_s_max - f_s_min)*rand(N, 1);
        pop(ind).f_CPU = f_CPU_min + (f_CPU_max - f_CPU_min)*rand(N, 1);
    end

    best_fitness = -inf; best_sol = [];

    for gen = 1:max_gen
        fitness = zeros(pop_size, 1);
        for ind = 1:pop_size
            sol = pop(ind);
            [f, pen] = calcFitness(sol.x_ij, sol.P_fiber, sol.f_s, sol.f_CPU, ...
                M, N, alpha_f, P_tx_arr, P_total, T_max, P_static, C_S, C_CPU, ...
                s, T0, n_flop, R, eta_oe, penalty_coeff);
            fitness(ind) = f - pen;
        end
        [max_fit, idx] = max(fitness);
        if max_fit > best_fitness
            best_fitness = max_fit;
            best_sol = pop(idx);
        end

        % 选择
        new_pop = struct('x_ij', {}, 'P_fiber', {}, 'f_s', {}, 'f_CPU', {});
        for i = 1:pop_size
            c = randperm(pop_size, 2);
            if fitness(c(1)) > fitness(c(2))
                new_pop(i) = pop(c(1));
            else
                new_pop(i) = pop(c(2));
            end
        end
        % 交叉 + 结构化修复
        for i = 1:2:pop_size-1
            if rand() < crossover_prob
                pt = randi(M*N + M + 2*N - 1);
                c1 = doCrossoverSFDR(new_pop(i), new_pop(i+1), pt, M, N, fiber_per_device);
                c2 = doCrossoverSFDR(new_pop(i+1), new_pop(i), pt, M, N, fiber_per_device);
                new_pop(i) = c1; new_pop(i+1) = c2;
            end
        end
        % 变异 + 结构化修复 + 交换变异
        for i = 1:pop_size
            if rand() < mutation_prob
                new_pop(i) = doMutateSFDR(new_pop(i), M, N, P_fiber_min, P_fiber_max, ...
                    f_s_min, f_s_max, f_CPU_min, f_CPU_max, fiber_per_device);
            end
        end
        pop = new_pop;
    end

    % 可行性检查
    feasible = true;
    for j = 1:N
        fp = sum(best_sol.x_ij(:,j) .* best_sol.P_fiber .* alpha_f(j));
        ep = fp * eta_oe;
        req = P_static + C_S*best_sol.f_s(j) + C_CPU*(best_sol.f_CPU(j)/100)^3 + P_tx_arr(j);
        if ep < req || fp < 1e-10
            feasible = false;
            break;
        end
    end
    if ~feasible
        best_ee = 0;
    else
        [best_ee, ~] = calcFitness(best_sol.x_ij, best_sol.P_fiber, best_sol.f_s, best_sol.f_CPU, ...
            M, N, alpha_f, P_tx_arr, P_total, T_max, P_static, C_S, C_CPU, ...
            s, T0, n_flop, R, eta_oe, 0);
    end
end

%% ========== SA ==========
function [best_ee, best_sol] = runSA(M, N, alpha_f, P_tx_arr, ...
    P_total, T_max, P_static, C_S, C_CPU, s, T0, n_flop, R, ...
    f_CPU_min, f_CPU_max, f_s_min, f_s_max, P_fiber_min, P_fiber_max, ...
    eta_oe, penalty_coeff, T_init, T_min, alpha_cool, L_inner)

    % 初始化
    x_ij = zeros(M, N);
    for i = 1:M
        if rand() < 0.5
            j = randi(N); x_ij(i, j) = 1;
        end
    end
    P_fiber = P_fiber_min + (P_fiber_max - P_fiber_min)*rand(M, 1);
    f_s = f_s_min + (f_s_max - f_s_min)*rand(N, 1);
    f_CPU = f_CPU_min + (f_CPU_max - f_CPU_min)*rand(N, 1);

    [cur_f, cur_pen] = calcFitness(x_ij, P_fiber, f_s, f_CPU, ...
        M, N, alpha_f, P_tx_arr, P_total, T_max, P_static, C_S, C_CPU, ...
        s, T0, n_flop, R, eta_oe, penalty_coeff);
    cur_fitness = cur_f - cur_pen;

    best_x = x_ij; best_P = P_fiber; best_fs = f_s; best_fc = f_CPU;
    best_fitness = cur_fitness;

    T = T_init;
    while T > T_min
        for k = 1:L_inner
            new_x = x_ij; new_P = P_fiber; new_fs = f_s; new_fc = f_CPU;
            pt = randi(4);
            switch pt
                case 1
                    i = randi(M); j = randi(N);
                    new_x(i,:) = 0; new_x(i,j) = 1;
                case 2
                    i = randi(M);
                    d = 0.3*(P_fiber_max-P_fiber_min)*randn();
                    new_P(i) = max(P_fiber_min, min(P_fiber_max, P_fiber(i)+d));
                case 3
                    j = randi(N);
                    d = 0.3*(f_s_max-f_s_min)*randn();
                    new_fs(j) = max(f_s_min, min(f_s_max, f_s(j)+d));
                case 4
                    j = randi(N);
                    d = 0.3*(f_CPU_max-f_CPU_min)*randn();
                    new_fc(j) = max(f_CPU_min, min(f_CPU_max, f_CPU(j)+d));
            end
            [nf, np] = calcFitness(new_x, new_P, new_fs, new_fc, ...
                M, N, alpha_f, P_tx_arr, P_total, T_max, P_static, C_S, C_CPU, ...
                s, T0, n_flop, R, eta_oe, penalty_coeff);
            new_fitness = nf - np;

            df = new_fitness - cur_fitness;
            if df > 0 || rand() < exp(df / T)
                x_ij = new_x; P_fiber = new_P; f_s = new_fs; f_CPU = new_fc;
                cur_fitness = new_fitness;
                if cur_fitness > best_fitness
                    best_fitness = cur_fitness;
                    best_x = x_ij; best_P = P_fiber; best_fs = f_s; best_fc = f_CPU;
                end
            end
        end
        T = T * alpha_cool;
    end

    best_sol.x_ij = best_x; best_sol.P_fiber = best_P;
    best_sol.f_s = best_fs; best_sol.f_CPU = best_fc;

    % 可行性检查
    feasible = true;
    for j = 1:N
        fp = sum(best_x(:,j) .* best_P .* alpha_f(j));
        ep = fp * eta_oe;
        req = P_static + C_S*best_fs(j) + C_CPU*(best_fc(j)/100)^3 + P_tx_arr(j);
        if ep < req || fp < 1e-10
            feasible = false;
            break;
        end
    end
    if ~feasible
        best_ee = 0;
    else
        [best_ee, ~] = calcFitness(best_x, best_P, best_fs, best_fc, ...
            M, N, alpha_f, P_tx_arr, P_total, T_max, P_static, C_S, C_CPU, ...
            s, T0, n_flop, R, eta_oe, 0);
    end
end

%% ========== 统一适应度函数 ==========
function [f, penalty_val] = calcFitness(x_ij, P_fiber, f_s, f_CPU, ...
    M, N, alpha_f, P_tx_arr, P_total, T_max, P_static, C_S, C_CPU, ...
    s, T0, n_flop, R, eta_oe, penalty_coeff)

    numerator = sum(f_s * s * T0);
    denominator = 0;
    v1 = 0; v2 = 0; v3 = 0;

    for j = 1:N
        fiber_power = sum(x_ij(:,j) .* P_fiber .* alpha_f(j));
        electrical_power = fiber_power * eta_oe;
        required = P_static + C_S*f_s(j) + C_CPU*(f_CPU(j)/100)^3 + P_tx_arr(j);

        if electrical_power < required
            v1 = v1 + (required - electrical_power);
        end

        term_j = fiber_power * (T0 + (f_s(j)*s*T0*n_flop)/f_CPU(j) + (f_s(j)*s*T0)/R);
        denominator = denominator + term_j;

        delay = T0 + (f_s(j)*s*T0*n_flop)/f_CPU(j) + (f_s(j)*s*T0)/R;
        if delay > T_max
            v3 = v3 + (delay - T_max);
        end
    end

    total_power = sum(sum(x_ij .* P_fiber));
    if total_power > P_total
        v2 = total_power - P_total;
    end

    if denominator == 0
        f = 0;
    else
        f = numerator / denominator;
    end

    penalty_val = penalty_coeff * (v1 + v2 + v3);
end

%% ========== 交叉（Standard GA） ==========
function child = doCrossover(p1, p2, point, M, N)
    vec1 = [p1.x_ij(:); p1.P_fiber; p1.f_s; p1.f_CPU];
    vec2 = [p2.x_ij(:); p2.P_fiber; p2.f_s; p2.f_CPU];
    cv = vec1; cv(point+1:end) = vec2(point+1:end);
    child.x_ij = reshape(cv(1:M*N), M, N);
    child.P_fiber = cv(M*N+1:M*N+M);
    child.f_s = cv(M*N+M+1:M*N+M+N);
    child.f_CPU = cv(M*N+M+N+1:end);
    for i = 1:M
        idx = find(child.x_ij(i,:) == 1);
        if numel(idx) > 1
            child.x_ij(i,:) = 0;
            child.x_ij(i, idx(randi(numel(idx)))) = 1;
        end
    end
end

%% ========== 变异（Standard GA） ==========
function m = doMutate(sol, M, N, Pf_min, Pf_max, fs_min, fs_max, fc_min, fc_max)
    m = sol;
    for i = 1:M
        if rand() < 0.1
            j = randi(N); m.x_ij(i,:) = 0; m.x_ij(i,j) = 1;
        end
    end
    i = randi(M); m.P_fiber(i) = Pf_min + (Pf_max-Pf_min)*rand();
    j = randi(N); m.f_s(j) = fs_min + (fs_max-fs_min)*rand();
    j = randi(N); m.f_CPU(j) = fc_min + (fc_max-fc_min)*rand();
end

%% ========== 交叉+结构化修复（GA-SFDR） ==========
function child = doCrossoverSFDR(p1, p2, point, M, N, fpd)
    vec1 = [p1.x_ij(:); p1.P_fiber; p1.f_s; p1.f_CPU];
    vec2 = [p2.x_ij(:); p2.P_fiber; p2.f_s; p2.f_CPU];
    cv = vec1; cv(point+1:end) = vec2(point+1:end);
    child.x_ij = reshape(cv(1:M*N), M, N);
    child.P_fiber = cv(M*N+1:M*N+M);
    child.f_s = cv(M*N+M+1:M*N+M+N);
    child.f_CPU = cv(M*N+M+N+1:end);
    % 修复：每根光纤最多连一个设备
    for i = 1:M
        idx = find(child.x_ij(i,:) == 1);
        if numel(idx) > 1
            child.x_ij(i,:) = 0;
            child.x_ij(i, idx(randi(numel(idx)))) = 1;
        end
    end
    % 修复：每设备连接数约束
    for j = 1:N
        t = fpd(j);
        idx = find(child.x_ij(:,j) == 1);
        if numel(idx) ~= t
            child.x_ij(:,j) = 0;
            avail = find(sum(child.x_ij, 2) == 0);
            if length(avail) < t; avail = 1:M; end
            chosen = randsample(avail, t);
            child.x_ij(chosen, j) = 1;
        end
    end
end

%% ========== 变异+结构化修复+交换变异（GA-SFDR） ==========
function m = doMutateSFDR(sol, M, N, Pf_min, Pf_max, fs_min, fs_max, fc_min, fc_max, fpd)
    m = sol;
    i = randi(M); m.P_fiber(i) = Pf_min + (Pf_max-Pf_min)*rand();
    j = randi(N); m.f_s(j) = fs_min + (fs_max-fs_min)*rand();
    j = randi(N); m.f_CPU(j) = fc_min + (fc_max-fc_min)*rand();
    % 修复x_ij
    for i = 1:M
        idx = find(m.x_ij(i,:) == 1);
        if numel(idx) > 1
            m.x_ij(i,:) = 0;
            m.x_ij(i, idx(randi(numel(idx)))) = 1;
        end
    end
    for j = 1:N
        t = fpd(j);
        idx = find(m.x_ij(:,j) == 1);
        if numel(idx) ~= t
            m.x_ij(:,j) = 0;
            avail = find(sum(m.x_ij, 2) == 0);
            if length(avail) < t; avail = 1:M; end
            chosen = randsample(avail, t);
            m.x_ij(chosen, j) = 1;
        end
    end
    % 结构化交换变异
    if rand() < 0.3
        i1 = randi(M); i2 = randi(M);
        while i2 == i1; i2 = randi(M); end
        j1 = find(m.x_ij(i1,:) == 1);
        j2 = find(m.x_ij(i2,:) == 1);
        if length(j1) == 1 && length(j2) == 1 && j1 ~= j2
            m.x_ij(i1,:) = 0; m.x_ij(i2,:) = 0;
            m.x_ij(i1, j2) = 1; m.x_ij(i2, j1) = 1;
        end
    end
end
