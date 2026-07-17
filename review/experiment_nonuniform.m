%% R2-3 非均匀部署场景对比实验
% 三种算法对比：EE-Opt / FixCPU / AvgPower
% N=5 设备，3 种部署场景
% 单位：功率 W，时间 s，频率 MHz，速率 Mbit/s

clear all;
cvx_clear;
cvx_solver mosek;
cvx_precision low;
cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', 120);
cvx_solver_settings('MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-3);
cvx_expert true;

%% ========== 公共参数 ==========
M = 10;                    % 光纤数量
N = 5;                    % 设备数量
alpha_loss = 0.2;         % 光纤损耗 (dB/km)
P_total = 25;             % 总供电功率 (W)
T_max = 1.4;              % 最大允许时延 (s)

P_static = 0.1;           % 静态功率 (W)
C_S = 0.2;                % 采样功率系数 (W/MHz)
C_CPU = 0.001;            % CPU功率系数 (W/GHz^3)
s = 2;                    % 每个样本大小 (bit)
T0 = 1;                   % 采样时间 (s)
n_flop = 100;             % 每个样本浮点操作数
R = 500;                  % 光纤速率 (Mbit/s)

f_CPU_min = 500; f_CPU_max = 1000;   % MHz
f_s_min = 1;      f_s_max = 2;        % MHz
P_fiber_min = 0;  P_fiber_max = 3;    % W

%% ========== 三组非均匀部署场景 ==========
%scenarios = {
%    [6, 7, 10, 12, 15],       ...  % 场景1：近端聚集
 %   [15, 18, 20, 24, 30],      ...  % 场景2：远端聚集
%    [5, 9, 16, 24, 30]              % 场景3：两端分散
%};
%scenario_names = {'近端聚集', '远端聚集', '两端分散'};

rng(42);  L1 = sort(5 + 10 * rand(1, N));   % 近端: [5, 15]
rng(43);  L2 = sort(15 + 15 * rand(1, N));  % 远端: [15, 30]
rng(44);  L3 = sort(5 + 25 * rand(1, N));   % 两端: [5, 30]

scenarios = {L1, L2, L3};
scenario_names = {'近端分布', '远端分布', '两端分布'};

%% ========== 结果存储 ==========
EE_results = zeros(3, 3);   % 行=场景, 列=算法(EE-Opt, FixCPU, AvgPower)

%% ========== 主实验循环 ==========
for sc = 1:3
    L = scenarios{sc};
    alpha_f = 10.^(-alpha_loss * L / 10);
    P_tx = 10.^((5 + 0.2 * L) / 10) / 1000;
    % = 10^((0.2·L_j + 5)/10) / 1000 

    fprintf('\n========================================\n');
    fprintf('场景 %d: %s\n', sc, scenario_names{sc});
    fprintf('距离 L = ['); fprintf('%.0f ', L); fprintf('] km\n');
    fprintf('========================================\n');

    %% ---------- 算法1: EE-Opt (Dinkelbach + MOSEK) ----------
    theta = 0;
    max_iter = 50;
    epsilon = 1e-4;
    ee_opt = 0;

    for iter = 1:max_iter
        cvx_begin quiet
            cvx_solver mosek
            variable x(M, N) binary
            variable P_fiber(M) nonnegative
            variable f_CPU(N) nonnegative
            variable f_s(N) nonnegative
            variable f(N) nonnegative
            variable w(M, N) nonnegative
            variable z(M, N) nonnegative
            variable u(M, N) nonnegative

            A_x = s * T0 * sum(f_s);
            B_x = 0;
            for j = 1:N
                a_j = alpha_f(j);
                B_x = B_x + a_j * T0 * sum(w(:,j)) ...
                    + a_j * s * T0 * n_flop * sum(z(:,j)) ...
                    + a_j * (s * T0 / R) * sum(u(:,j));
            end

            maximize(A_x - theta * B_x)
            subject to
                for i = 1:M
                    sum(x(i, :)) <= 1;
                end
                sum(sum(w)) <= P_total;

                for j = 1:N
                    T0 + s * T0 * n_flop * f(j) + (s * T0 * f_s(j)) / R <= T_max;
                    P_static + C_S * f_s(j) + C_CPU * pow_pos(f_CPU(j)/100, 3) + P_tx(j) <= alpha_f(j) * sum(w(:, j)) * 0.3;
                end

                f_CPU_min <= f_CPU <= f_CPU_max;
                f_s_min <= f_s <= f_s_max;
                P_fiber_min <= P_fiber <= P_fiber_max;

                f_min_val = f_s_min / f_CPU_max;
                f_max_val = f_s_max / f_CPU_min;
                f_min_val <= f <= f_max_val;

                for j = 1:N
                    f(j) <= f_s(j) / f_CPU_min;
                    f(j) >= (f_s_max/f_CPU_max) * (f_CPU_max - f_CPU(j))/(f_CPU_max - f_CPU_min);
                end

                for i = 1:M
                    for j = 1:N
                        w(i,j) >= P_fiber(i) - P_fiber_max * (1 - x(i,j));
                        w(i,j) <= P_fiber(i);
                        w(i,j) <= P_fiber_max * x(i,j);

                        z(i,j) >= 0 * f(j) + w(i,j) * f_min_val - 0 * f_min_val;
                        z(i,j) >= P_fiber_max * f(j) + w(i,j) * f_max_val - P_fiber_max * f_max_val;
                        z(i,j) <= 0 * f(j) + w(i,j) * f_max_val - 0 * f_max_val;
                        z(i,j) <= P_fiber_max * f(j) + w(i,j) * f_min_val - P_fiber_max * f_min_val;

                        u(i,j) >= P_fiber_min * f_s(j) + f_s_min * w(i,j) - P_fiber_min * f_s_min;
                        u(i,j) >= P_fiber_max * f_s(j) + f_s_max * w(i,j) - P_fiber_max * f_s_max;
                        u(i,j) <= P_fiber_min * f_s(j) + f_s_max * w(i,j) - P_fiber_min * f_s_max;
                        u(i,j) <= P_fiber_max * f_s(j) + f_s_min * w(i,j) - P_fiber_max * f_s_min;
                    end
                end
        cvx_end

        if strcmp(cvx_status, 'Infeasible') || strcmp(cvx_status, 'Failed')
            fprintf('  EE-Opt: 问题无解 (iter=%d, status=%s)\n', iter, cvx_status);
            break;
        end

        if abs(A_x - theta * B_x) < epsilon
            ee_opt = A_x / max(B_x, 1e-10);
            fprintf('  EE-Opt: EE = %.4f (迭代 %d 次收敛)\n', ee_opt, iter);
            fprintf('  f_CPU = '); disp(f_CPU');
            fprintf('  f_s   = '); disp(f_s');
            fprintf('  x 分配:\n'); disp(full(x));
            fprintf('  P_fiber = '); disp(P_fiber');
            break;
        end
        theta = A_x / max(B_x, 1e-10);

        if iter == max_iter
            ee_opt = A_x / max(B_x, 1e-10);
            fprintf('  EE-Opt: EE = %.4f (达到最大迭代 %d 次)\n', ee_opt, max_iter);
        end
    end
    EE_results(sc, 1) = ee_opt;

    %% ---------- 算法2: FixCPU (固定f_CPU=600, 优化f_s和功率分配) ----------
    % 参照 fixed_fCPU.m: f_CPU固定后 f=f_s/f_CPU 是f_s的线性函数
    % 不需要 f 和 z 辅助变量, 只需 u = w * f_s 的 McCormick
    f_CPU_fix = 600;   % 固定CPU频率 (MHz)
    ee_fix = 0;

    theta_fix = 0;
    for iter = 1:max_iter
        cvx_begin quiet
            cvx_solver mosek
            variable x2(M, N) binary
            variable P_fiber2(M) nonnegative
            variable f_s2(N) nonnegative
            variable w2(M, N) nonnegative
            variable u2(M, N) nonnegative

            A_fix = s * T0 * sum(f_s2);
            B_fix = 0;
            for j = 1:N
                a_j = alpha_f(j);
                B_fix = B_fix + a_j * T0 * sum(w2(:,j)) ...
                    + a_j * s * T0 * n_flop / f_CPU_fix * sum(u2(:,j)) ...
                    + a_j * (s * T0 / R) * sum(u2(:,j));
            end

            maximize(A_fix - theta_fix * B_fix)
            subject to
                for i = 1:M
                    sum(x2(i, :)) <= 1;
                end
                sum(sum(w2)) <= P_total;

                % 时延约束: f_CPU固定后变为f_s的线性约束
                for j = 1:N
                    T0 + s * T0 * n_flop * f_s2(j) / f_CPU_fix + (s * T0 * f_s2(j)) / R <= T_max;
                end

                % 功率约束: C_CPU*(f_CPU_fix/100)^3 为常数
                for j = 1:N
                    P_static + C_S * f_s2(j) + C_CPU * pow_pos(f_CPU/100, 3) + P_tx(j) <= alpha_f(j) * sum(w2(:, j)) * 0.3;
                end

                f_s_min <= f_s2 <= f_s_max;
                P_fiber_min <= P_fiber2 <= P_fiber_max;

                for i = 1:M
                    for j = 1:N
                        % Big-M for w = x * P_fiber
                        w2(i,j) >= P_fiber2(i) - P_fiber_max * (1 - x2(i,j));
                        w2(i,j) <= P_fiber2(i);
                        w2(i,j) <= P_fiber_max * x2(i,j);
                        % McCormick for u = w * f_s
                        u2(i,j) >= P_fiber_min * f_s2(j) + f_s_min * w2(i,j) - P_fiber_min * f_s_min;
                        u2(i,j) >= P_fiber_max * f_s2(j) + f_s_max * w2(i,j) - P_fiber_max * f_s_max;
                        u2(i,j) <= P_fiber_min * f_s2(j) + f_s_max * w2(i,j) - P_fiber_min * f_s_max;
                        u2(i,j) <= P_fiber_max * f_s2(j) + f_s_min * w2(i,j) - P_fiber_max * f_s_min;
                    end
                end
        cvx_end

        if strcmp(cvx_status, 'Infeasible') || strcmp(cvx_status, 'Failed')
            fprintf('  FixCPU: 问题无解\n');
            break;
        end

        if abs(A_fix - theta_fix * B_fix) < epsilon
            ee_fix = A_fix / max(B_fix, 1e-10);
            fprintf('  FixCPU: EE = %.4f (迭代 %d 次)\n', ee_fix, iter);
            fprintf('  f_CPU = %.0f (固定)\n', f_CPU_fix);
            fprintf('  f_s   = '); disp(f_s2');
            fprintf('  x 分配:\n'); disp(full(x2));
            fprintf('  P_fiber = '); disp(P_fiber2');
            break;
        end
        theta_fix = A_fix / max(B_fix, 1e-10);

        if iter == max_iter
            ee_fix = A_fix / max(B_fix, 1e-10);
            fprintf('  FixCPU: EE = %.4f (达到最大迭代)\n', ee_fix);
        end
    end
    EE_results(sc, 2) = ee_fix;

    %% ---------- 算法3: AvgPower (固定光纤功率=P_total/M, 优化f_CPU和f_s) ----------
    % 参照 fixed_Pfiber.m: 功率固定后 w=x*Pf 变为精确等式
    % f_CPU 和 f_s 仍为变量, z 和 u 用已知 Pf 做 Big-M
    P_fiber_raw = (P_total / M) * ones(M, 1);
    P_fiber_capped = min(P_fiber_raw, P_fiber_max);
    f_min_avg = f_s_min / f_CPU_max;
    f_max_avg = f_s_max / f_CPU_min;
    ee_avg = 0;

    theta_avg = 0;
    for iter = 1:max_iter
        cvx_begin quiet
            cvx_solver mosek
            variable x3(M, N) binary
            variable f_CPU3(N) nonnegative
            variable f_s3(N) nonnegative
            variable f3(N) nonnegative
            variable w3(M, N) nonnegative
            variable z3(M, N) nonnegative
            variable u3(M, N) nonnegative

            A_avg = s * T0 * sum(f_s3);
            B_avg = 0;
            for j = 1:N
                a_j = alpha_f(j);
                B_avg = B_avg + a_j * T0 * sum(w3(:,j)) ...
                    + a_j * s * T0 * n_flop * sum(z3(:,j)) ...
                    + a_j * (s * T0 / R) * sum(u3(:,j));
            end

            maximize(A_avg - theta_avg * B_avg)
            subject to
                for i = 1:M
                    sum(x3(i, :)) <= 1;
                end

                % 时延约束 (与EE-Opt相同)
                for j = 1:N
                    T0 + s * T0 * n_flop * f3(j) + (s * T0 * f_s3(j)) / R <= T_max;
                end

                % 功率约束 (与EE-Opt相同)
                for j = 1:N
                    P_static + C_S * f_s3(j) + C_CPU * pow_pos(f_CPU3(j)/100, 3) + P_tx(j) <= alpha_f(j) * sum(w3(:, j)) * 0.3;
                end

                f_CPU_min <= f_CPU3 <= f_CPU_max;
                f_s_min <= f_s3 <= f_s_max;
                f_min_avg <= f3 <= f_max_avg;

                for j = 1:N
                    f3(j) <= f_s3(j) / f_CPU_min;
                    f3(j) >= (f_s_max/f_CPU_max) * (f_CPU_max - f_CPU3(j))/(f_CPU_max - f_CPU_min);
                end

                % 功率固定: w = x * Pf (精确等式), z/u 用已知 Pf 做 Big-M
                for i = 1:M
                    Pf = P_fiber_capped(i);
                    for j = 1:N
                        w3(i,j) == Pf * x3(i,j);

                        % Big-M for z = w * f
                        z3(i,j) <= Pf * f3(j);
                        z3(i,j) >= Pf * f3(j) - (1 - x3(i,j)) * Pf * f_max_avg;
                        z3(i,j) <= x3(i,j) * Pf * f_max_avg;
                        z3(i,j) >= x3(i,j) * Pf * f_min_avg;

                        % Big-M for u = w * f_s
                        u_U = Pf * f_s_max;
                        u_L = Pf * f_s_min;
                        u3(i,j) <= u_U * x3(i,j);
                        u3(i,j) >= u_L * x3(i,j);
                        u3(i,j) <= f_s3(j) * Pf;
                        u3(i,j) >= f_s3(j) * Pf - (1 - x3(i,j)) * u_U;
                    end
                end
        cvx_end

        if strcmp(cvx_status, 'Infeasible') || strcmp(cvx_status, 'Failed')
            fprintf('  AvgPower: 问题无解\n');
            break;
        end

        if abs(A_avg - theta_avg * B_avg) < epsilon
            ee_avg = A_avg / max(B_avg, 1e-10);
            fprintf('  AvgPower: EE = %.4f (迭代 %d 次)\n', ee_avg, iter);
            fprintf('  f_CPU = '); disp(f_CPU3');
            fprintf('  f_s   = '); disp(f_s3');
            fprintf('  x 分配:\n'); disp(full(x3));
            fprintf('  P_fiber (固定) = '); disp(P_fiber_capped');
            break;
        end
        theta_avg = A_avg / max(B_avg, 1e-10);

        if iter == max_iter
            ee_avg = A_avg / max(B_avg, 1e-10);
            fprintf('  AvgPower: EE = %.4f (达到最大迭代)\n', ee_avg);
        end
    end
    EE_results(sc, 3) = ee_avg;
end

%% ========== 汇总结果 ==========
fprintf('\n\n========================================\n');
fprintf('          非均匀部署实验结果汇总\n');
fprintf('========================================\n');
fprintf('%-12s  %10s  %10s  %10s\n', '场景', 'EE-Opt', 'FixCPU', 'AvgPower');
fprintf('----------------------------------------\n');
for sc = 1:3
    fprintf('%-12s  %10.4f  %10.4f  %10.4f\n', ...
        scenario_names{sc}, EE_results(sc,1), EE_results(sc,2), EE_results(sc,3));
end
fprintf('----------------------------------------\n');
fprintf('%-12s  %10.4f  %10.4f  %10.4f\n', '平均', ...
    mean(EE_results(:,1)), mean(EE_results(:,2)), mean(EE_results(:,3)));

% 计算相对提升
if mean(EE_results(:,3)) > 0
    imp_fix = (mean(EE_results(:,1)) - mean(EE_results(:,2))) / mean(EE_results(:,2)) * 100;
    imp_avg = (mean(EE_results(:,1)) - mean(EE_results(:,3))) / mean(EE_results(:,3)) * 100;
    fprintf('\nEE-Opt vs FixCPU:  +%.1f%%\n', imp_fix);
    fprintf('EE-Opt vs AvgPower: +%.1f%%\n', imp_avg);
end

%% ========== 可视化 ==========
figure('Name', '非均匀部署场景对比', 'Position', [200, 200, 700, 450]);

% 柱状图
bar_data = EE_results;
b = bar(bar_data);
b(1).FaceColor = [0.2 0.4 0.8];    % EE-Opt: 蓝色
b(2).FaceColor = [0.9 0.5 0.2];    % FixCPU: 橙色
b(3).FaceColor = [0.6 0.6 0.6];    % AvgPower: 灰色

set(gca, 'XTickLabel', scenario_names, 'FontSize', 12);
xlabel('部署场景');
ylabel('能效 (bit/J)');
title('非均匀部署场景下的能效对比 (N=5)');
legend('EE-Opt', 'FixCPU', 'AvgPower', 'Location', 'best');
grid on;

% 在柱顶标注数值
hold on;
for i = 1:3
    for j = 1:3
        if bar_data(i,j) > 0
            text(j + 0.27*(i-2), bar_data(i,j) + max(bar_data(:))*0.02, ...
                sprintf('%.2f', bar_data(i,j)), ...
                'HorizontalAlignment', 'center', 'FontSize', 9);
        end
    end
end
hold off;
