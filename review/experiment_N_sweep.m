%% N=3~8 参数扫描实验
% 基于 new_time_constraint_mosek_OE.m，循环运行不同设备数量
% 所有结果保存到 results_N_sweep.txt

clear all;
cvx_clear;
cvx_solver mosek;
cvx_precision low;
cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', 300);
cvx_solver_settings('MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-3);
cvx_expert true;

%% 公共参数
M = 10;
L0 = 5;
r = 1.5;
alpha_loss = 0.2;
P_total = 20;
T_max = 1.4;
oe_eta = 0.3;
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

%% 打开结果文件
fid = fopen('D:\OneDrive\文档\Post_student\code\paper_matlab\review\results_N_sweep.txt', 'w');
if fid == -1
    error('无法创建文件，请检查路径权限');
end
fprintf(fid, '============================================================\n');
fprintf(fid, '  EE-Opt Results: N = 3 to 8\n');
fprintf(fid, '  M=%d, P_total=%dW, P_fiber_max=%dW, eta_oe=%.1f\n', M, P_total, P_fiber_max, oe_eta);
fprintf(fid, '============================================================\n\n');

%% 汇总表格
fprintf(fid, '%-4s  %10s  %8s  %12s  %s\n', 'N', 'EE', 'Iter', 'Status', 'f_s');
fprintf(fid, '------------------------------------------------------------\n');

%% 主循环
for N = 3:8
    fprintf('\n===== N = %d =====\n', N);

    % 距离相关参数
    alpha_f = arrayfun(@(j) 10^(-alpha_loss * (L0 + r*(2*j-1)) / 10), 1:N);
    P_tx = arrayfun(@(j) 10^((5 + 0.2 * (L0 + r*(2*j-1))) / 10) / 1000, 1:N);

    % Dinkelbach 参数
    theta = 0;
    epsilon = 1e-4;
    max_iter = 50;
    iter = 0;
    ee_result = 0;
    cvx_status_str = '';

    while iter < max_iter
        iter = iter + 1;

        cvx_begin quiet
            cvx_solver mosek;

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
                end

                for j = 1:N
                    P_static + C_S * f_s(j) + C_CPU * pow_pos(f_CPU(j)/100, 3) + P_tx(j) <= oe_eta * alpha_f(j) * sum(w(:, j));
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

        cvx_status_str = cvx_status;

        if strcmp(cvx_status, 'Infeasible') || strcmp(cvx_status, 'Failed')
            fprintf('  N=%d: Infeasible/Failed at iter %d\n', N, iter);
            break;
        end

        if abs(A_x - theta * B_x) < epsilon || abs(B_x) < 1e-6
            ee_result = A_x / max(B_x, 1e-6);
            break;
        end
        theta = A_x / max(B_x, 1e-6);

        if iter == max_iter
            ee_result = A_x / max(B_x, 1e-6);
        end
    end

    %% 写入详细结果
    fprintf(fid, '\n===== N = %d =====\n', N);
    fprintf(fid, 'Status: %s\n', cvx_status_str);
    fprintf(fid, 'Dinkelbach iterations: %d\n', iter);
    fprintf(fid, 'EE (theta): %.6f\n', ee_result);

    if ~strcmp(cvx_status_str, 'Infeasible') && ~strcmp(cvx_status_str, 'Failed')
        fprintf(fid, '\nf_s   = [');
        for j = 1:N
            fprintf(fid, '%.4f ', f_s(j));
        end
        fprintf(fid, ']\n');

        fprintf(fid, 'f_CPU = [');
        for j = 1:N
            fprintf(fid, '%.2f ', f_CPU(j));
        end
        fprintf(fid, ']\n');

        fprintf(fid, 'P_fiber = [');
        for i = 1:M
            fprintf(fid, '%.4f ', P_fiber(i));
        end
        fprintf(fid, ']\n');

        fprintf(fid, 'x allocation:\n');
        x_full = full(x);
        for i = 1:M
            fprintf(fid, '  Fiber %2d -> ', i);
            assigned = find(x_full(i,:) == 1);
            if isempty(assigned)
                fprintf(fid, '(none)');
            else
                fprintf(fid, 'Device %d', assigned);
            end
            fprintf(fid, '\n');
        end

        % 约束检查
        fprintf(fid, '\nConstraint check:\n');
        for j = 1:N
            Lj = L0 + r * (2*j - 1);
            power_received = sum(full(w(:,j)) .* alpha_f(j)) * oe_eta;
            power_required = P_static + C_S * f_s(j) + C_CPU * (f_CPU(j)/100)^3 + P_tx(j);
            delay_j = T0 + f(j) * s * T0 * n_flop + (f_s(j) * s * T0) / R;
            fprintf(fid, '  Device %d (L=%.1fkm): P_recv=%.4fW, P_req=%.4fW [%s] | Delay=%.4fs [%s]\n', ...
                j, Lj, power_received, power_required, ...
                iif(power_received >= power_required - 1e-4, 'OK', 'VIOLATED'), ...
                delay_j, iif(delay_j <= T_max + 1e-4, 'OK', 'VIOLATED'));
        end

        total_power = sum(sum(full(w)));
        fprintf(fid, '  Total power: %.4fW / %dW [%s]\n', total_power, P_total, ...
            iif(total_power <= P_total + 1e-4, 'OK', 'VIOLATED'));
    end

    fprintf(fid, '------------------------------------------------------------\n');

    % 写入汇总表
    f_s_str = sprintf('[%.2f ', f_s(1));
    for j = 2:N
        f_s_str = [f_s_str, sprintf('%.2f ', f_s(j))];
    end
    f_s_str = [f_s_str, ']'];
    fprintf(fid, '%-4d  %10.6f  %8d  %12s  %s\n', N, ee_result, iter, cvx_status_str, f_s_str);

    fprintf('  N=%d: EE = %.6f (%s, %d iter)\n', N, ee_result, cvx_status_str, iter);
end

%% 关闭文件
fprintf(fid, '\n============================================================\n');
fprintf(fid, '  Done.\n');
fprintf(fid, '============================================================\n');
fclose(fid);

fprintf('\nAll results saved to results_N_sweep_18W.txt\n');


%% 辅助函数
function out = iif(cond, trueVal, falseVal)
    if cond
        out = trueVal;
    else
        out = falseVal;
    end
end
