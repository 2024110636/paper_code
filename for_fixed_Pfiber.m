%% 使用 MOSEK 求解器（CVX接口）
cvx_clear;
cvx_solver mosek;
cvx_precision low
cvx_solver_settings('MSK_DPAR_OPTIMIZER_MAX_TIME', 300);
cvx_solver_settings('MSK_DPAR_INTPNT_TOL_REL_GAP', 1e-3);
cvx_expert true;

%% 设置结果文件路径
target_dir = 'D:\OneDrive\文档\Post_student\model\论文';
result_filename = 'results1.txt';

% 转换为标准路径格式
target_dir = strrep(target_dir, '/', '\');

% 创建目录（如果不存在）
[status, msg] = mkdir(target_dir);  
if ~status
    error('目录创建失败: %s\n原因: %s', target_dir, msg);
end

% 检查目录可写性
test_file = fullfile(target_dir, 'write_test.tmp');
fid = fopen(test_file, 'w');
if fid == -1
    error('目录不可写: %s', target_dir);
else
    fclose(fid);
    delete(test_file);
end

% 打开结果文件
result_file_path = fullfile(target_dir, result_filename);
[result_file, errmsg] = fopen(result_file_path, 'w', 'n', 'UTF-8');
if result_file == -1
    error('文件创建失败: %s\n错误: %s', result_file_path, errmsg);
end

fprintf(result_file, '不同N值下的优化结果比较\n');
fprintf(result_file, '========================\n\n');

%% 参数设置（N将在循环中变化）
M = 10;        % 光纤数量
L0 = 5;        % 基础距离 (km)
r = 1.5;       % 每个设备增加距离 (km)
alpha_loss = 0.2;   % 光纤损耗 (dB/km)
P_total = 15;       % 总供电功率 (W)
T_max = 1.4;        % 最大允许时延 (s)

% 固定功率并裁剪
P_fiber_raw = (P_total / M) * ones(M, 1);
P_fiber_capped = min(P_fiber_raw, 2);  % 裁剪上限为 2W

% 设备功率模型
P_static = 0.1;       % 静态功率 (W)
C_S = 0.5;            % 采样功率系数 (W/MHz)
C_CPU = 0.001;        % CPU功率系数 (W/GHz^3)
s = 2;                % 每个样本大小 (bit)
T0 = 1;               % 采样时间 (s)
n_flop = 100;         % 每个样本浮点操作数
R = 500;              % 光纤速率 (Mbit/s)

% CPU & 采样频率 (单位：Hz)
f_CPU_min = 500; f_CPU_max = 1000;   % MHz
f_s_min = 1;      f_s_max = 2;        % MHz

%% 循环遍历不同的N值
for N = 3:8
    fprintf(result_file, '\n===== N = %d =====\n', N);
    fprintf('正在计算 N = %d...\n', N);
    
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
            variable f_CPU(N) nonnegative
            variable f_s(N) nonnegative
            
            % 辅助变量
            variable t(N) nonnegative
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
                % 每根光纤最多连接一个设备
                for i = 1:M
                    sum(x(i, :)) <= 1;
                end
                
                % 时延约束
                for j = 1:N
                    T0 + s * T0 * n_flop * f(j) + (s * T0 * f_s(j)) / R <= T_max;
                end
                
                % 设备功率约束
                for j = 1:N
                    P_static + C_S * f_s(j) + C_CPU * pow_pos(f_CPU(j)/100, 3) + P_tx(j) <= alpha_f(j) * sum(w(:, j));
                end
                
                % 变量范围
                f_CPU_min <= f_CPU <= f_CPU_max;
                f_s_min <= f_s <= f_s_max;
                f_min = f_s_min / f_CPU_max;
                f_max = f_s_max / f_CPU_min;
                f_min <= f <= f_max;
                
                for j = 1:N
                    f(j) <= f_s(j)/(f_CPU_min);
                    f(j) >= (f_s_max/f_CPU_max) * (f_CPU_max - f_CPU(j))/(f_CPU_max - f_CPU_min);
                end
                
                % 线性化 w = x * capped_power，z = w * f，u = w * f_s
                for i = 1:M
                    Pf = P_fiber_capped(i);
                    for j = 1:N
                        % w(i,j) = x(i,j) * P_fiber_capped(i)
                        w(i,j) == Pf * x(i,j);
                        
                        %% Big-M for z(i,j) = w(i,j) * f(j)
                        z(i,j) <= Pf * f(j);
                        z(i,j) >= Pf * f(j) - (1 - x(i,j)) * Pf * f_max;
                        z(i,j) <= x(i,j) * Pf * f_max;
                        z(i,j) >= x(i,j) * Pf * f_min;
                        
                        %% Big-M for u(i,j) = w(i,j) * f_s(j)
                        u_U = Pf * f_s_max;
                        u_L = Pf * f_s_min;
                        u(i,j) <= u_U * x(i,j);
                        u(i,j) >= u_L * x(i,j);
                        u(i,j) <= f_s(j) * Pf;
                        u(i,j) >= f_s(j) * Pf - (1 - x(i,j)) * u_U;
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
    
    %% 将结果写入文件
    fprintf(result_file, 'Dinkelbach 迭代收敛:\n');  
    fprintf(result_file, '最终 theta: %.6f\n', theta);  
    fprintf(result_file, '迭代次数: %d\n', iter);  
    
    fprintf(result_file, '光纤分配矩阵 x:\n');  
    % 将稀疏矩阵转换为完整矩阵
    x_full = full(x);
    for i = 1:size(x_full,1)
        fprintf(result_file, '%d ', x_full(i,:));
        fprintf(result_file, '\n');
    end
    
    fprintf(result_file, '光纤功率 P_fiber: %.4f \n', P_fiber_capped);
    fprintf(result_file, '\n CPU 频率 f_CPU (MHz):%.2f \n ', f_CPU);
    fprintf(result_file, '\n 采样频率 f_s (MHz): %.2f \n', f_s);
    
    fprintf(result_file, '\n 最优目标值: %.6f \n', cvx_optval);  
    fprintf(result_file, '分子 A_x: %.6f \n 分母 B_x: %.6f \n', A_x, B_x);
    
    %% 频率估算误差分析
    %% 修改后的结果输出部分（仅修正fprintf错误）

%% 频率估算误差分析
f_CPU_approx = f_s ./ f;
relative_error = (f_CPU - f_CPU_approx) ./ f_CPU;

fprintf(result_file, '\nCPU频率估算误差:\n真实值\t估算值\t相对误差\n');
for j = 1:N
    fprintf(result_file, '%.2f \t %.2f \t %.4f \n', ...
        full(f_CPU(j)), full(f_CPU_approx(j)), full(relative_error(j)));
end

%% 检查每个设备的时延约束
fprintf(result_file, '\n时延约束检查 (T_max=%.1fs):\n', T_max);
for j = 1:N
    delay = T0 + full(f(j)) * s * T0 * n_flop + (full(f_s(j)) * s * T0)/R;
    status = ternary(delay <= T_max, '满足', '不满足');
    fprintf(result_file, '设备 %d: %.4fs --> %s\n', j, full(delay), status);
end

%% 检查每个设备的功率约束
fprintf(result_file, '\n功率约束检查:\n');
for j = 1:N
    power_lhs = P_static + C_S * full(f_s(j)) + C_CPU * (full(f_CPU(j))/100)^3 + full(P_tx(j));
    power_rhs = sum(full(w(:, j))) * full(alpha_f(j));
    status = ternary(power_lhs <= power_rhs, '满足', '不满足');
    fprintf(result_file, '设备 %d: LHS=%.4fW, RHS=%.4fW --> %s\n', ...
        j, full(power_lhs), full(power_rhs), status);
end
    
    fprintf(result_file, '\n%s\n', repmat('-', 1, 60));
end

%% 关闭文件和完成提示
fclose(result_file);
fprintf('计算完成！结果已保存至:\n%s\n', result_file_path);

%% 辅助函数
function out = ternary(cond, trueStr, falseStr)
    if cond
        out = trueStr;
    else
        out = falseStr;
    end
end