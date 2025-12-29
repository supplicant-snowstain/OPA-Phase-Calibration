function opa_ga_peakfocus_pointing_live_imaq_2D_AutoStop()
% ============================================================
% OPA 初始相位校准 —— GA（遗传算法）2D版 (智能截止优化)
% 
% 全图滑动匹配 (Sliding Window / Convolution)
%              不依赖亮度峰值，寻找形状最佳匹配点。
% ============================================================
try
    oldObj = instrfindall;
    if ~isempty(oldObj)
        fclose(oldObj);
        delete(oldObj);
        clear oldObj;
    end
catch
end
%% ============ 参数区 ============
% ---------- 硬件与通道 ----------
COMPORT     = 'COM3';
CH_LIST     = 1:16;
% ---------- 电压范围 ----------
VMIN        = 0;
VMAX        = 5.5; 
V_STEP      = 0.001;        
% ---------- 写入 ramp ----------
RAMP_DV     = 0.05;
RAMP_DT     = 0.08;
SETTLE_T    = 0.20;
% ---------- 2D ROI (搜索大区域) ----------
RECT = [0, 580, 1280, 80]; 
LOCK_RECT_FROM_FIRST = true;

% ---------- 2D 模板参数 (滑动窗口大小) ----------
HY_MODE      = 'gaussian';   
% 这里定义了模版的大小 (例如 200x50)，也是滑动窗口的大小
SUBRECT_SIZE = [200, 50]; % [Width, Height]
HY_SIGMA_PX  = 5;         

% ---------- 预处理 ----------
PP_REMOVE_BASE  = true;
PP_L2_NORM      = true; % 注意：滑动模式下，归一化在内部计算

% ---------- GA 基础参数 ----------
GA_POP        = 50;       
GA_ELITE      = 2;
GA_CROSS_FRAC = 0.80;
GA_MUT_FCN    = @mutationadaptfeasible;
GA_AVG_FRAMES = 1;

% ---------- 初始电压 ----------
V_INIT = zeros(1, numel(CH_LIST));
% ---------- 物理参数 ----------
lambda        = 1550e-9;
d             = 775e-9;
V_pi_val      = 3.08; 
V_2pi_val     = V_pi_val * sqrt(2);

% ---------- 指向模块（保留可选） ----------
DO_POINTING   = false;
THETA_DEG     = 10;
lambda        = 1550e-9;
d             = 775e-9;
P_pi_mW       = 21.9135;
R_ohm         = 434.78;
WRAP_PHASE    = true;
USE_PUSH_PULL = true;
% ---------- 日志 ----------
SAVE_DIR = fullfile(pwd,'ga_2D_Sliding');
if ~exist(SAVE_DIR,'dir'), mkdir(SAVE_DIR); end
LOG_CSV  = fullfile(SAVE_DIR,'progress.csv');

%% ============ 打开串口 & 初值写入 ============
FPGA = serial(COMPORT);
set(FPGA,'BaudRate',115200,'StopBits',1,'DataBits',8);
fopen(FPGA);
fprintf('> Port %s opened.\n', COMPORT);
N = numel(CH_LIST);
V = clamp_vec(V_INIT(:), VMIN, VMAX);
V = quantize_voltage(V, VMIN, V_STEP);
write_vec_soft(FPGA, CH_LIST, V, VMIN, VMAX, RAMP_DV, RAMP_DT);
fprintf('> All channels set to 0.0 V\n');

%% ============ 相机初始化 ============
g          = gigecam;
if isprop(g,'Timeout'), g.Timeout = 2.0; end
hImg       = [];
fail_count = 0;

%% ============ 首帧采图 & 构建模版 ============
[I0, g, hImg, fail_count] = capture_frame_live_safe_minimal(g, hImg, fail_count);
I0 = double_gray(I0);
[H0,W0] = size(I0);
rect = crop_rect_to_image(RECT, W0, H0);
if LOCK_RECT_FROM_FIRST
    rect_locked = rect;
else
    rect_locked = [];
end

% 显示首帧搜索区
show_first_frame_with_roi(I0, rect);

% ====== 【滑动模式修改】构造 2D 参考模板 H2D ======
% 我们直接生成一个指定大小 (200x50) 的模版矩阵
Lx = SUBRECT_SIZE(1); % 200
Ly = SUBRECT_SIZE(2); % 50
[H2D, Hinfo] = build_H2D_pure_gaussian(Lx, Ly, HY_SIGMA_PX);
fprintf('[H2D] Template built: %dx%d, Sigma=%d. Ready for sliding.\n', Lx, Ly, HY_SIGMA_PX);

% 可视化模版
figure('Name','Template','Color','w');
imagesc(H2D); axis image; title('Sliding Template'); colorbar;

% ====== 收敛曲线初始化 ======
[figC, axC, hCost, hK] = init_convergence_plot_K();
write_csv_header_if_needed(LOG_CSV, 'eval,cost,K');

% 初始化实时截面监控窗口
[hProfLines] = init_profile_monitor();  
assignin('base', 'hProfLines', hProfLines); % 存入工作区供函数调用

write_csv_header_if_needed(LOG_CSV, 'eval,cost,K');

eval_step = 0;

% 初始化运行时结构体
RT = struct();
RT.FPGA = FPGA; RT.g = g; RT.hImg = hImg; RT.fail_count = fail_count;
RT.CH_LIST = CH_LIST;
RT.VMIN = VMIN; RT.VMAX = VMAX; RT.V_STEP = V_STEP;
RT.RAMP_DV = RAMP_DV; RT.RAMP_DT = RAMP_DT; RT.SETTLE_T =SETTLE_T;

% 【滑动模式修改】初始评价 - 传入 H2D 矩阵而非一维 refProfile
[c0, K0] = evaluate_once_sliding(RT, V, H2D, ...
    rect, rect_locked, LOCK_RECT_FROM_FIRST, ...
    PP_REMOVE_BASE, ...
    GA_AVG_FRAMES);

eval_step = eval_step + 1;
addpoints(hCost, eval_step, c0);
addpoints(hK,    eval_step, K0);
title(axC, sprintf('Init: cost=%.4f | K=%.4f', c0, K0), 'FontWeight','normal');
drawnow limitrate;
append_csv(LOG_CSV, eval_step, c0, K0);

%% ============ GA 初始种群与参数配置 ============
fprintf('\n==== Entering GA optimization (Sliding Mode) ====\n');
lb = VMIN * ones(1, N);
ub = VMAX * ones(1, N);
% 【滑动模式修改】传入 H2D
fitnessFcn = @(Vrow) ga_cost(RT, Vrow, H2D, rect, rect_locked, LOCK_RECT_FROM_FIRST, PP_REMOVE_BASE, GA_AVG_FRAMES); 

% 生成物理种子
num_seeds = 20;
PopInit = zeros(num_seeds, N);
max_slope = V_pi_val;
possible_slopes = linspace(-max_slope, max_slope, num_seeds);
for i = 1:num_seeds
    d_volts = possible_slopes(i);
    vec = (1:N) * d_volts;
    vec = mod(vec, V_2pi_val);
    PopInit(i, :) = vec;
end

% GA 选项
target_K = 0.88; 
opts = optimoptions('ga', ...
    'PopulationSize', GA_POP, ...
    'InitialPopulationMatrix', PopInit, ...
    'MaxGenerations', 50, ...
    'MaxStallGenerations', 8, ...
    'FunctionTolerance', 1e-4, ...
    'FitnessLimit', 1 - target_K, ...
    'EliteCount', GA_ELITE, ...
    'CrossoverFraction', GA_CROSS_FRAC, ...
    'MutationFcn', GA_MUT_FCN, ...
    'UseParallel', false, ...
    'Display', 'iter', ...
    'PlotFcn', {@gaplotbestf});

% 运行 GA
[Vbest_row, bestCost] = ga(fitnessFcn, N, [],[],[],[], lb, ub, [], opts);

% ... 后续处理 ...
Vbest = Vbest_row(:);
Vbest = clamp_vec(Vbest, VMIN, VMAX);
Vbest = quantize_voltage(Vbest, VMIN, V_STEP);
fprintf('\n[GA] Done. bestCost=%.6f (Final K=%.6f)\n', bestCost, 1-bestCost);

% 最终结果验证
% 1. 使用主程序中现有的 RT 结构体写入最佳电压
write_vec_soft(RT.FPGA, RT.CH_LIST, Vbest, RT.VMIN, RT.VMAX, RT.RAMP_DV, RT.RAMP_DT);

% 2. 等待热稳定
pause(RT.SETTLE_T);

% 3. 采最终帧 (直接使用 RT 中的句柄)
[I_final, RT.g, RT.hImg, RT.fail_count] = capture_frame_live_safe_minimal(RT.g, RT.hImg, RT.fail_count);

% 4. 图像处理
I_final = double_gray(I_final);
if LOCK_RECT_FROM_FIRST
    rect_eval = rect_locked;
else
    rect_eval = crop_rect_to_image(RECT, size(I_final,2), size(I_final,1));
end

% 5. 【滑动模式】最终计算 K 值
% 注意：这里直接调用 sliding_window_max_K
[K_final, best_patch, best_loc] = sliding_window_max_K( ...
    I_final(rect_eval(2):rect_eval(2)+rect_eval(4)-1, rect_eval(1):rect_eval(1)+rect_eval(3)-1), ...
    H2D, ...
    PP_REMOVE_BASE);

% 6. 可视化对比
figure('Name','Final Sliding Result','Color','w');
subplot(1,3,1); imagesc(H2D); axis image; title('Template');
subplot(1,3,2); imagesc(best_patch); axis image; title(sprintf('Best Match\nK=%.4f', K_final));
subplot(1,3,3); imagesc(I_final); axis image; hold on;
% 画出最终匹配到的位置
rx = rect_eval(1) + best_loc(1) - 1;
ry = rect_eval(2) + best_loc(2) - 1;
rectangle('Position',[rx, ry, size(H2D,2), size(H2D,1)], 'EdgeColor','g','LineWidth',2);
title('Final Frame Location');

% 7. 保存结果
save(fullfile(SAVE_DIR,'final_ga_2D_Sliding.mat'), 'Vbest','bestCost','K_final','RT');
writematrix(Vbest, fullfile(SAVE_DIR,'final_voltage_ga.csv'));
% 【滑动模式修改】最终计算
% [K_final, best_patch, best_loc] = sliding_window_max_K(I_final(rect_eval(2):rect_eval(2)+rect_eval(4)-1, rect_eval(1):rect_eval(1)+rect_eval(3)-1), H2D, PP_REMOVE_BASE);

%% ============ 指向模块（保留可选） ============
if DO_POINTING
    [V_set, dphi_adj, dphi_n, dP_mW] = apply_theta_pointing( ...
        Vbest, THETA_DEG, lambda, d, P_pi_mW, R_ohm, WRAP_PHASE, USE_PUSH_PULL, ...
        VMIN, VMAX, RAMP_DV, RAMP_DT, FPGA, CH_LIST);

    fprintf('> Pointing applied: theta=%.3f deg | dphi_adj=%.6f rad | mean(V)=%.4f V\n', ...
            THETA_DEG, dphi_adj, mean(V_set));
    writematrix(V_set, fullfile(SAVE_DIR, sprintf('voltage_pointing_%+ddeg.csv', round(THETA_DEG))));
    save(fullfile(SAVE_DIR, sprintf('pointing_%+ddeg.mat', round(THETA_DEG))), ...
         'THETA_DEG','lambda','d','P_pi_mW','R_ohm','WRAP_PHASE','Vbest','V_set','dphi_adj','dphi_n','dP_mW');
end
%% ============ 收尾 ============
try clear g; end %#ok<TRYNC>
try fclose(FPGA); end %#ok<TRYNC>
try delete(instrfindall); end %#ok<TRYNC>
fprintf('> Port closed.\n');


%% ============================================================
%                 GA 代价函数
% ============================================================
function cost = ga_cost(RT, Vrow, H2D, rect, rect_locked, is_locked, pp_remove, avg_frames) 
    Vcand = clamp_vec(Vrow(:), RT.VMIN, RT.VMAX);
    Vcand = quantize_voltage(Vcand, RT.VMIN, RT.V_STEP);
    
    [cost, K] = evaluate_once_sliding(RT, Vcand, H2D, ...
        rect, rect_locked, is_locked, ...
        pp_remove, ...
        avg_frames);
    
    % 更新图表 (需访问主函数变量，这里简化处理，实际运行时如果报错可去掉绘图部分或用 globals)
    % 简单起见，这里假设 handle 可用，如果不可用请注释掉绘图更新
    try
        hCost = evalin('base','hCost');
        hK    = evalin('base','hK');
        axC   = evalin('base','axC');
        LOG_CSV = evalin('base','LOG_CSV');
        eval_step = evalin('base','eval_step') + 1;
        assignin('base','eval_step', eval_step);
        
        addpoints(hCost, eval_step, cost);
        addpoints(hK,    eval_step, K);
        title(axC, sprintf('eval=%d | cost=%.4f | K=%.4f', eval_step, cost, K), 'FontWeight','normal');
        drawnow limitrate;
        append_csv(LOG_CSV, eval_step, cost, K);
    catch
    end
    
    if ~isfinite(cost), cost = 1.0; end
end

end 

%% ============================================================
%           【核心修改】单次评价：滑动窗口版 (含实时绘图)
% ============================================================
function [cost, K, RT] = evaluate_once_sliding( ...
    RT, Vvec, H2D, rect, rect_locked, lock_rect, ...
    pp_remove_base, avg_frames)

    write_vec_soft(RT.FPGA, RT.CH_LIST, Vvec, RT.VMIN, RT.VMAX, RT.RAMP_DV, RT.RAMP_DT);
    pause(RT.SETTLE_T);
    
    K_acc = 0;
    
    % 准备变量用于最后一次绘图
    last_K_Map = [];
    last_best_patch = [];
    
    for t = 1:avg_frames
        [I, RT.g, RT.hImg, RT.fail_count] = capture_frame_live_safe_minimal(RT.g, RT.hImg, RT.fail_count);
        I = double(I);
        
        % 1. 裁剪搜索区域
        if lock_rect
            rect_eval = rect_locked;
        else
            rect_eval = crop_rect_to_image(rect, size(I,2), size(I,1));
        end
        x1 = rect_eval(1); y1 = rect_eval(2); w = rect_eval(3); h = rect_eval(4);
        I_ROI = I(y1:y1+h-1, x1:x1+w-1);
        
        % 2. 执行滑动计算 (注意：这里多接收一个 map 输出)
        [K_best, best_patch, ~, K_Map] = sliding_window_max_K(I_ROI, H2D, pp_remove_base);
        
        K_acc = K_acc + K_best;
        
        % 记录最后一次的数据用于显示
        if t == avg_frames
            last_K_Map = K_Map;
            last_best_patch = best_patch;
        end
    end
    
    K = K_acc / avg_frames;
    cost = 1 - K;
    
    % ====== 【新增】实时更新截面曲线 ======
    try
        h = evalin('base', 'hProfLines'); % 从主程序获取句柄
        if isvalid(h.fig)
            % 1. 更新滑动搜索曲线 (把2D K_Map 沿Y轴取最大值，变成1D曲线)
            % 这代表：在X轴的每一个位置，纵向上能找到的最好K值是多少
            scan_curve = max(last_K_Map, [], 1); 
            set(h.lineScan, 'XData', 1:numel(scan_curve), 'YData', scan_curve);
            
            % 2. 准备截面数据
            % 实际光斑的中心切片
            mid_y = round(size(last_best_patch,1)/2);
            mid_x = round(size(last_best_patch,2)/2);
            profX_Real = last_best_patch(mid_y, :);
            profY_Real = last_best_patch(:, mid_x);
            
            % 模版的中心切片 (归一化到和实际光斑同高度，方便对比)
            profX_Ref = H2D(mid_y, :);
            profY_Ref = H2D(:, mid_x);
            scale_factor = max(last_best_patch(:)); % 简单缩放以便重叠显示
            
            % 3. 更新截面图
            set(h.lineX_Real, 'XData', 1:numel(profX_Real), 'YData', profX_Real);
            set(h.lineX_Ref,  'XData', 1:numel(profX_Ref),  'YData', profX_Ref * scale_factor);
            
            set(h.lineY_Real, 'XData', 1:numel(profY_Real), 'YData', profY_Real);
            set(h.lineY_Ref,  'XData', 1:numel(profY_Ref),  'YData', profY_Ref * scale_factor);
            
            drawnow limitrate;
        end
    catch
        
    end
end
% % %% ============================================================
% % %           【核心修改】单次评价：滑动窗口版
% % % ============================================================
% % function [cost, K, RT] = evaluate_once_sliding( ...
% %     RT, Vvec, H2D, rect, rect_locked, lock_rect, ...
% %     pp_remove_base, avg_frames)
% % 
% %     
% %     write_vec_soft(RT.FPGA, RT.CH_LIST, Vvec, RT.VMIN, RT.VMAX, RT.RAMP_DV, RT.RAMP_DT);
% %     pause(RT.SETTLE_T);
% %     
% %     K_acc = 0;
% %     for t = 1:avg_frames
% %         [I, RT.g, RT.hImg, RT.fail_count] = capture_frame_live_safe_minimal(RT.g, RT.hImg, RT.fail_count);
% %         I = double(I);
% %         
% %         % 1. 裁剪出搜索大区域 (1280x100)
% %         if lock_rect
% %             rect_eval = rect_locked;
% %         else
% %             rect_eval = crop_rect_to_image(rect, size(I,2), size(I,1));
% %         end
% %         x1 = rect_eval(1); y1 = rect_eval(2); w = rect_eval(3); h = rect_eval(4);
% %         I_ROI = I(y1:y1+h-1, x1:x1+w-1);
% %         
% %         % 2. 【核心】执行全图滑动计算 K 值，返回最大的那个 K
% %         [K_best, ~] = sliding_window_max_K(I_ROI, H2D, pp_remove_base);
% %         
% %         K_acc = K_acc + K_best;
% %     end
% %     
% %     K = K_acc / avg_frames;
% %     cost = 1 - K;
% % end

%% ============================================================
%           【核心算法】滑动窗口计算最大 K 值
%   原理：利用卷积 (conv2) 快速实现“模版在图像上滑动挪1000次”
% ============================================================
function [max_K, best_patch, best_loc, K_Map] = sliding_window_max_K(Img, Template, remove_base)
    % Img: 搜索区域图像 (例如 100x1280)
    % Template: 模版 (例如 50x200)
    
    % 1. 预处理
    if remove_base
        Img = Img - min(Img(:));
        Img(Img<0) = 0;
    end
    
    % 2. 准备分母的一部分：模版能量 (常数)
    norm_T = norm(Template(:));
    
    % 3. 快速滑动计算
    % ---------------------------------------------------------
    % 我们要算的是： Σ(G*H) / (√ΣG² * √ΣH²)
    
    % A. 分子：Σ(G*H) —— 这就是卷积！
    % 使用 'valid' 模式，只计算模版完全落在图像内的区域，相当于不填充边缘
    % Template 也是高斯对称的，所以卷积=互相关，不用旋转180度
    Numerator_Map = conv2(Img, Template, 'valid');
    
    % B. 分母：√ΣG² —— 这是局部能量图
    % 技巧：计算 Img.^2 和 全1矩阵 的卷积，就能得到每个窗口内的平方和
    Ones_Mask = ones(size(Template));
    Local_Energy_Sq_Map = conv2(Img.^2, Ones_Mask, 'valid');
    Denominator_G_Map = sqrt(Local_Energy_Sq_Map);
    
    % C. 计算 K 值地图 (K Map)
    % 此时 K_Map 的每一个点，代表模版滑动到该位置时的 K 值
    K_Map = Numerator_Map ./ (Denominator_G_Map * norm_T + 1e-12);
    
    % ---------------------------------------------------------
    
    % 4. 寻找最大值 (The Best Match)
    [max_K, max_idx] = max(K_Map(:));
    
    % 5. (可选) 提取最佳匹配的位置和图像块，用于显示
    [y_peak, x_peak] = ind2sub(size(K_Map), max_idx);
    
    % 注意：conv2 'valid' 的结果坐标需要换算回原图坐标
    % 卷积结果的 (1,1) 对应模版左上角在原图 (1,1) 的情况
    % 所以最佳匹配块的左上角就是 (x_peak, y_peak)
    best_loc = [x_peak, y_peak]; 
    
    H_h = size(Template, 1);
    H_w = size(Template, 2);
    
    best_patch = Img(y_peak : y_peak + H_h - 1, ...
                     x_peak : x_peak + H_w - 1);
end


%% ============================================================
%           【新增】初始化实时截面监控窗口
% ============================================================
function h = init_profile_monitor()
    fig = figure('Name','Real-time Profile Monitor','Color','w');
    
    % 1. 上图：滑动搜索结果 (K值分布)
    ax1 = subplot(2,2,[1,2], 'Parent', fig); hold(ax1,'on'); grid(ax1,'on');
    h.lineScan = plot(ax1, NaN, NaN, 'b-', 'LineWidth', 1.5);
    title(ax1, 'Sliding Search Result (Max K along X-axis)');
    ylabel(ax1, 'K Value'); xlabel(ax1, 'X Position in RECT');
    ylim(ax1, [0, 1]);
    
    % 2. 左下图：X轴截面对比
    ax2 = subplot(2,2,3, 'Parent', fig); hold(ax2,'on'); grid(ax2,'on');
    h.lineX_Real = plot(ax2, NaN, NaN, 'b-', 'LineWidth', 1.5, 'DisplayName','Real');
    h.lineX_Ref  = plot(ax2, NaN, NaN, 'r--', 'LineWidth', 1.5, 'DisplayName','Template');
    title(ax2, 'X-Axis Profile (Best Match)');
    legend(ax2, 'Location','best');
    
    % 3. 右下图：Y轴截面对比
    ax3 = subplot(2,2,4, 'Parent', fig); hold(ax3,'on'); grid(ax3,'on');
    h.lineY_Real = plot(ax3, NaN, NaN, 'b-', 'LineWidth', 1.5, 'DisplayName','Real');
    h.lineY_Ref  = plot(ax3, NaN, NaN, 'r--', 'LineWidth', 1.5, 'DisplayName','Template');
    title(ax3, 'Y-Axis Profile (Best Match)');
    
    h.fig = fig;
end

%% ============================================================
%                 指向模块（保留）
% ============================================================
function [V_set, dphi_adj, dphi_n, dP_mW] = apply_theta_pointing( ...
        V_base, THETA_DEG, lambda, d, P_pi_mW, R_ohm, WRAP_PHASE, USE_PUSH_PULL, ...
        VMIN, VMAX, RAMP_DV, RAMP_DT, FPGA, CH_LIST)

N = numel(CH_LIST);
theta_rad  = deg2rad(THETA_DEG);
dphi_adj   = (2*pi/lambda) * d * sin(theta_rad);
c          = (N+1)/2;
m_idx      = (1:N) - c;
dphi_n     = m_idx(:) .* dphi_adj;
if WRAP_PHASE
    dphi_n = atan2(sin(dphi_n), cos(dphi_n));
end

dP_mW      = (P_pi_mW/pi) .* dphi_n;

P0_W       = (V_base(:).^2) ./ R_ohm;
P_target_W = P0_W + dP_mW(:) * 1e-3;
P_target_W = max(0, P_target_W);
V_set      = sqrt(R_ohm .* P_target_W);

if USE_PUSH_PULL
    V_set = V_set - mean(V_set) + mean(V_base(:));
end

V_set = clamp_vec(V_set, VMIN, VMAX);
for i = 1:N
    write_vec_soft(FPGA, CH_LIST(i), V_set(i), VMIN, VMAX, RAMP_DV, RAMP_DT);
end
end


%% ============================================================
%                 辅助函数 
% ============================================================
function [H2D, info] = build_H2D_pure_gaussian(Lx, Ly, sigma_px)
    cx = (Lx + 1) / 2;
    cy = (Ly + 1) / 2;
    [xx, yy] = meshgrid(1:Lx, 1:Ly);
    dist_sq = (xx - cx).^2 + (yy - cy).^2;
    H2D = exp(- dist_sq / (2 * sigma_px^2));
    H2D = H2D / max(H2D(:)); % 峰值归一化为1
    info.Lx = Lx; info.Ly = Ly; info.sigma = sigma_px;
end

function [I, g, hImg, fail_count] = capture_frame_live_safe_minimal(g, hImg, fail_count)
    while true
        try
            if isempty(g) || ~isvalid(g)
                g = gigecam; if isprop(g,'Timeout'), g.Timeout=2.0; end
            end
            Iraw = snapshot(g);
            if isempty(hImg) || ~isvalid(hImg)
                figure('Name','Live','Color','k','NumberTitle','off');
                hImg = imshow(double_gray(Iraw), [], 'Border','tight'); title('Live','Color','w'); drawnow;
            else
                set(hImg, 'CData', Iraw); drawnow limitrate;
            end
            I = Iraw; return;
        catch
            fail_count = fail_count + 1; 
            try clear g; 
            end; g=[]; pause(0.2);
        end
    end
end
function D = double_gray(I)
if ndims(I) == 3, I = rgb2gray(I); end
if isa(I,'double'), D = I; else, D = double(I); end
end
function rect = crop_rect_to_image(RECT, W, H)
x1 = max(1, round(RECT(1))); y1 = max(1, round(RECT(2)));
x2 = min(W, x1 + round(RECT(3)) - 1); y2 = min(H, y1 + round(RECT(4)) - 1);
w = max(0, x2-x1+1); h = max(0, y2-y1+1);
rect = [x1, y1, w, h];
end
function V = clamp_vec(V, vmin, vmax), V = min(vmax, max(vmin, V)); end
function Vq = quantize_voltage(V, vmin, vstep)
if vstep <= 0, Vq = double(V(:)); return; end
Vq = vmin + round((double(V(:)) - vmin)/vstep) * vstep;
end
function write_vec_soft(FPGA, CH_LIST, V_target, vmin, vmax, dv, dt)
persistent Vcache
if isempty(Vcache), Vcache = zeros(1, 128); end


V_target = clamp_vec(V_target(:), vmin, vmax);
for i = 1:numel(CH_LIST)
    ch = CH_LIST(i); vt = V_target(i);
    VC96_setVone(FPGA, ch, vt);
    Vcache(ch) = vt;
end
end
function show_first_frame_with_roi(I, rect)
figure('Name','First frame','Color','w','NumberTitle','off');
imshow(I, [], 'Border','tight'); hold on;
rectangle('Position', rect, 'EdgeColor', 'r', 'LineWidth', 1.8);
hold off; drawnow;
end
function [figC, axC, hCost, hK] = init_convergence_plot_K()
figC = figure('Name','GA Convergence','Color','w');
axC = axes('Parent',figC); hold(axC,'on'); grid(axC,'on');
hCost = animatedline(axC,'LineWidth',1.6);
hK    = animatedline(axC,'LineWidth',1.6);
xlabel(axC,'Evaluation step'); ylabel(axC,'Value');
legend(axC, {'cost','K'}, 'Location','best');
end
function write_csv_header_if_needed(fp, header)
if ~exist(fp,'file'), fid = fopen(fp,'w'); fprintf(fid,'%s\n', header); fclose(fid); end
end
function append_csv(fp, eval_step, cost, K)
fid = fopen(fp,'a'); fprintf(fid,'%d,%.8f,%.8f\n', eval_step, cost, K); fclose(fid);
end
