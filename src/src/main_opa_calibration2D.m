function opa_ga_peakfocus_pointing_live_imaq_2D()
% ============================================================
% OPA 初始相位校准 —— GA（遗传算法）闭环（实时采图）
%
% 本版本（按你的最新要求）：
%   1) 评价函数仅使用 K（余弦相似度），GA 最小化 cost = 1 - K
%   2) 使用 2D 窗口（矩形 ROI 或 ROI 内子矩形 SUBRECT）计算 K
%   3) 参考模板 H(x,y) 来自“相位差=0° 的物理曲线”扩展到 2D：
%        - x 方向：H_x(x) = 0°物理远场强度曲线（角度域）重采样到 ROI 宽度 Lx
%        - y 方向：H_y(y) 默认均匀（uniform），可选 gaussian（工程上可调）
%        - H2D(y,x) = H_y(y) * H_x(x)^T   （外积得到2D模板）
%   4) 中心不写死；模板由对称角度轴映射自动居中（Lx变化自动适配）
%   5) 电压范围缩小（默认0~5V），并支持电压量化步进 V_STEP
%
% 依赖：
%   - 你已有的 VC96_setVone(FPGA, ch, v) 串口写电压函数
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

% ---------- 电压范围（建议按你说的2π≈5V） ----------
VMIN        = 0;
VMAX        = 5.0;

% ---------- 电压量化步进（避免“错过极值”的工程化控制） ----------
V_STEP      = 0.02;        % 20 mV，可按噪声改成 0.05

% ---------- 写入 ramp ----------
RAMP_DV     = 0.05;
RAMP_DT     = 0.08;
SETTLE_T    = 0.12;

% ---------- 2D ROI（整幅图坐标） ----------
% 注意：MATLAB 1-based，建议 x 从 1 开始写
RECT = [1, 250, 1280, 120];   % [x,y,w,h] 你现在横向1280；h按你感兴趣区域调
LOCK_RECT_FROM_FIRST = true;

% ---------- 2D K 计算窗口模式 ----------
K2D_MODE = 'roi';      % 'roi' 或 'subroi'
% 若用 subroi：SUBRECT 是 ROI 内坐标（1-based） [sx, sy, sw, sh]
% 例如只取 ROI 中部区域（你可按兴趣改）
SUBRECT  = [1, 1, 1280, 120];

% ---------- 2D 模板 H 的 y 向扩展 ----------
HY_MODE      = 'uniform';   % 'uniform' 或 'gaussian'
HY_SIGMA_PX  = 40;          % HY_MODE='gaussian' 时生效（像素）

% ---------- 预处理（对 G 与 H 同一处理，以实现 K 定义一致） ----------
PP_REMOVE_BASE  = true;
PP_L2_NORM      = true;

% ---------- GA 参数 ----------
GA_POP        = 30;
GA_GENS       = 25;
GA_ELITE      = 2;
GA_CROSS_FRAC = 0.80;
GA_MUT_FCN    = @mutationadaptfeasible;
GA_AVG_FRAMES = 1;        % 每次评估平均帧数：2~3更稳但更慢

% ---------- 初始电压 ----------
V_INIT = zeros(1, numel(CH_LIST));

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
SAVE_DIR = fullfile(pwd,'ga_peakfocus_logs_K2D_physH');
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

%% ============ 首帧采图（0V） ============
[I0, g, hImg, fail_count] = capture_frame_live_safe_minimal(g, hImg, fail_count);
I0 = double_gray(I0);
[H0,W0] = size(I0);

rect = crop_rect_to_image(RECT, W0, H0);
if LOCK_RECT_FROM_FIRST
    rect_locked = rect;
else
    rect_locked = [];
end

% 显示：首帧 + ROI 红框
show_first_frame_with_roi(I0, rect);

% 提取首帧 patch（用于展示，不参与模板）
[patch0, rect_used0] = extract_roi_patch_from_frame(I0, rect, K2D_MODE, SUBRECT);

% ====== 构造 2D 参考模板 H2D（物理0°曲线扩展到2D） ======
Lx = size(patch0, 2);  % ROI(或subROI)宽度
Ly = size(patch0, 1);  % ROI(或subROI)高度
[H2D, Hinfo] = build_H2D_from_phys0curve(Lx, Ly, HY_MODE, HY_SIGMA_PX);

% 预处理成向量，作为 refProfile
refProfile = preprocess_for_cosine_2d(H2D, PP_REMOVE_BASE, PP_L2_NORM);

fprintf('[H2D] built from phys(phase=0): Lx=%d, Ly=%d, HY_MODE=%s\n', ...
        Lx, Ly, HY_MODE);

% ====== 可视化：首帧 patch 与 H2D ======
figPatch = figure('Name','2D Patch & H2D Template (Init)','Color','w');
subplot(1,2,1); imagesc(patch0); axis image; title('Init patch (G)'); colorbar;
subplot(1,2,2); imagesc(H2D);    axis image; title('Template H2D (phys 0° expanded)'); colorbar;

% ====== 收敛曲线（cost & K） ======
[figC, axC, hCost, hK] = init_convergence_plot_K();
write_csv_header_if_needed(LOG_CSV, 'eval,cost,K');
eval_step = 0;

% 初始化运行时句柄
init_runtime_handles(FPGA, g, hImg, fail_count, CH_LIST, ...
                     VMIN, VMAX, V_STEP, RAMP_DV, RAMP_DT, SETTLE_T);

% 初始评价
[c0, K0] = evaluate_once_2d(V, refProfile, ...
    rect, rect_locked, LOCK_RECT_FROM_FIRST, ...
    K2D_MODE, SUBRECT, ...
    PP_REMOVE_BASE, PP_L2_NORM, ...
    GA_AVG_FRAMES);

eval_step = eval_step + 1;
addpoints(hCost, eval_step, c0);
addpoints(hK,    eval_step, K0);
title(axC, sprintf('Init: cost=%.4f | K=%.4f', c0, K0), 'FontWeight','normal');
drawnow limitrate;
append_csv(LOG_CSV, eval_step, c0, K0);

fprintf('[INIT] cost=%.4f | K=%.4f\n', c0, K0);

%% ============ GA 主优化 ============
fprintf('\n==== Entering GA optimization (2D cost = 1 - K) ====\n');

lb = VMIN * ones(1, N);
ub = VMAX * ones(1, N);

fitnessFcn = @(Vrow) ga_cost(Vrow);

opts = optimoptions('ga', ...
    'PopulationSize', GA_POP, ...
    'MaxGenerations', GA_GENS, ...
    'EliteCount', GA_ELITE, ...
    'CrossoverFraction', GA_CROSS_FRAC, ...
    'MutationFcn', GA_MUT_FCN, ...
    'UseParallel', false, ...
    'Display', 'iter', ...
    'PlotFcn', {@gaplotbestf});

[Vbest_row, bestCost] = ga(fitnessFcn, N, [],[],[],[], lb, ub, [], opts);
Vbest = Vbest_row(:);
Vbest = clamp_vec(Vbest, VMIN, VMAX);
Vbest = quantize_voltage(Vbest, VMIN, V_STEP);

fprintf('\n[GA] Done. bestCost=%.6f\n', bestCost);
fprintf('[GA] Vbest=['); fprintf(' %.5f', Vbest); fprintf(' ]\n');

% 写入最佳电压并采最终帧
Hrun = get_runtime_handles();
write_vec_soft(Hrun.FPGA, Hrun.CH_LIST, Vbest, Hrun.VMIN, Hrun.VMAX, Hrun.RAMP_DV, Hrun.RAMP_DT);
pause(Hrun.SETTLE_T);

[I_final, Hrun.g, Hrun.hImg, Hrun.fail_count] = capture_frame_live_safe_minimal(Hrun.g, Hrun.hImg, Hrun.fail_count);
set_runtime_handles(Hrun);

I_final = double_gray(I_final);

if LOCK_RECT_FROM_FIRST
    rect_eval = rect_locked;
else
    rect_eval = crop_rect_to_image(RECT, size(I_final,2), size(I_final,1));
end

% 最终 patch（与模板维度一致的窗口）
[patchF, rect_usedF] = extract_roi_patch_from_frame(I_final, rect_eval, K2D_MODE, SUBRECT);

% 最终 K（单次）
GnF = preprocess_for_cosine_2d(patchF, PP_REMOVE_BASE, PP_L2_NORM);
K_final = cosine_similarity_Kvec(GnF, refProfile);

% 可视化：Init/Final patch 对比
figure('Name','2D Patch Compare (Init vs Final)','Color','w');
subplot(1,2,1); imagesc(patch0); axis image; title('Init patch (G)');  colorbar;
subplot(1,2,2); imagesc(patchF); axis image; title(sprintf('Final patch (G), K=%.4f', K_final)); colorbar;

% 最终帧图像+ROI红框 + 若subroi则显示子框
figure('Name','Final frame + ROI','Color','w');
imshow(I_final, [], 'Border','tight'); hold on;
rectangle('Position', rect_eval, 'EdgeColor','r', 'LineWidth', 1.8);
if strcmpi(K2D_MODE,'subroi')
    % 在 full-image 坐标画 subroi
    rectangle('Position', rect_usedF, 'EdgeColor','y', 'LineWidth', 1.6, 'LineStyle','--');
    legend({'ROI','subROI'},'TextColor','w','Location','southoutside');
end
title('Final frame (ROI red, subROI yellow dashed)','FontWeight','normal');
hold off;

% 保存结果
Hrun = get_runtime_handles();
fprintf('\n[Camera] snapshot 失败总次数：%d\n', Hrun.fail_count);

save(fullfile(SAVE_DIR,'final_ga_K2D_physH.mat'), ...
    'Vbest','bestCost','RECT','rect_eval','K2D_MODE','SUBRECT', ...
    'VMIN','VMAX','V_STEP','RAMP_DV','RAMP_DT','SETTLE_T', ...
    'PP_REMOVE_BASE','PP_L2_NORM','Hinfo','HY_MODE','HY_SIGMA_PX', ...
    'K0','K_final','rect_used0','rect_usedF');

writematrix(Vbest, fullfile(SAVE_DIR,'final_voltage_ga_K2D_physH.csv'));

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
try, clear g; end %#ok<TRYNC>
try, fclose(FPGA); end %#ok<TRYNC>
try, delete(instrfindall); end %#ok<TRYNC>
fprintf('> Port closed.\n');

%% ============================================================
%                 GA 代价函数（嵌套）
% ============================================================
function cost = ga_cost(Vrow)
    Vcand = clamp_vec(Vrow(:), VMIN, VMAX);
    Vcand = quantize_voltage(Vcand, VMIN, V_STEP);

    [cost, K] = evaluate_once_2d(Vcand, refProfile, ...
        rect, rect_locked, LOCK_RECT_FROM_FIRST, ...
        K2D_MODE, SUBRECT, ...
        PP_REMOVE_BASE, PP_L2_NORM, ...
        GA_AVG_FRAMES);

    eval_step = eval_step + 1;
    addpoints(hCost, eval_step, cost);
    addpoints(hK,    eval_step, K);
    title(axC, sprintf('eval=%d | cost=%.4f | K=%.4f', eval_step, cost, K), 'FontWeight','normal');
    drawnow limitrate;

    append_csv(LOG_CSV, eval_step, cost, K);

    if ~isfinite(cost), cost = 1.0; end
end

end % end main function

%% ============================================================
%           单次评价（写电压 + 采图 + 2D K）
% ============================================================
function [cost, K] = evaluate_once_2d( ...
    Vvec, refProfile, rect, rect_locked, lock_rect, ...
    k2d_mode, subrect_in_roi, ...
    pp_remove_base, pp_l2_norm, ...
    avg_frames)

H = get_runtime_handles();

Vvec = clamp_vec(Vvec(:), H.VMIN, H.VMAX);
Vvec = quantize_voltage(Vvec, H.VMIN, H.V_STEP);

write_vec_soft(H.FPGA, H.CH_LIST, Vvec, H.VMIN, H.VMAX, H.RAMP_DV, H.RAMP_DT);
pause(H.SETTLE_T);

K_acc = 0;

for t = 1:avg_frames
    [I, H.g, H.hImg, H.fail_count] = capture_frame_live_safe_minimal(H.g, H.hImg, H.fail_count);
    I = double_gray(I);

    if lock_rect
        rect_eval = rect_locked;
    else
        rect_eval = crop_rect_to_image(rect, size(I,2), size(I,1));
    end

    [patch, ~] = extract_roi_patch_from_frame(I, rect_eval, k2d_mode, subrect_in_roi);
    Gn = preprocess_for_cosine_2d(patch, pp_remove_base, pp_l2_norm);

    Kk = cosine_similarity_Kvec(Gn, refProfile);
    K_acc = K_acc + Kk;
end

K = K_acc / avg_frames;
cost = 1 - K;

set_runtime_handles(H);
end

%% ============================================================
%                 运行时句柄管理
% ============================================================
function H = get_runtime_handles()
persistent HH
if isempty(HH)
    error('Runtime handles not initialized.');
end
H = HH;
end

function set_runtime_handles(H)
persistent HH
HH = H;
end

function init_runtime_handles(FPGA, g, hImg, fail_count, CH_LIST, ...
                              VMIN, VMAX, V_STEP, RAMP_DV, RAMP_DT, SETTLE_T)
H = struct();
H.FPGA = FPGA;
H.g = g;
H.hImg = hImg;
H.fail_count = fail_count;
H.CH_LIST = CH_LIST;
H.VMIN = VMIN; H.VMAX = VMAX;
H.V_STEP = V_STEP;
H.RAMP_DV = RAMP_DV; H.RAMP_DT = RAMP_DT;
H.SETTLE_T = SETTLE_T;
set_runtime_handles(H);
end

%% ============================================================
%                 收敛曲线 UI（cost & K）
% ============================================================
function [figC, axC, hCost, hK] = init_convergence_plot_K()
figC = figure('Name','GA Convergence (2D cost & K)','Color','w');
axC = axes('Parent',figC); hold(axC,'on'); grid(axC,'on');
hCost = animatedline(axC,'LineWidth',1.6);
hK    = animatedline(axC,'LineWidth',1.6);
xlabel(axC,'Evaluation step');
ylabel(axC,'Value');
legend(axC, {'cost = 1-K (lower better)','K (higher better)'}, 'Location','best');
title(axC,'GA convergence');
end

%% ============================================================
%                 CSV 日志
% ============================================================
function write_csv_header_if_needed(fp, header)
if ~exist(fp,'file')
    fid = fopen(fp,'w');
    fprintf(fid,'%s\n', header);
    fclose(fid);
end
end

function append_csv(fp, eval_step, cost, K)
fid = fopen(fp,'a');
fprintf(fid,'%d,%.8f,%.8f\n', eval_step, cost, K);
fclose(fid);
end

%% ============================================================
%   2D 模板：0°物理曲线 Hx(x) 扩展到二维 H2D(y,x)
% ============================================================
function [H2D, info] = build_H2D_from_phys0curve(Lx, Ly, hy_mode, hy_sigma_px)
    % 1) 生成 0°物理曲线（角度域）并重采样到 x 像素长度 Lx
    Hx = build_Hx_from_slit_model_phase0(Lx);

    % 2) y 向扩展
    switch lower(hy_mode)
        case 'uniform'
            Hy = ones(Ly,1);
        case 'gaussian'
            c = round((Ly+1)/2);
            yy = (1:Ly).';
            Hy = exp(-0.5*((yy - c)/(hy_sigma_px+1e-12)).^2);
        otherwise
            error('Unknown HY_MODE: %s', hy_mode);
    end

    % 3) 外积得到 2D 模板：Ly x Lx
    H2D = Hy * (Hx(:).');  % (Ly,1)*(1,Lx)

    % 4) 归一化到峰值=1（可选），之后预处理时还会去基线/L2
    H2D = H2D / (max(H2D(:)) + 1e-12);

    info.Lx = Lx;
    info.Ly = Ly;
    info.HY_MODE = hy_mode;
    info.center_x = round((Lx+1)/2);
    info.center_y = round((Ly+1)/2);
end

function Hx = build_Hx_from_slit_model_phase0(L)
    % 物理参数（按你给的参考代码）
    lambda = 1550e-9;
    slit_width = 500e-9;
    d = 775e-9;
    N = 32;

    theta_max = 90;       % deg
    theta_points = 5000;

    theta = linspace(-theta_max, theta_max, theta_points) * pi/180; % rad
    theta_deg = theta * 180/pi;
    k = 2*pi/lambda;

    % 单缝包络
    beta = (k * slit_width/2) .* sin(theta);
    single_slit_factor = (sin(beta)./beta).^2;
    single_slit_factor(isnan(single_slit_factor)) = 1;

    % 多缝干涉：phase_diff=0°
    field = zeros(size(theta));
    for j = 1:N
        y_pos = (j-1) * d;
        field = field + exp(1i * (k * y_pos * sin(theta))); % phase=0
    end
    multi_slit_factor = abs(field).^2 / N^2;

    total_intensity = single_slit_factor .* multi_slit_factor;
    total_intensity = total_intensity / (max(total_intensity) + 1e-12);

    % 映射到像素：-90..90 等间隔对应 1..L
    theta_pix = linspace(-theta_max, theta_max, L);
    Hx = interp1(theta_deg, total_intensity, theta_pix, 'linear', 'extrap');
    Hx = Hx(:);
    Hx = Hx / (max(Hx) + 1e-12);
end

%% ============================================================
%                 2D patch 提取 / 预处理 / K
% ============================================================
function [patch, rect_used] = extract_roi_patch_from_frame(I, rectROI, k2d_mode, subrect_in_roi)
    I = double_gray(I);

    x1 = rectROI(1); y1 = rectROI(2); w = rectROI(3); h = rectROI(4);
    ROI = I(y1:y1+h-1, x1:x1+w-1);

    if strcmpi(k2d_mode,'roi')
        patch = ROI;
        rect_used = rectROI;
        return;
    end

    sr = subrect_in_roi;
    sx = round(sr(1)); sy = round(sr(2)); sw = round(sr(3)); sh = round(sr(4));

    sx = max(1, sx); sy = max(1, sy);
    sx2 = min(w, sx + sw - 1);
    sy2 = min(h, sy + sh - 1);

    patch = ROI(sy:sy2, sx:sx2);
    rect_used = [x1+sx-1, y1+sy-1, size(patch,2), size(patch,1)];
end

function v = preprocess_for_cosine_2d(P, remove_base, do_l2)
    P = double(P);
    if remove_base
        P = P - min(P(:));
        P(P < 0) = 0;
    end
    v = P(:);
    if do_l2
        v = v / (norm(v) + 1e-12);
    end
end

function K = cosine_similarity_Kvec(a,b)
    a = a(:); b = b(:);
    L = min(numel(a), numel(b));
    a = a(1:L); b = b(1:L);
    K = (a' * b) / ((norm(a)+1e-12) * (norm(b)+1e-12));
    K = max(min(K,1),-1);
end

%% ============================================================
%                 指向模块（保留）
% ============================================================
function [V_set, dphi_adj, dphi_n, dP_mW] = apply_theta_pointing( ...
        V_base, THETA_DEG, lambda, d, P_pi_mW, R_ohm, WRAP_PHASE, USE_PUSH_PULL, ...
        VMIN, VMAX, RAMP_DV, RAMP_DT, FPGA, CH_LIST)

N = numel(CH_LIST);
theta_rad  = deg2rad_local(THETA_DEG);
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
    soft_ramp_set(FPGA, CH_LIST(i), V_set(i), VMIN, VMAX, RAMP_DV, RAMP_DT);
end
end

%% ============================================================
%                 相机采图（含重连）
% ============================================================
function [I, g, hImg, fail_count] = capture_frame_live_safe_minimal(g, hImg, fail_count)
    while true
        try
            if isempty(g) || ~isvalid(g)
                fprintf('[Camera] 创建 gigecam 对象...\n');
                g = gigecam;
                if isprop(g,'Timeout'), g.Timeout = 2.0; end
            end

            Iraw = snapshot(g);

            if isempty(hImg) || ~isvalid(hImg)
                figure('Name','Live Camera (gigecam)', 'Color','k','NumberTitle','off');
                hImg = imshow(double_gray(Iraw), [], 'Border','tight');
                title('Live','Color','w');
                drawnow;
            else
                set(hImg, 'CData', Iraw);
                drawnow limitrate;
            end

            I = Iraw;
            return;

        catch ME
            fail_count = fail_count + 1;
            fprintf('[Camera] snapshot 失败 %d 次：%s\n', fail_count, ME.message);
            try, clear g; end %#ok<TRYNC>
            g = [];
            pause(0.2);
            fprintf('[Camera] 将重新创建 gigecam 并重试 snapshot...\n');
        end
    end
end

%% ============================================================
%                 基础工具函数
% ============================================================
function D = double_gray(I)
if ndims(I) == 3, I = rgb2gray(I); end
if isa(I,'double'), D = I; else, D = double(I); end
end

function rect = crop_rect_to_image(RECT, W, H)
x1 = round(RECT(1)); y1 = round(RECT(2)); w = round(RECT(3)); h = round(RECT(4));
x1 = max(1, x1); y1 = max(1, y1);
x2 = min(W, x1 + w - 1); y2 = min(H, y1 + h - 1);
w  = max(0, x2 - x1 + 1); h  = max(0, y2 - y1 + 1);
if w <= 0 || h <= 0, error('ROI 越界，请调整 RECT。'); end
rect = [x1, y1, w, h];
end

function V = clamp_vec(V, vmin, vmax)
V = min(vmax, max(vmin, V));
end

function Vq = quantize_voltage(V, vmin, vstep)
V = double(V(:));
if vstep <= 0
    Vq = V;
    return;
end
Vq = vmin + round((V - vmin)/vstep) * vstep;
end

function v_final = soft_ramp_set(FPGA, ch, v_target, vmin, vmax, dv, dt)
v_target = min(vmax, max(vmin, v_target));
persistent Vcache
if isempty(Vcache), Vcache = zeros(1, 128); end
v0 = Vcache(ch);

if abs(v_target - v0) <= dv
    VC96_setVone(FPGA, ch, v_target);
    Vcache(ch) = v_target;
    v_final = v_target;
    return;
end

n  = ceil(abs(v_target - v0)/dv);
vv = linspace(v0, v_target, n+1);
for k = 2:numel(vv)
    VC96_setVone(FPGA, ch, vv(k));
    Vcache(ch) = vv(k);
    pause(dt);
end
v_final = Vcache(ch);
end

function write_vec_soft(FPGA, CH_LIST, V_target, vmin, vmax, dv, dt)
V_target = clamp_vec(V_target(:), vmin, vmax);
for i = 1:numel(CH_LIST)
    soft_ramp_set(FPGA, CH_LIST(i), V_target(i), vmin, vmax, dv, dt);
end
end

function y = deg2rad_local(x_deg)
y = pi/180 * x_deg;
end

function show_first_frame_with_roi(I, rect)
if ndims(I) == 3, I = rgb2gray(I); end
if ~isa(I,'double'), I = double(I); end
figure('Name','First frame + ROI (red)', 'Color','w', 'NumberTitle','off');
imshow(I, [], 'Border','tight'); hold on;
rectangle('Position', rect, 'EdgeColor', 'r', 'LineWidth', 1.8);
title(sprintf('ROI=[%d,%d,%d,%d]', rect(1), rect(2), rect(3), rect(4)), 'FontWeight','normal');
hold off; drawnow;
end
