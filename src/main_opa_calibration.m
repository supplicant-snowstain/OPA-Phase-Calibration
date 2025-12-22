function opa_ga_peakfocus_pointing_live_imaq()
% ============================================================
% OPA 初始相位校准 —— GA（遗传算法）闭环（实时采图）
%
% 本版本（按你的要求更新）：
%   1) 评价函数只使用截图中的 K（余弦相似度）
%        K = (Gn' * Hn) / (||Gn||2 * ||Hn||2)
%      GA 最小化 cost = 1 - K
%
%   2) 理想参考剖面 H(n)：
%      - 先用多狭缝远场物理模型（phase_diff=0°）计算主瓣 FWHM(角度)
%      - 由 FWHM 得到 sigma(角度)
%      - 在像素域生成高斯模板，并将高斯中心固定到 ROI 横向中心
%            c0 = round((L+1)/2)   % L=640 -> 320
%
%   3) 不再在标题/输出中显示 Ewin/sigma2/sideRatio（完全移除这些监控项）
%
%   4) 电压搜索范围缩小（默认 VMAX=5V），并支持电压量化步进 V_STEP
%      - GA 变量做量化：Vcand = quantize_voltage(Vcand, VMIN, V_STEP)
%      - 写入仍用 ramp 平滑
%
% 注意：
%   - 仍依赖你已有的底层函数 VC96_setVone(FPGA, ch, v)
%   - 相机：gigecam + snapshot（含失败重连）
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

% ---------- 电压范围（缩小） ----------
VMIN        = 0;
VMAX        = 5.5;

% ---------- 电压量化步进 ----------
V_STEP      = 0.05;     % 50 mV；

% ---------- 写入 ramp ----------
RAMP_DV     = 0.05;
RAMP_DT     = 0.08;
SETTLE_T    = 0.12;

% ---------- ROI ----------
RECT = [0, 250, 640, 80];   % [x,y,w,h]
LOCK_RECT_FROM_FIRST = true;

% ---------- 剖面（横线） ----------
PROFILE_MODE    = 'line_ratio';
PROFILE_Y_RATIO = 0.33;
SMOOTH_WIN      = 3;

% ---------- K 的预处理 ----------
PP_REMOVE_BASE  = true;  % v = v - min(v), v<0置0
PP_L2_NORM      = true;  % v = v / ||v||2

% ---------- GA 参数 ----------
GA_POP        = 30;
GA_GENS       = 25;
GA_ELITE      = 2;
GA_CROSS_FRAC = 0.80;
GA_MUT_FCN    = @mutationadaptfeasible;
GA_AVG_FRAMES = 1;       % 每次评估平均帧数（抗噪，可设2~3）

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
SAVE_DIR = fullfile(pwd,'ga_peakfocus_logs_K_centerFixed');
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

% 显示：背景 + ROI红框 + 绿线
show_first_frame_with_roi_and_line_ratio(I0, rect, PROFILE_Y_RATIO);

% ====== 首帧：横线灰度曲线 ======
[prof_init_raw, x_axis, y_line_abs] = extract_line_profile_absY(I0, rect, PROFILE_Y_RATIO, SMOOTH_WIN);

figProf = figure('Name','Line Profile Compare (Init vs Final)','Color','w');
axProf  = axes('Parent',figProf); hold(axProf,'on'); grid(axProf,'on');
plot(axProf, x_axis, prof_init_raw, 'LineWidth', 1.6);
xlabel(axProf,'Pixel position (column index)');
ylabel(axProf,'Gray level (ADU)');
title(axProf, sprintf('Line profile at y=%d (ROI line): init vs final', y_line_abs), 'FontWeight','normal');

[~, pk_init_idx] = max(prof_init_raw);
pk_init_x = x_axis(pk_init_idx);
xline(axProf, pk_init_x, '--', sprintf('init peak @ x=%d', pk_init_x), ...
    'LabelVerticalAlignment','bottom');

% ====== 初始化运行时句柄（让 evaluate_once 能访问 FPGA/相机） ======
init_runtime_handles(FPGA, g, hImg, fail_count, CH_LIST, VMIN, VMAX, V_STEP, RAMP_DV, RAMP_DT, SETTLE_T);

% ====== 构造参考剖面 Hn（物理模型 -> sigma(角度) -> 像素域高斯；中心固定到 ROI 中心） ======
Lprof = rect(3);  % ROI 宽度
% % [refProfile, refInfo] = build_reference_profile_physics_gaussian_fixed_center(Lprof, PP_REMOVE_BASE, PP_L2_NORM);
[refProfile, refInfo] = build_reference_profile_physics_zero_curve_as_H(Lprof, PP_REMOVE_BASE, PP_L2_NORM);

fprintf('[REF] Center fixed: L=%d, centerPix=%d, FWHM=%.4f deg, sigma=%.4f deg\n', ...
        refInfo.L, refInfo.center_pix, refInfo.fwhm_deg, refInfo.sigma_deg);

% ====== 收敛曲线（画 cost 与 K） ======
[figC, axC, hCost, hK] = init_convergence_plot_K();
write_csv_header_if_needed(LOG_CSV, 'eval,cost,K');

eval_step = 0;

% 初始评价
[c0, K0] = evaluate_once(V, refProfile, ...
    rect, rect_locked, LOCK_RECT_FROM_FIRST, ...
    PROFILE_MODE, PROFILE_Y_RATIO, SMOOTH_WIN, ...
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
fprintf('\n==== Entering GA optimization (cost = 1 - K) ====\n');

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

% 最终剖面曲线 + 叠加对比
[prof_final_raw, x_axis2, y_line_abs2] = extract_line_profile_absY(I_final, rect_eval, PROFILE_Y_RATIO, SMOOTH_WIN);

figure(figProf);
plot(axProf, x_axis2, prof_final_raw, 'LineWidth', 1.6);

[~, pk_fin_idx] = max(prof_final_raw);
pk_fin_x = x_axis2(pk_fin_idx);
xline(axProf, pk_fin_x, '--', sprintf('final peak @ x=%d', pk_fin_x), ...
    'LabelVerticalAlignment','top');

legend(axProf, {'Init','Final'}, 'Location','best');
title(axProf, sprintf('Line profile compare (y=%d). init peak=%d, final peak=%d', y_line_abs2, pk_init_x, pk_fin_x), ...
      'FontWeight','normal');

% 最终帧图像+ROI+绿线
figure('Name','Final frame + ROI + Line','Color','w');
imshow(I_final, [], 'Border','tight'); hold on;
rectangle('Position', rect_eval, 'EdgeColor','r', 'LineWidth', 1.8);
yline_abs = rect_eval(2) + round(PROFILE_Y_RATIO*(rect_eval(4)-1));
plot([rect_eval(1), rect_eval(1)+rect_eval(3)-1], [yline_abs,yline_abs], 'g-', 'LineWidth', 1.8);
title(sprintf('Final frame, line y=%d', yline_abs), 'FontWeight','normal');
hold off;

% 保存结果
Hrun = get_runtime_handles();
fprintf('\n[Camera] snapshot 失败总次数：%d\n', Hrun.fail_count);

save(fullfile(SAVE_DIR,'final_ga_K_centerFixed.mat'), ...
    'Vbest','bestCost','rect_eval','RECT','PROFILE_Y_RATIO', ...
    'VMIN','VMAX','V_STEP','RAMP_DV','RAMP_DT','SETTLE_T', ...
    'PP_REMOVE_BASE','PP_L2_NORM','refInfo', ...
    'pk_init_x','pk_fin_x');

writematrix(Vbest, fullfile(SAVE_DIR,'final_voltage_ga_K_centerFixed.csv'));

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

    [cost, K] = evaluate_once(Vcand, refProfile, ...
        rect, rect_locked, LOCK_RECT_FROM_FIRST, ...
        PROFILE_MODE, PROFILE_Y_RATIO, SMOOTH_WIN, ...
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
%                 单次评价（写电压 + 采图 + 算 K）
% ============================================================
function [cost, K] = evaluate_once( ...
    Vvec, refProfile, rect, rect_locked, lock_rect, ...
    profile_mode, profile_y_ratio, smooth_win, ...
    pp_remove_base, pp_l2_norm, ...
    avg_frames)

H = get_runtime_handles();

% 电压量化（确保 evaluate_once 与 ga_cost 一致）
Vvec = clamp_vec(Vvec(:), H.VMIN, H.VMAX);
Vvec = quantize_voltage(Vvec, H.VMIN, H.V_STEP);

% 写入电压
write_vec_soft(H.FPGA, H.CH_LIST, Vvec, H.VMIN, H.VMAX, H.RAMP_DV, H.RAMP_DT);
pause(H.SETTLE_T);

K_acc = 0;

for k = 1:avg_frames
    [I, H.g, H.hImg, H.fail_count] = capture_frame_live_safe_minimal(H.g, H.hImg, H.fail_count);
    I = double_gray(I);

    if lock_rect
        rect_eval = rect_locked;
    else
        rect_eval = crop_rect_to_image(rect, size(I,2), size(I,1));
    end

    [prof_raw, ~] = extract_profile_from_frame(I, rect_eval, profile_mode, profile_y_ratio, smooth_win);

    Gn = preprocess_for_cosine(prof_raw, pp_remove_base, pp_l2_norm);
    Kk = cosine_similarity_K(Gn, refProfile);

    K_acc = K_acc + Kk;
end

K = K_acc / avg_frames;

% 评价函数：cost = 1 - K
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

function init_runtime_handles(FPGA, g, hImg, fail_count, CH_LIST, VMIN, VMAX, V_STEP, RAMP_DV, RAMP_DT, SETTLE_T)
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
figC = figure('Name','GA Convergence (cost & K)','Color','w');
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
%   参考剖面：物理模型 -> 主瓣 FWHM(角度) -> sigma(角度) -> 像素域高斯
%   且中心固定到 ROI 中心像素 c0 = round((L+1)/2)
% ============================================================
function [refProfile, info] = build_reference_profile_physics_zero_curve_as_H(L, rm_base, l2norm)
    % 用“相位差=0°”的物理曲线 total_intensity(theta) 直接作为 H(n)
    [Hn, info] = build_Hn_from_slit_model_phase0(L);

    % 与 Gn 预处理一致（建议：去基线 + L2）
    if rm_base
        Hn = Hn - min(Hn);
        Hn(Hn<0) = 0;
    end
    if l2norm
        Hn = Hn / (norm(Hn) + 1e-12);
    end

    refProfile = Hn;
end

function [Hn, info] = build_Hn_from_slit_model_phase0(L)
    % ===== 物理参数（按你给的）=====
    lambda = 1550e-9;      
    slit_width = 500e-9;   
    d = 775e-9;            
    N = 32;

    % 角度范围（±90°）
    theta_max = 90;        % deg
    theta_points = 5000;

    theta = linspace(-theta_max, theta_max, theta_points) * pi/180; % rad
    theta_deg = theta * 180/pi;

    k = 2*pi/lambda;

    % 单缝衍射包络
    beta = (k * slit_width/2) .* sin(theta);
    single_slit_factor = (sin(beta)./beta).^2;
    single_slit_factor(isnan(single_slit_factor)) = 1;

    % 多缝干涉：phase_diff = 0°
    field = zeros(size(theta));
    for j = 1:N
        y_pos = (j-1) * d;
        field = field + exp(1i * (k * y_pos * sin(theta))); % phase=0
    end
    multi_slit_factor = abs(field).^2 / N^2;

    % 0°曲线：总光强
    total_intensity = single_slit_factor .* multi_slit_factor;
    total_intensity = total_intensity / max(total_intensity); % 归一化到峰值=1

    % 将角度域曲线映射/重采样到像素长度 L
    % 这里使用等间隔角度 -> 等间隔像素（与之前一致）
    theta_pix = linspace(-theta_max, theta_max, L);  % deg，长度 L
    Hn = interp1(theta_deg, total_intensity, theta_pix, 'linear', 'extrap');
    Hn = Hn(:);

    % info（可选：记录中心像素与角度映射）
    info.L = L;
    info.theta_max = theta_max;
    info.center_pix = round((L+1)/2);
end

% % function [refProfile, info] = build_reference_profile_physics_gaussian_fixed_center(L, rm_base, l2norm)
% %     [Hn, info] = build_Hn_gaussian_from_slit_model_fixed_center(L);
% % 
% %     % 与 Gn 的预处理保持一致（去基线+L2）
% %     if rm_base
% %         Hn = Hn - min(Hn);
% %         Hn(Hn<0) = 0;
% %     end
% %     if l2norm
% %         Hn = Hn / (norm(Hn)+1e-12);
% %     end
% % 
% %     refProfile = Hn;
% % end
% % 
% % function [Hn, info] = build_Hn_gaussian_from_slit_model_fixed_center(L)
% %     % 物理参数（按你给的）
% %     lambda = 1550e-9;
% %     slit_width = 500e-9;
% %     d = 775e-9;
% %     N = 32;
% % 
% %     theta_max = 90;        % deg
% %     theta_points = 5000;
% % 
% %     theta = linspace(-theta_max, theta_max, theta_points) * pi/180; % rad
% %     theta_deg = theta * 180/pi;
% % 
% %     k = 2*pi/lambda;
% % 
% %     % 单缝包络
% %     beta = (k * slit_width/2) .* sin(theta);
% %     single_slit_factor = (sin(beta)./beta).^2;
% %     single_slit_factor(isnan(single_slit_factor)) = 1;
% % 
% %     % 多缝干涉：phase_diff=0°
% %     field = zeros(size(theta));
% %     for j = 1:N
% %         y_pos = (j-1) * d;
% %         field = field + exp(1i * (k * y_pos * sin(theta))); % phase_diff=0
% %     end
% %     multi_slit_factor = abs(field).^2 / N^2;
% % 
% %     total_intensity = single_slit_factor .* multi_slit_factor;
% %     total_intensity = total_intensity / max(total_intensity);
% % 
% %     % 主瓣 FWHM(角度)
% %     [~, pk] = max(total_intensity);
% %     halfmax = 0.5;
% % 
% %     li = find(total_intensity(1:pk) <= halfmax, 1, 'last');
% %     if isempty(li), li = 1; end
% %     ri_rel = find(total_intensity(pk:end) <= halfmax, 1, 'first');
% %     if isempty(ri_rel)
% %         ri = numel(theta_deg);
% %     else
% %         ri = pk + ri_rel - 1;
% %     end
% % 
% %     if li < pk
% %         left_deg = interp1(total_intensity(li:li+1), theta_deg(li:li+1), halfmax, 'linear','extrap');
% %     else
% %         left_deg = theta_deg(li);
% %     end
% %     if ri > pk
% %         right_deg = interp1(total_intensity(ri-1:ri), theta_deg(ri-1:ri), halfmax, 'linear','extrap');
% %     else
% %         right_deg = theta_deg(ri);
% %     end
% % 
% %     fwhm_deg  = right_deg - left_deg;
% %     sigma_deg = fwhm_deg / (2*sqrt(2*log(2)));
% % 
% %     % 像素域角度映射：-90..90 对应 1..L（等间隔）
% %     theta_pix = linspace(-theta_max, theta_max, L).';
% % 
% %     % 固定中心像素（L=640 -> 320）
% %     c0 = round((L+1)/2);
% %     mu_deg = theta_pix(c0);
% % 
% %     % 像素域高斯模板
% %     Hn = exp(-0.5 * ((theta_pix - mu_deg)/sigma_deg).^2);
% %     Hn = Hn(:);
% % 
% %     info.L = L;
% %     info.center_pix = c0;
% %     info.fwhm_deg = fwhm_deg;
% %     info.sigma_deg = sigma_deg;
% % end


%% ============================================================
%                 剖面提取 / 预处理 / K
% ============================================================
function [prof, x_axis] = extract_profile_from_frame(I, rectROI, profile_mode, profile_y_ratio, smooth_win)
    I = double_gray(I);

    x1 = rectROI(1); y1 = rectROI(2); w = rectROI(3); h = rectROI(4);
    ROI = I(y1:y1+h-1, x1:x1+w-1);

    switch lower(profile_mode)
        case 'mean'
            prof = mean(ROI,1);
        case 'sum'
            prof = sum(ROI,1);
        case 'max'
            prof = max(ROI,[],1);
        case 'centerline'
            prof = double(ROI(round(size(ROI,1)/2),:));
        case 'line_ratio'
            profile_y_ratio = max(0, min(1, profile_y_ratio));
            yy = round(profile_y_ratio * (size(ROI,1)-1)) + 1;
            yy = max(1,min(size(ROI,1),yy));
            prof = double(ROI(yy,:));
        otherwise
            error('Unknown profile_mode');
    end
    prof = double(prof(:));

    if smooth_win > 1
        k = ones(smooth_win,1)/smooth_win;
        prof = conv(prof, k, 'same');
    end

    x_axis = (x1:(x1+w-1))';
end

function v = preprocess_for_cosine(v, remove_base, do_l2)
    v = double(v(:));
    if remove_base
        v = v - min(v);
        v(v < 0) = 0;
    end
    if do_l2
        v = v / (norm(v) + 1e-12);
    end
end

function K = cosine_similarity_K(Gn, Hn)
    Gn = Gn(:); Hn = Hn(:);
    L = min(numel(Gn), numel(Hn));
    Gn = Gn(1:L); Hn = Hn(1:L);
    num = sum(Gn .* Hn);
    den = (norm(Gn)+1e-12) * (norm(Hn)+1e-12);
    K = num / den;
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

function show_first_frame_with_roi_and_line_ratio(I, rect, profile_y_ratio)
if ndims(I) == 3, I = rgb2gray(I); end
if ~isa(I,'double'), I = double(I); end

x1 = rect(1); y1 = rect(2); w = rect(3); h = rect(4);

profile_y_ratio = max(0, min(1, profile_y_ratio));
y_line = y1 + round(profile_y_ratio * (h-1));
y_line = max(1, min(size(I,1), y_line));

figure('Name','First frame + ROI (red) + profile line (green)', ...
    'Color','w', 'NumberTitle','off');
imshow(I, [], 'Border','tight'); hold on;
rectangle('Position', rect, 'EdgeColor', 'r', 'LineWidth', 1.8);
plot([x1, x1+w-1], [y_line, y_line], 'g-', 'LineWidth', 1.8);
title(sprintf('ROI=[%d,%d,%d,%d], profile\\_y\\_ratio=%.2f (y=%d)', ...
    x1, y1, w, h, profile_y_ratio, y_line), 'FontWeight','normal');
hold off; drawnow;
end

function [prof, x_axis, y_line_abs] = extract_line_profile_absY(I, rect, profile_y_ratio, smooth_win)
I = double_gray(I);
x1 = rect(1); y1 = rect(2); w = rect(3); h = rect(4);
x2 = x1 + w - 1;

profile_y_ratio = max(0, min(1, profile_y_ratio));
yy = round(profile_y_ratio * (h-1)) + 1;
yy = max(1, min(h, yy));
y_line_abs = y1 + yy - 1;

prof = double(I(y_line_abs, x1:x2));
prof = prof(:)';

if smooth_win > 1
    k = ones(1,smooth_win)/smooth_win;
    prof = conv(prof, k, 'same');
end

x_axis = x1:x2;
end
