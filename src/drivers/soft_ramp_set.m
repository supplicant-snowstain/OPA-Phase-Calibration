function v_final = soft_ramp_set(FPGA, ch, v_target, vmin, vmax, dv, dt)
% 平滑斜坡设压，兼容 serialport 的 VC96_setVone
    % 钳位
    v_target = min(vmax, max(vmin, v_target));

    % 记忆每路当前电压（用于平滑）
    persistent Vcache
    if isempty(Vcache), Vcache = zeros(1, 128); end
    v0 = Vcache(ch);

    % 一步到位
    if abs(v_target - v0) <= dv
        VC96_setVone(FPGA, ch, v_target);
        Vcache(ch) = v_target;
        v_final = v_target;
        return;
    end

    % 分段斜坡
    n  = ceil(abs(v_target - v0)/dv);
    vv = linspace(v0, v_target, n+1);
    for k = 2:numel(vv)
        VC96_setVone(FPGA, ch, vv(k));
        Vcache(ch) = vv(k);
        pause(dt);
    end
    v_final = Vcache(ch);
end
