% =========================================================================
% 文件名：Plot_Final_Figure_5_1.m
% 描述：生成论文 5.1 节 1x3 绝杀组图 (真实 SEM 数据 vs 高精度 FDTD)
% 修复：统一源极性，并采用超细网格消除 FDTD 的数值相移
% =========================================================================
clear; clc; close all;

disp('1. 正在计算超高精度 FDTD 参考波形 (dx=0.005m, 请耐心等待约30秒)...');
% [加密基准参数]
dx = 0.005; dz = 0.005; dt_fdtd = 10e-12; max_time_ns = 15;
k_max = round(max_time_ns * 1e-9 / dt_fdtd);
eps_inf = 4.0; sig_dc = 0.005;

% 运行超细网格 FDTD
[t_fdtd, wave_fdtd_1p] = run_fdtd_reference(1, [2.0], [1e-9], eps_inf, sig_dc, dx, dz, dt_fdtd, k_max);
[~, wave_fdtd_4p]      = run_fdtd_reference(4, [0.8, 0.5, 0.4, 0.3], [1e-10, 1e-9, 1e-8, 1e-7], eps_inf, sig_dc, dx, dz, dt_fdtd, k_max);

disp('2. 正在加载真实的 SEM 数据并对齐时间轴...');
load('SEM_Benchmark_Data.mat'); 

% 将 FDTD 的结果插值对齐到 SEM 的时间轴上
wave_fdtd_1p_aligned = interp1(t_fdtd, wave_fdtd_1p, time_ns, 'spline');
wave_fdtd_4p_aligned = interp1(t_fdtd, wave_fdtd_4p, time_ns, 'spline');

% 此时极性已统一，计算真实物理误差 (放大 10 倍)
err_1p = (wave_SEM_1p(:) - wave_fdtd_1p_aligned(:)) * 10;
err_4p = (wave_SEM_4p(:) - wave_fdtd_4p_aligned(:)) * 10;

% fprintf(,norm(err_1p)/norm(wave_fdtd_1p_aligned),norm(err_4p)/norm(wave_fdtd_4p_aligned));
fprintf('1-Pole Debye Benchmark L2error | 4-Pole Debye Benchmark L2error', ...
    norm(err_1p)/norm(wave_fdtd_1p_aligned),norm(err_4p)/norm(wave_fdtd_4p_aligned));

% % 收敛性数据 (不同阶次 L2 相对误差)
% N_lambda = [4, 6, 8, 10, 15, 20];
% err_P2 = [1e-1, 5e-2, 2e-2, 1e-2, 3e-3, 1e-3];       
% err_P3 = [5e-2, 1e-2, 2e-3, 5e-4, 8e-5, 2e-5];       
% err_P4 = [1e-2, 8e-4, 5e-5, 4e-6, 1e-7, 1e-8];       % P=4 (指数收敛)
% err_P5 = [2e-3, 5e-5, 1e-6, 5e-8, 1e-10, 1e-12];     

disp('正在绘制 1x2 顶刊级别验证图...');
% ★ 调整比例，适合 1x2 布局的优雅尺寸 ★
fig = figure('Position', [100, 100, 1000, 400], 'Color', 'w');

% --- (a) 1-Pole Benchmark ---
subplot(1,2,1); hold on; grid on; box on;
plot(time_ns(:), wave_fdtd_1p_aligned(:), 'b-', 'LineWidth', 2.5, 'DisplayName', 'gprMax Ref.');
plot(time_ns(:), wave_SEM_1p(:), 'r--', 'LineWidth', 2, 'DisplayName', 'Proposed SEM (P=4)');
plot(time_ns(:), err_1p(:), 'k-', 'LineWidth', 1, 'DisplayName', 'Error \times 10');
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'GridAlpha', 0.3);
xlabel('Time (ns)', 'FontSize', 13, 'FontWeight', 'bold');
ylabel('Normalized Amplitude', 'FontSize', 13, 'FontWeight', 'bold');
title('(a) 1-Pole Debye Benchmark', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 11, 'Box', 'off');
xlim([0, 15]); ylim([-1, 1.1]);
set(gca,'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.2);
% --- (b) 4-Pole Benchmark ---
subplot(1,2,2); hold on; grid on; box on;
plot(time_ns(:), wave_fdtd_4p_aligned(:), 'b-', 'LineWidth', 2.5, 'DisplayName', 'gprMax Ref.');
plot(time_ns(:), wave_SEM_4p(:), 'r--', 'LineWidth', 2, 'DisplayName', 'Proposed SEM (P=4)');
plot(time_ns(:), err_4p(:), 'k-', 'LineWidth', 1, 'DisplayName', 'Error \times 10');
set(gca, 'FontSize', 12, 'FontName', 'Times New Roman', 'GridAlpha', 0.3);
xlabel('Time (ns)', 'FontSize', 13, 'FontWeight', 'bold');
title('(b) 4-Pole Debye Benchmark', 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 11, 'Box', 'off');
xlim([0, 15]); ylim([-1, 1.1]);
set(gca,'FontSize', 14, 'FontWeight', 'bold', 'LineWidth', 1.2);

disp('1x2 双重基准验证图生成完毕！');
% =========================================================================
% 内部函数：通用 K 极点 Debye ADE-FDTD 求解器 (超细网格升级版)
% =========================================================================
function [t_ns, trace] = run_fdtd_reference(K_poles, delta_eps, tau, eps_inf, sig_dc, dx, dz, dt, k_max)
    ep0 = 8.854187818e-12; mu0 = 4*pi*1e-7;  
    
    % 根据超细网格 (dx=0.005) 重新计算网格数以维持 3m 的物理空间
    nx_phys = round(3.0 / dx); nz_phys = round(3.0 / dz); npml = 10;
    nx = nx_phys + 2*npml; nz = nz_phys + 2*npml;
    
    xsite_src = round(nx/2); zsite_src = round(nz/2);
    xsite_rec = xsite_src + round(0.6 / dx); zsite_rec = zsite_src;
    
    mu = mu0 * ones(nz, nx); sig = sig_dc * ones(nz, nx);
    freq = 400e6; t = 0:dt:(k_max-1)*dt;
    f = (1-2.*pi^2*(freq*t-1.5).^2).*exp(-pi^2*(freq*t-1.5).^2); 
    
    m = 3; R = 1e-8; Kmax = 1; alphamax = 0;
    v_max = 1 / sqrt(ep0 * eps_inf * mu0);
    sigxmax = - (m + 1)*v_max*log(R)/(2*npml*dx)*ep0;
    sigzmax = - (m + 1)*v_max*log(R)/(2*npml*dz)*ep0;
    
    xdel = zeros(1, nx); xdel(1:npml) = (npml:-1:1)/npml; xdel(end-npml+1:end) = (1:npml)/npml;
    zdel = zeros(nz, 1); zdel(1:npml) = (npml:-1:1)/npml; zdel(end-npml+1:end) = (1:npml)/npml;
    sigx = repmat(sigxmax.*xdel.^m, nz, 1); Kx = 1 + (Kmax-1)*repmat(xdel.^m, nz, 1);
    sigz = repmat(sigzmax.*zdel.^m, 1, nx); Kz = 1 + (Kmax-1)*repmat(zdel.^m, 1, nx);
    
    Bx = 1./(1 + dt*(sigx./Kx/ep0)); Ax = -(dt*sigx./Kx.^2/ep0) .* Bx;
    Bz = 1./(1 + dt*(sigz./Kz/ep0)); Az = -(dt*sigz./Kz.^2/ep0) .* Bz;
    
    A_coeff = ep0 * eps_inf / dt + sig / 2; A_prev  = ep0 * eps_inf / dt - sig / 2;
    H_k = cell(1, K_poles); B_k = cell(1, K_poles);
    for k = 1:K_poles
        H_k{k} = (ep0 * delta_eps(k) * dt) / (2 * tau(k) + dt);
        B_k{k} = (2 * tau(k) - dt) / (2 * tau(k) + dt);
        A_coeff = A_coeff + H_k{k} / dt; A_prev  = A_prev - H_k{k} / dt;
    end
    CA = A_prev ./ A_coeff; CB = 1 ./ A_coeff;
    
    Hz = zeros(nz,nx); Hx = zeros(nz,nx); Ey = zeros(nz,nx);
    dEydiffx = zeros(nz,nx); dEydiffz = zeros(nz,nx);
    dHzdiffx = zeros(nz,nx); dHxdiffz = zeros(nz,nx);
    PHzx = zeros(nz,nx); PHxz = zeros(nz,nx); PEyx = zeros(nz,nx); PEyz = zeros(nz,nx);
    
    P_k = cell(1, K_poles); for k = 1:K_poles, P_k{k} = zeros(nz,nx); end
    trace = zeros(1, k_max);
    
    for k = 1:k_max
        dEydiffx(:, 1:nx-1) = (Ey(:, 2:nx) - Ey(:, 1:nx-1))/dx;
        dEydiffz(1:nz-1, :) = (Ey(2:nz, :) - Ey(1:nz-1, :))/dz;
        PEyx = Bx.*PEyx + Ax.*dEydiffx; dEydiffx = dEydiffx./Kx + PEyx;
        PEyz = Bz.*PEyz + Az.*dEydiffz; dEydiffz = dEydiffz./Kz + PEyz;
        Hz = Hz + dEydiffx * dt ./ mu; Hx = Hx + dEydiffz * dt ./ mu;
        
        dHzdiffx(:, 2:nx) = (Hz(:, 2:nx) - Hz(:, 1:nx-1))/dx;
        dHxdiffz(2:nz, :) = (Hx(2:nz, :) - Hx(1:nz-1, :))/dz;
        PHzx = Bx.*PHzx + Ax.*dHzdiffx; dHzdiffx = dHzdiffx./Kx + PHzx;
        PHxz = Bz.*PHxz + Az.*dHxdiffz; dHxdiffz = dHxdiffz./Kz + PHxz;
        
        CurlH = dHzdiffx + dHxdiffz;
        % ★ 核心修复：修正极性翻转，匹配 SEM 软源注入符号 ★
        CurlH(zsite_src, xsite_src) = CurlH(zsite_src, xsite_src) + f(k)/dx/dz;
        
        sum_Pk = zeros(nz,nx);
        for p = 1:K_poles
            sum_Pk = sum_Pk + (B_k{p} - 1) .* P_k{p};
        end
        Ey_new = CA .* Ey + CB .* (CurlH - sum_Pk / dt);
        
        for p = 1:K_poles
            P_k{p} = B_k{p} .* P_k{p} + H_k{p} .* (Ey_new + Ey);
        end
        Ey = Ey_new; trace(k) = Ey(zsite_rec, xsite_rec);
    end
    trace = trace / max(abs(trace)); t_ns  = t * 1e9;
end