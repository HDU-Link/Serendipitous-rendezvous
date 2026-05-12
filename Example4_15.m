clear; close all; clc;
%% 参数设置
N = 10;                                                                     % 跟随者数量
alpha = [1, 1, zeros(1, N-2)];                                              % 控制增益
T = 30;                                                                     % 仿真时间
dt = 0.01;                                                                  % 时间步长
t = 0:dt:T;                                                                 % 时间向量
phi = @(s) exp(-s);                                                         % 权重函数
r = pi; zeta = 100;
Psi_func = @(s) s^2 / (r-s+r^2/zeta);                                       % 势函数
Psi_prime_func = @(s) (2*s*(r+r^2/zeta)-s^2)/(r-s+r^2/zeta)^2;
%% 初始条件设置
% 领导者初始角度和角速度
theta_0 = @(t) 0.2*t + pi;                                                  % 领导者轨迹
dtheta0 = 0.2;                                                              % 领导者角速度
d2theta0 = 0;                                                               % 领导者角加速度
theta_init = pi*[1/6, 1/4, 1/2, 2/3, 4/5, 6/5, 19/10, 13/8, 14/9, 5/4];     % 跟随者初始角度
theta_dot = [1, 1, 2, -2, 1.7, 1, 2, 1, 1/2, -1];                           % 跟随者初始角速度
% 初始化状态变量
theta = zeros(N, length(t)); dtheta = zeros(N, length(t));                  % 跟随者角度/角速度
theta(:,1) = theta_init'; dtheta(:,1) = theta_dot';
%% 主仿真循环
for k = 1:length(t)-1
    theta0_k = theta_0(t(k));                                               % 当前时刻领导者状态
    a = zeros(N,1);
    for i = 1:N
        d_theta_i0 = theta0_k - theta(i,k);                                 % 计算与领导者的角度差
        d_i = abs(d_theta_i0);
        psi = Psi_prime_func(d_i) / d_i * d_theta_i0;
        inter = 0;                                                          % 与其他跟随者的交互项
        for j = 1:N
            if i ~= j                  
                d_ij = abs(theta(j,k) - theta(i,k));                        % 计算与邻居的角度差
                inter = inter + phi(d_ij) * (dtheta(j,k) - dtheta(i,k));
            end
        end
        a(i) = alpha(i) * (dtheta0 - dtheta(i,k)) + d2theta0 + psi + inter; % 控制输入 + 交互项
    end
    % 更新跟随者状态
    dtheta(:,k+1) = dtheta(:,k) + dt * a;
    theta(:,k+1) = theta(:,k) + dt * dtheta(:,k);
end
%% 计算领导者轨迹
theta_leader = theta_0(t);
%% 可视化结果
% 图1：角速度收敛情况
figure;colors = ['b','g','k','m','c','g','m','k','b','c'];
set(gcf,'color',[1 1 1 0]);hold on;box on;
for i = 1:N
    if(i<=5)
        plot(t, dtheta(i,:), '-', 'Color', colors(i), 'LineWidth', 1);
    else
        plot(t, dtheta(i,:), '--', 'Color', colors(i), 'LineWidth', 1);
    end
end
plot(t, dtheta0 * ones(size(t)), 'r-', 'LineWidth', 1);
xlabel('$t$','interpreter','latex');
ylabel('$\dot{\theta}_i$','interpreter','latex','rotation',0);
legend('$\dot{\theta}_1$','$\dot{\theta}_2$','$\dot{\theta}_3$',...
    '$\dot{\theta}_4$','$\dot{\theta}_5$','$\dot{\theta}_6$',...
    '$\dot{\theta}_7$','$\dot{\theta}_8$','$\dot{\theta}_9$',...
    '$\dot{\theta}_{10}$','$\dot{\theta}_0$',...
       'Interpreter', 'latex','NumColumns',3,'location','se');
% 图2：单位圆上的轨迹动画
figure;set(gcf,'color',[1 1 1 0]);
theta_circle = linspace(0, 2*pi, 100);
x_circle = cos(theta_circle);y_circle = sin(theta_circle);
time_points = [0, 10, 20, 30];
markers = {'o', 'd'};colors = ['b','g','k','m','c','b','g','k','m','c'];
for idx = 1:length(time_points)
    subplot(2, 2, idx);
    plot(x_circle, y_circle, 'k-', 'LineWidth', 1);    hold on;
    [~, time_idx] = min(abs(t - time_points(idx)));
    % 绘制跟随者
    for i = 1:N
        plot(cos(theta(i, time_idx)), sin(theta(i, time_idx)), ...
             markers{(theta_init(i)>pi)+1}, 'Color', colors(i), ...
             'MarkerSize', 4, 'linewidth', 1);
    end
    % 绘制领导者
    plot(cos(theta_leader(time_idx)), sin(theta_leader(time_idx)), ...
         'rp', 'MarkerSize', 6, 'linewidth', 1);
    xlabel('$x$','interpreter','latex');
    ylabel('$y$','interpreter','latex','rotation',0);
    axis equal;    xlim([-1.2, 1.2]); ylim([-1.2, 1.2]);
    title(sprintf('t = %.0f', time_points(idx)),'interpreter','latex');
end