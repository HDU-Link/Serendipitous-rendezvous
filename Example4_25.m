clear; close all; clc;
%% 参数设置
N = 6;T = 25;dt = 0.01;tspan = 0:dt:T;
alpha = [0.1; 2; 0.05; 0.9; 0; 0];K = 1;beta = 0.5;
%% 初始条件
global theta0_0 dz0_0 psi_prime phi
theta0_0 = pi*(1.25); dtheta0_0 = 0.5; z0_0 = 0; dz0_0 = 0.5;
theta0 = [3*pi/2; 7*pi/6; 11*pi/8; 5*pi/3; 7*pi/4; -pi/6+2*pi];
dtheta0 = [1;-1;2;-2;3;1];
z0 = [-1/2;1/2;-1;1;3/4;-2/3];
dz0 = [0.5;-0.1;0.2;-0.3;-0.2;1];
%% 初始化状态向量
y0 = [theta0; dtheta0; z0; dz0];
state_history = zeros(length(y0), length(tspan));
%% 定义函数
r = pi; psi_prime = @(s) s;
phi = @(s) K ./ (1 + s.^2).^beta;
%% 数值积分
y = y0; state_history(:, 1) = y0;
for k = 1:length(tspan)-1
    t = tspan(k);
    k1 = dt * system_dynamics(t, y, alpha);
    k2 = dt * system_dynamics(t + dt/2, y + k1/2, alpha);
    k3 = dt * system_dynamics(t + dt/2, y + k2/2, alpha);
    k4 = dt * system_dynamics(t + dt, y + k3, alpha);
    y = y + (k1 + 2*k2 + 2*k3 + k4)/6;
    state_history(:, k+1) = y;
end
%% 提取结果
theta = state_history(1:N, :);dtheta = state_history(N+1:2*N, :);
z = state_history(2*N+1:3*N, :);dz = state_history(3*N+1:4*N, :);
theta_leader = 0.5*tspan + theta0_0;z_leader = tspan/2;
%% 可视化结果
figure('Position', [100, 100, 1000, 400]);set(gcf,'color',[1 1 1 0]);
time_points = [6, 12, 25];
colors = {'red', 'blue', 'magenta', 'green', 'cyan', [255 96 0]/255};
for i = 1:3
    subplot(1, 3, i);hold on;
    r = 1;z_min = -5*i-1;z_max = 5*i+1;    tmp = linspace(0, 2*pi, 80);
    z_side = linspace(z_min, z_max, 3);zlim([-5,z_max-1]);
    [Theta, Z_side] = meshgrid(tmp, z_side);
    X_side = r * cos(Theta);Y_side = r * sin(Theta);
    surf(X_side, Y_side, Z_side, 'FaceColor', 'white', 'Linewidth',1,...
    'FaceAlpha', 0.8, 'EdgeColor', [11,188,182]/255,'EdgeAlpha',0.6);
    T_end = time_points(i);
    idx_end = find(tspan <= T_end, 1, 'last');
    for j = 1:N
        x_follower = cos(theta(j, 1:idx_end));
        y_follower = sin(theta(j, 1:idx_end));
        plot3(x_follower, y_follower, z(j, 1:idx_end), 'Color', colors{j}, 'LineWidth', 1.5);
        plot3(x_follower(1), y_follower(1), z(j, 1), 'k*', 'MarkerSize', 4, 'Color', colors{j});
        plot3(x_follower(end), y_follower(end), z(j, idx_end), 'ko', 'MarkerSize', 4, 'Color', colors{j});
    end
    x_leader = cos(theta_leader(1:idx_end));y_leader = sin(theta_leader(1:idx_end));
    plot3(x_leader, y_leader, z_leader(1:idx_end), 'k:', 'LineWidth', 2);
    plot3(x_leader(1), y_leader(1), z_leader(1), 'k*', 'MarkerSize', 4);
    plot3(x_leader(end), y_leader(end), z_leader(idx_end), 'ko', 'MarkerSize', 4, 'Color', 'k');
    xlabel('$x$','interpreter','latex');ylabel('$y$','interpreter','latex');
    zlabel('$z$','interpreter','latex');grid off;
    title(strcat('$t \in [0,$', num2str(T_end),'$]$'),'interpreter','latex');view(340, 30);
end
%% 角速度和z方向变化率图像
figure('Position', [100, 100, 800, 400]);set(gcf,'color',[1 1 1 0]);
subplot(1, 2, 1);
for i = 1:6
    plot(tspan, dtheta(i,:), 'LineWidth', 1, 'color', colors{i});hold on;
end
plot(tspan, dtheta0_0*ones(size(tspan)), 'k-', 'LineWidth', 1);
xlabel('$t$','interpreter','latex');
ylabel('$\dot{\theta}_i$','interpreter','latex','rotation',0);
legend('$\dot{\theta}_1$','$\dot{\theta}_2$','$\dot{\theta}_3$',...
    '$\dot{\theta}_4$','$\dot{\theta}_5$','$\dot{\theta}_6$',...
    '$\dot{\theta}_0$', 'Interpreter', 'latex','NumColumns',2,'location','se');
subplot(1, 2, 2);
for i = 1:6
    plot(tspan, dz(i,:), 'LineWidth', 1, 'color', colors{i});hold on;
end
plot(tspan, dz0_0*ones(size(tspan)), 'k-', 'LineWidth', 1);
xlabel('$t$','interpreter','latex');
ylabel('$\dot{z}_i$','interpreter','latex','rotation',0);
legend('$\dot{z}_1$','$\dot{z}_2$','$\dot{z}_3$',...
    '$\dot{z}_4$','$\dot{z}_5$','$\dot{z}_6$',...
    '$\dot{z}_0$', 'Interpreter', 'latex','NumColumns',2,'location','se');
%% 系统动力学函数
function dydt = system_dynamics(t, y, alpha)
    global theta0_0 dz0_0 psi_prime phi;
    N = length(alpha);
    theta = y(1:N);    dtheta = y(N+1:2*N);
    z = y(2*N+1:3*N);    dz = y(3*N+1:4*N);
    theta_leader = 0.5*t + theta0_0;    dtheta_leader = 0.5;
    ddtheta_leader = 0;    z_leader = t*0.5;    ddz_leader = 0;
    ddtheta = zeros(N, 1);    ddz = zeros(N, 1);
    for i = 1:N
        d_i = sqrt((theta_leader - theta(i))^2 + (z_leader - z(i))^2);
        coupling_theta = 0;        coupling_z = 0;
        for j = 1:N
            if i ~= j
                d_ij = sqrt((theta(j) - theta(i))^2 + (z(j) - z(i))^2);
                phi_ij = phi(d_ij);
                coupling_theta = coupling_theta + phi_ij * (dtheta(j) - dtheta(i));
                coupling_z = coupling_z + phi_ij * (dz(j) - dz(i));
            end
        end
        u_theta = alpha(i) * (dtheta_leader - dtheta(i)) + ddtheta_leader + ...
                 (psi_prime(d_i)/d_i) * (theta_leader - theta(i));
        u_z = alpha(i) * (dz0_0 - dz(i)) + ddz_leader + ...
              (psi_prime(d_i)/d_i) * (z_leader - z(i));
        ddtheta(i) = coupling_theta + u_theta;
        ddz(i) = coupling_z + u_z;
    end
    dydt = [dtheta; ddtheta; dz; ddz];
end