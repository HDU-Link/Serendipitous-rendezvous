clear; clc;
%% 参数设置
N = 5;                                                                      % 跟随者数量
T = 40;                                                                     % 仿真时间
dt = 0.01;                                                                  % 时间步长
tspan = 0:dt:T;                                                             % 时间向量
nt = length(tspan);                                                         % 时间步数
% 控制参数
alpha = [0.1, 1, 1, 0, 0];                                                  % 控制增益
%alpha = [0, 0, 0, 0, 0];
phi = @(s) 1./(1 + s.^2);                                                   % 权重函数
Psi = @(s) s.^2;                                                            % 势函数
Psi_prime = @(s) 2*s;
%% 初始条件
x0_0 = [0; 0];                                                              % 领导者初始位置
v0_0 = [cos(0); sin(0)];                                                    % 领导者初始速度
x_init = [10,9;6,-1;6,-4;-3,8;2,8]';                                        % 跟随者的初始位置
v_init = [-8,15;18,2;5,16;5,-4;-3, 6]';                                     % 跟随者的初始速度
%% 初始化状态变量
x0 = zeros(2, nt); v0 = zeros(2, nt); a0 = zeros(2, nt);                    % 领导者位置、速度、加速度
x = zeros(2, N, nt); v = zeros(2, N, nt);                                   % 跟随者位置、速度
x0(:,1) = x0_0;v0(:,1) = v0_0;x(:,:,1) = x_init;v(:,:,1) = v_init;          % 设置初始条件
%% 主仿真循环
for k = 1:nt-1
    t = tspan(k);
    x0(:,k) = [t*cos(t/2); t*sin(t/2)];
    v0(:,k) = [cos(t/2) - (t/2)*sin(t/2); 
               sin(t/2) + (t/2)*cos(t/2)];
    a0(:,k) = [ -sin(t/2) - (1/2)*sin(t/2) - (t/4)*cos(t/2);
                cos(t/2) + (1/2)*cos(t/2) - (t/4)*sin(t/2)];    
    % 计算每个跟随者的控制输入和加速度
    for i = 1:N
        inter_term = zeros(2,1);                                            % 计算与其他跟随者的交互项
        for j = 1:N
            if i ~= j
                phi_ij = phi(norm(x(:,i,k) - x(:,j,k)));
                inter_term = inter_term + phi_ij * (v(:,j,k) - v(:,i,k));
            end
        end
        % 计算控制输入
        di0 = norm(x(:,i,k) - x0(:,k));
        attract_term = (Psi_prime(di0)/di0) * (x0(:,k) - x(:,i,k));
        u_i = alpha(i) * (v0(:,k) - v(:,i,k)) + a0(:,k) + attract_term;
        % 更新跟随者状态
        a_i = inter_term + u_i;
        v(:,i,k+1) = v(:,i,k) + a_i * dt;
        x(:,i,k+1) = x(:,i,k) + v(:,i,k) * dt;
    end
end
%% 可视化结果
% 图1: 轨迹图 
figure('Position', [400, 230, 400, 400]);hold on; axis equal;box on;
% 绘制跟随者轨迹
colors = ['g', 'b', 'r', 'c', 'm'];
for i = 1:N
    plot(squeeze(x(1,i,:)), squeeze(x(2,i,:)), [colors(i) '-'], 'LineWidth', 1.5);
    plot(x(1,i,1), x(2,i,1), [colors(i) 'o'], 'MarkerSize', 6);
end
plot(x0(1,1:end-1), x0(2,1:end-1), 'k-', 'LineWidth', 1.5);
plot(x0(1,1), x0(2,1), 'ko', 'MarkerSize', 6);
xlabel('$$x$$','Interpreter', 'latex');ylabel('$$y$$','Interpreter', 'latex');
xlim([-40 40]);ylim([-40 40]);set(gcf,'color',[1 1 1 0]);
% 图2: 速度分量随时间变化
figure('Position', [50, 50, 1200, 200]);set(gcf,'color',[1 1 1 0]);
subplot(1,2,1);hold on; box on;
for i = 1:N
    plot(tspan, squeeze(v(1,i,:)), [colors(i) '-'], 'LineWidth', 1);
end
plot(tspan(1:end-1), v0(1,1:end-1), 'k-', 'LineWidth', 1);
xlabel('$$t$$','Interpreter', 'latex'); ylabel('$$v_x^i$$','Interpreter', 'latex');
title('The $x$-components velocities', 'Interpreter', 'latex');
legend('$$v_x^1$$','$$v_x^2$$','$$v_x^3$$','$$v_x^4$$','$$v_x^5$$','$$v_x^0$$',...
       'Interpreter', 'latex','NumColumns',2,'location','se');
subplot(1,2,2);hold on; box on;
for i = 1:N
    plot(tspan, squeeze(v(2,i,:)), [colors(i) '-'], 'LineWidth', 1);
end
plot(tspan(1:end-1), v0(2,1:end-1), 'k-', 'LineWidth', 1);
xlabel('$$t$$','Interpreter', 'latex'); ylabel('$$v_y^i$$','Interpreter', 'latex');
title('The $y$-components velocities', 'Interpreter', 'latex');
legend('$$v_y^1$$','$$v_y^2$$','$$v_y^3$$','$$v_y^4$$','$$v_y^5$$','$$v_y^0$$',...
       'Interpreter', 'latex','NumColumns',2,'location','se');