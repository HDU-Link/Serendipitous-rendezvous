clear;clc;
%% 参数设置
N = 3; 	alpha = [1, 0, 0];  K = 10; beta = 1;
phi = @(s) K / (1 + s^2)^beta;                                              % 权函数
r = pi/2;   zeta = 20;
Psi = @(s) s^2 / (r - s + r^2/zeta);                                        % 势函数
dPsi = @(s) s*(2*r+2*r^2/zeta-s)/(r-s+r^2/zeta)^2;                          % 势函数求导
dt = 0.01;  T = 20;  steps = T/dt;                                          % 时间步长、总时间、总步数
%% 初始条件
% 领导者初始状态
x0 = @(t) [cos(t); sqrt(2)/2*sin(t); sqrt(2)/2*sin(t)];
v0 = @(t) [-sin(t); sqrt(2)/2*cos(t); sqrt(2)/2*cos(t)];
dv0 = @(t) [-cos(t); -sqrt(2)/2*sin(t); -sqrt(2)/2*sin(t)];
% 跟随者初始位置
x_init = [sqrt(3)/3, sqrt(3)/2, 1/2;
          sqrt(2)/3, 0, sqrt(3)/2;
          2/3, -1/2, 0];
% 跟随者初始速度  
v_init = [sqrt(3)/3, -1/30, sqrt(3)/6;
          0, 1/3, -1/6;
          -1/2, sqrt(3)/30, 0];
%% 存储变量
X = zeros(3, N, steps+1);    V = zeros(3, N, steps+1);                      % 位置、速度
X(:,:,1) = x_init; V(:,:,1) = v_init;
%% 数值积分
for k = 1:steps
    t = (k-1)*dt;
    x = X(:,:,k);    v = V(:,:,k);                                          % 当前状态
    x0t = x0(t);     v0t = v0(t);     dv0t = dv0(t);                        % 领导者状态
    a = zeros(3, N);
    for i = 1:N
        A = dot(v0t, x0t);                   
        B = 1 + dot(x(:,i), x0t);              
        C = x(:,i) + x0t;                      
        A_prime = dot(dv0t, x0t) + dot(v0t, v0t);  
        B_prime = dot(v(:,i), x0t) + dot(x(:,i), v0t); 
        C_prime = v(:,i) + v0t;                   
        Pi0v0 = para(v0t, x0t, x(:,i));                                     % 球面上的平行移动
        dPi0v0dt = dv0t - ((A_prime*B - A*B_prime)/(B^2) * C + (A/B) * C_prime);
        DPi0v0 = dPi0v0dt - dot(dPi0v0dt, x(:,i)) .* x(:,i);                % 协变导数
        d_i0 = acos(dot(x(:,i), x0t));
        if d_i0 >1e-5
            grad_Psi = dPsi(d_i0)/d_i0 * log(x(:,i), x0t);                  % 势函数项
        else
            grad_Psi = 0;
        end
        coupling = zeros(3,1);                                              % 耦合项
        for j = 1:N
            if i ~= j
                dij = acos(dot(x(:,i), x(:,j)));
                Pijvj = para(v(:,j), x(:,j), x(:,i));
                PijPj0v0 = para(para(v0t, x0t, x(:,j)), x(:,j), x(:,i));
                coupling = coupling + phi(dij) * (Pi0v0 - PijPj0v0);
                coupling = coupling + phi(dij) * (Pijvj - v(:,i));
            end
        end
        a(:,i) = coupling + DPi0v0 + alpha(i)*(Pi0v0 - v(:,i)) + grad_Psi;
    end
    % 更新状态（使用投影法保持球面约束）
    for i = 1:N
%        ds = norm(v(:,i)*dt);
%        x_new = cos(ds) * x(:,i) + sin(ds) * v(:,i)/norm(v(:,i));
       x_new = x(:,i) + v(:,i)*dt;                                         % 位置更新
       x_new = x_new / norm(x_new);                                        % 投影回球面
       v_new = v(:,i) + a(:,i)*dt - x(:,i)*norm(v(:,i))^2*dt;              % 速度更新        
       X(:,i,k+1) = x_new; V(:,i,k+1) = v_new;
    end
end
%% 可视化结果
% 绘制球面
figure('Position', [200, 200, 400, 300]);
[sx,sy,sz] = sphere(50);
surf(sx,sy,sz, 'FaceAlpha', 0.3, 'EdgeColor', 'none');colormap(cool);
xlabel('$$x$$','Interpreter', 'latex');ylabel('$$y$$','Interpreter', 'latex');
zlabel('$$z$$','Interpreter', 'latex');set(gcf,'color',[1 1 1 0]);hold on;
% 绘制速度、位置轨迹
v0_traj = zeros(3, steps+1);    x0_traj = zeros(3, steps+1);
for k = 1:steps+1
    v0_traj(:,k) = v0((k-1)*dt);
    x0_traj(:,k) = x0((k-1)*dt);
end
colors = ['b', 'g', 'k'];
for i = 1:N
    plot3(squeeze(X(1,i,:)), squeeze(X(2,i,:)), squeeze(X(3,i,:)), colors(i), 'LineWidth', 1);
    plot3(X(1,i,1), X(2,i,1), X(3,i,1), '*', 'Color', colors(i), 'MarkerSize', 4, 'MarkerFaceColor', colors(i));
    plot3(X(1,i,end), X(2,i,end), X(3,i,end), 'o', 'Color', colors(i),'MarkerSize', 4);
end
% 领导者轨迹
plot3(x0_traj(1,:), x0_traj(2,:), x0_traj(3,:), 'r', 'LineWidth', 2);
plot3(x0_traj(1,1), x0_traj(2,1), x0_traj(3,1), 'r*', 'MarkerSize', 4, 'MarkerFaceColor', 'r');
plot3(x0_traj(1,end), x0_traj(2,end), x0_traj(3,end), 'ro', 'MarkerSize', 4);
xticks(-1:1:1);yticks(-1:1:1);zticks(-1:1:1);
axis equal; grid on;view(160,10);
%% 绘制速度收敛情况
figure('Position', [660, 100, 400, 500]);
set(gcf,'color',[1 1 1 0]);t_vec = 0:dt:T;
subplot(3,1,1);hold on;box on;
for i = 1:N
    plot(t_vec, squeeze(V(1,i,:)), colors(i));
end
plot(t_vec, v0_traj(1,:), 'r', 'LineWidth', 1);
xlabel('$$t$$','Interpreter', 'latex');ylabel('$$v_x^i$$','Interpreter', 'latex');
title('The $x$-components velocities','Interpreter', 'latex');
legend('$$v_x^1$$','$$v_x^2$$','$$v_x^3$$','$$v_x^0$$','Interpreter', 'latex','NumColumns',2,'location','se');
subplot(3,1,2);hold on;box on;
for i = 1:N
    plot(t_vec, squeeze(V(2,i,:)), colors(i));
end
plot(t_vec, v0_traj(2,:), 'r', 'LineWidth', 1);
xlabel('$$t$$','Interpreter', 'latex');ylabel('$$v_y^i$$','Interpreter', 'latex');
title('The $y$-components velocities','Interpreter', 'latex');
legend('$$v_x^1$$','$$v_x^2$$','$$v_x^3$$','$$v_x^0$$','Interpreter', 'latex','NumColumns',2,'location','se');
subplot(3,1,3);hold on;box on;
for i = 1:N
    plot(t_vec, squeeze(V(3,i,:)), colors(i));
end
plot(t_vec, v0_traj(3,:), 'r', 'LineWidth', 1);
xlabel('$$t$$','Interpreter', 'latex');ylabel('$$v_z^i$$','Interpreter', 'latex');
legend('$$v_x^1$$','$$v_x^2$$','$$v_x^3$$','$$v_x^0$$','Interpreter', 'latex','NumColumns',2,'location','se');
title('The $z$-components velocities','Interpreter', 'latex');
legend('$$v_x^1$$','$$v_x^2$$','$$v_x^3$$','$$v_x^0$$','Interpreter', 'latex','NumColumns',2,'location','se');
%% 辅助函数
function v_log = log(x, y)
    % 球面上的对数映射
    d = acos(dot(x, y));
    if d < 1e-5
        v_log = (y - dot(x,y)*x);
    else
        v_log = (d / sin(d)) * (y - dot(x,y)*x);
    end
end
function pijvj = para(vj, xj, xi)
    % 球面上沿最短测地线从i到j的平行移动近似
    pijvj = vj - (dot(vj, xi) / (1 + dot(xj, xi))) * (xj + xi);
end