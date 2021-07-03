%% Задача линейного управления
%\cdot x = A(t)x + Bu + f(t)
A = {@(x) 0 , @(x) -1; @(x) 1, @(x) 0};
B = [5,0;0,0];
f = {@(x) 0 , @(x) sin(x)}';
t0 = 0; % начальное время

% P - множество допустимых управлений
% P = {x | a(x-p1)^2+b(x-p2)^2<=1}, a > 0, b > 0;
a = 1;  
b = 1;
p1=0;
p2=0;
params = [a,b,p1,p2];
save params_P params;
% \mathcal{X}_0 - начальное множество значений фазового вектора
% \mathcal{X}_0 = {x_0} - точка
x0=[3,3]';

% \mathcal{X}_1 - целевое множество значений фазового вектора
% \mathcal{X}_1 = {x | max { \alpha |x_1-k|, \beta |x_2-m| } <= q } q>0
% если альфа и бета положительны то прямоугольник, если одна не
% положительна, то полоса; если обе то вся плоскость. 
k=0;
m=0;
q=1;
alpha = 1;
beta = 1;
params1 = [k,m,q,alpha,beta];
g = @(x) params1(3)-max(params1(4)*abs(x(1)-params1(1)),params1(5)*abs(x(2)-params1(2)));
save params_X1 params1;
%сохранение параметов X1

plot(x0(1),x0(2),'r.','MarkerSize',30);
hold on
drawSet(@Square_Lebesgue,10);
% Отрисовка множеств X1 и X2

t_max=2*pi;
n = 25;
start_angle = 0;
fin_angle = 2*pi;
phi_t0 = [sin(start_angle:(fin_angle-start_angle)/n:fin_angle-(fin_angle-start_angle)/n); cos(start_angle:(fin_angle-start_angle)/n:fin_angle-(fin_angle-start_angle)/n)];
A_sp = A';

n_min=1;
t_min = t_max;

for i=1:n
    tspan_0 = [t0 t_max];
    [t,phi] = ode45(@(t,phi) odefcn(t,phi,A_sp), tspan_0, phi_t0(:,i));
    size_t = length(t);
     
    u_opt = zeros(2,size_t);
    for j=1:size_t
        u_1 = zeros(1,length(t));
        u_2 = zeros(1,length(t));
        u_1(1:length(t)) = B(1,1)*phi(1:length(t),1)+B(1,2)*phi(1:length(t),2);
        u_2(1:length(t)) = B(2,1)*phi(1:length(t),1)+B(2,2)*phi(1:length(t),2);
        
        [pho,supp_vec]=Ellipse_Lebesgue([u_1(j),u_2(j)]);
        u_opt(1,j)=supp_vec(1);
        u_opt(2,j)=supp_vec(2);
    end% Находим оптимальное управление в момент времени t_i
    tspan_1 = [t0 t_max];
    opts = odeset('Events',@StopEvents,'RelTol',1e-6,'AbsTol',1e-7,'Refine',6);
    [t_opt,psi_opt,te,ye,ie] = ode45(@(t_opt,psi_opt) odefcn_OC(t_opt,psi_opt,A,B,f,t,u_opt), tspan_1, x0,opts);
    if length(te)==1
        if te<t_min
            t_min = te;
            n_min = i;
        end
    end
    if (g(x0)>=0)
        t_min=t0;
        psi_opt = zeros(length(t_opt));
    end
    plot(psi_opt(:,1),psi_opt(:,2),'g');
end% Решаем диффур %\cdot x = A(t)x + Bu + f(t) и получаем кандидатов на оптимальную траекторию



tspan_0 = [t0 t_max];
[t_uopt,phi] = ode45(@(t,phi) odefcn(t,phi,A_sp), tspan_0, phi_t0(:,n_min));
size_t = length(t_uopt);

u_opt = zeros(2,size_t);
for j=1:size_t
    u_1 = zeros(1,length(t_uopt));
    u_2 = zeros(1,length(t_uopt));
    u_1(1:length(t_uopt)) = B(1,1)*phi(1:length(t_uopt),1)+B(1,2)*phi(1:length(t_uopt),2);
    u_2(1:length(t_uopt)) = B(2,1)*phi(1:length(t_uopt),1)+B(2,2)*phi(1:length(t_uopt),2);

    [pho,supp_vec]=Ellipse_Lebesgue([u_1(j),u_2(j)]);
    u_opt(1,j)=supp_vec(1);
    u_opt(2,j)=supp_vec(2);
end% Находим оптимальное управление в момент времени t_i


tspan_1 = [t0 t_max];
opts_ = odeset('Events',@StopEvents,'RelTol',1e-8,'AbsTol',1e-9,'Refine',5);
[t_opt,psi_opt,te,ye,ie] = ode45(@(t_opt,psi_opt) odefcn_OC(t_opt,psi_opt,A,B,f,t,u_opt), tspan_1, x0,opts_);
if (t0==t_min)
    psi_opt = zeros(length(t_opt));
end

if t_min < t_max-0.1
    if t_min==t0
        t = t0;
        phi = phi_t0(:,n_min);
        length(t)
        sk_pr = -(phi(1,1)*psi_opt(length(t_opt),1)+phi(2,1)*psi_opt(length(t_opt),2));
        [op_func,hhhhh]=Square_Lebesgue(-[phi(length(t),1),phi(2,1)]);        
    else
        [t,phi] = ode45(@(t,phi) odefcn(t,phi,A_sp), [t0 t_min], phi_t0(:,n_min));
        sk_pr = -(phi(length(t),1)*psi_opt(length(t_opt),1)+phi(length(t),2)*psi_opt(length(t_opt),2));
        [op_func,hhhhh]=Square_Lebesgue(-[phi(length(t),1),phi(length(t),2)]);        
    end
    plot(psi_opt(:,1),psi_opt(:,2),'b','LineWidth',2);
    title('Графики компонент кандидатов на оптимальные траектории x2(x1)')

    xlabel('x1');
    ylabel('x2');
    disp('Задача разрешима');
    disp('Оптимальное время');
    disp(t_min);
    disp('Номер угла');
    disp(n_min);
    disp('Вектор psi(t0) оптимальный');
    disp(phi_t0(:,n_min));
    disp('Погреность из условия трансвенсальности для X1');
    disp(abs(sk_pr-op_func));
else
    disp('Задача не разрешима');
    disp('Время просчета');
    disp(t_max);
end
% Просчет оптимального управления 
axis equal;
hold off
%% Графики компонент оптимального управления
hold on

plot(u_opt(1,:),u_opt(2,:),'b','LineWidth',2);
title('Графики компонент оптимального управления u2(u1)')
drawSet(@Ellipse_Lebesgue,20);
lgd = legend('u2(u1)','NumColumns',2);
xlabel('u opt_1');
ylabel('u opt_2');
axis equal;
hold off
%% Компонент оптимальной траектории

plot(psi_opt(:,1),psi_opt(:,2),'b','LineWidth',2);
title('Графики компонент оптимального траектории x2(x1)')

hold on
plot(x0(1),x0(2),'r.','MarkerSize',30);
drawSet(@Square_Lebesgue,10);

xlabel('x1');
ylabel('x2');
lgd = legend('x^* opt','NumColumns',2);

axis equal;
hold off

%% Графики сопряженных переменных
plot(phi(:,1),phi(:,2),'b','LineWidth',2);
title('Графики сопряженных переменных psi2(psi1)')

hold on

xlabel('psi1');
ylabel('psi2');
lgd = legend('psi^* opt','NumColumns',2);
axis equal;
hold off

%% Графики компонент оптимального управления от времени
hold on
subplot(2,1,1);
plot(t_uopt,u_opt(1,:),'b','LineWidth',2);
title('График компоненты оптимального управления  u1(t)')
legend('u1(t)');
xlabel('t');
ylabel('u1(t)');
subplot(2,1,2);
plot(t_uopt,u_opt(2,:),'b','LineWidth',2);
title('График компоненты оптимального управления u2(t)')
legend('u2(t)');
xlabel('t');
ylabel('u2(t)');
hold off
%% Графики сопряженных переменных от времени
hold on
subplot(2,1,1);
plot(t,phi(:,1),'b','LineWidth',2);
axis([t0 t_min min(phi(:,1)) max(phi(:,1))]);
title('График сопряженной переменной от времени  psi1(t)')
legend('psi1(t)');
xlabel('t');
ylabel('psi1(t)');
subplot(2,1,2);
plot(t,phi(:,2),'b','LineWidth',2);
axis([t0 t_min min(phi(:,2)) max(phi(:,2))]);
title('Графики компонент оптимального управления psi2(t)')
legend('psi2(t)');
xlabel('t');
ylabel('psi2(t)');
hold off
%%

function dphidt = odefcn_OC(t,phi,A,B,f,setka,u)
    u_1 = zeros(1,length(setka));
    u_2 = zeros(1,length(setka));
    u_1(1:length(setka)) = B(1,1)*u(1,1:length(setka))+B(1,2)*u(2,1:length(setka));
    u_2(1:length(setka)) = B(2,1)*u(1,1:length(setka))+B(2,2)*u(2,1:length(setka));

    U_1 = interp1(setka,u_1,t); % Interpolate the data set (ft,f) at time t
    U_2 = interp1(setka,u_2,t); % Interpolate the data set (ft,f) at time t

    dphidt = zeros(2,1);
    dphidt(1) = A{1,1}(t)*phi(1)+A{1,2}(t)*phi(2)+f{1}(t)+U_1;
    dphidt(2) = A{2,1}(t)*phi(1)+A{2,2}(t)*phi(2)+f{2}(t)+U_2;
end
function dphidt = odefcn(t,phi,A)
    dphidt = zeros(2,1);
    dphidt(1) = A{1,1}(t)*phi(1)+A{1,2}(t)*phi(2);
    dphidt(2) = A{2,1}(t)*phi(1)+A{2,2}(t)*phi(2);
end
function [val, point]  = SupportLebesgue_2(f,opts)
    l1=  opts.lx;
    l2 = opts.ly;
    A = [];
    b = [];
    fun = @(x) -(x(1)*l1+x(2)*l2);
    x0 = opts.x0;
    Aeq = [];
    beq = [];
    lb=[];
    ub=[];
    save params_2 f;
    nonlcon = @unitdisk_2;
    [x,fval] = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
    val = -fval;
    point = x;
end
function [c,ceq] = unitdisk_2(x)
    load params_2 f;
    ceq = f(x);
    c = [];
end
function [val, point] = Ellipse_Lebesgue(y)
    load params_P params;
    g = @(x) params(1)*(x(1)-params(3))^2-1+params(2)*(x(2)-params(4))^2;
    s = struct('lx',y(1),'ly',y(2),'x0',[params(3)+0.1,params(4)]);
    [val, point] = SupportLebesgue_2(g,s);
end
function [val, point] = Square_Lebesgue(y)
% \mathcal{X}_1 - целевое множество значений фазового вектора
% \mathcal{X}_1 = {x | max { \alpha |x_1-k|, \beta |x_2-m| } <= q } q>0
% если альфа и бета положительны то прямоугольник, если одна не
% положительна, то полоса; если обе то вся плоскость
%params1 = [k,m,q,alpha,beta];
    load params_X1 params1;
    g = @(x) params1(3)-max(params1(4)*abs(x(1)-params1(1)),params1(5)*abs(x(2)-params1(2)));
    s = struct('lx',y(1),'ly',y(2),'x0',[params1(1),params1(2)]);
    [val, point] = SupportLebesgue_2(g,s);
end
function res = drawSet(rho,N)
    %t = linspace(0,2*pi,400);
    %x_r = cos(t)*2;
    %y_r = sin(t)-2;
    %plot(x_r,y_r);
    hold on
    p = linspace(0,2*pi-2*pi/N,N);
    x_t = cos(p);
    y_t = sin(p);
    [val, point] = rho([x_t(1),y_t(1)]);
    point1 = point;
    point_last = point;

    for i = 2:N
        [val, point] = rho([x_t(i),y_t(i)]);
        % строим внутр прямую 
        alf = linspace(point_last(1),point(1),100);
        bet = linspace(point_last(2),point(2),100);
        plot(alf,bet,'r');
        % строим внешнюю прямую
        c1 = -(point_last(1)*x_t(i-1)+point_last(2)*y_t(i-1));
        c2 = -(point(1)*x_t(i)+point(2)*y_t(i));
        x0 = (c1*y_t(i)-c2*y_t(i-1))/(x_t(i)*y_t(i-1)-x_t(i-1)*y_t(i));
        y0 = (c2*x_t(i-1)-c1*x_t(i))/(x_t(i)*y_t(i-1)-x_t(i-1)*y_t(i));
        
        alf_1 = linspace(point_last(1),x0,100);
        bet_1 = linspace(point_last(2),y0,100);
        
        alf_2 = linspace(x0,point(1),100);
        bet_2 = linspace(y0,point(2),100);
        %plot(alf_1,bet_1,'g',alf_2,bet_2,'g');
        point_last = point;
    end
    % построение последней внутренней прямой
    alf = linspace(point_last(1),point1(1),100);
    bet = linspace(point_last(2),point1(2),100);
    plot(alf,bet,'r');
    % построение последней внешней прямой
    c1 = -(point_last(1)*x_t(N)+point_last(2)*y_t(N));
    c2 = -(point1(1)*x_t(1)+point1(2)*y_t(1));
    x0 = (c1*y_t(1)-c2*y_t(N))/(x_t(1)*y_t(N)-x_t(N)*y_t(1));
    y0 = (c2*x_t(N)-c1*x_t(1))/(x_t(1)*y_t(N)-x_t(N)*y_t(1));

    alf_1 = linspace(point_last(1),x0,100);
    bet_1 = linspace(point_last(2),y0,100);

    alf_2 = linspace(x0,point1(1),100);
    bet_2 = linspace(y0,point1(2),100);
    %plot(alf_1,bet_1,'g',alf_2,bet_2,'g');
end
function [value,isterminal,direction] = StopEvents(t,y)
load params_X1 params1;
%params1 = [k,m,q,alpha,beta];

g = @(x) params1(3)-max(params1(4)*abs(x(1)-params1(1)),params1(5)*abs(x(2)-params1(2)));
value = 0+(g([y(1), y(2)])>=0);     % Detect height = 0
isterminal = 1;   % Stop the integration
direction = 0;   % Negative direction only
end





