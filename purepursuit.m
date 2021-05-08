clear;

%%%%%%%%%%%%%%处理GPCHC文件中的轨迹数据%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = load('GPCHC.mat');
noproave_lon = mean(m.GPCHC(:,13));     %求取经度的平均值

%筛除错误数据%
new_GPCHC = [];
a = 1;
for i = 1:3431

    if abs(m.GPCHC(i,13)-noproave_lon)<0.04
       new_GPCHC(a,1) = m.GPCHC(i,13);
       new_GPCHC(a,2) = m.GPCHC(i,14);
    a = a+1;
    end
    
end
position = new_GPCHC(:,1:2);
%筛除错误数据结束%

%MATLAB程序实现经纬度转换成笛卡尔坐标%：
M_PI=3.14159265358979323846;
L = 6381372 * M_PI * 2;     %地球周长  
W = L;                      % 平面展开后，x轴等于周长  
H = L / 2;                  % y轴约等于周长一半  
mill = 2.3;                 % 米勒投影中的一个常数，范围大约在正负2.3之间  

n=size(position,1);
new_position=[];
for i =1:n
    lon=position(i,1);
    lat=position(i,2);
    x = lon * M_PI / 180; % 将经度从度数转换为弧度  
    y = lat * M_PI / 180; %将纬度从度数转换为弧度  
    y1 = 1.25 * log(tan(0.25 * M_PI + 0.4 * y)); % 米勒投影的转换  
    % 弧度转为实际距离  
    dikaerX = (W / 2) + (W / (2 * M_PI)) * x ; %笛卡尔坐标x
    dikaerY = (H / 2) - (H / (2 * mill)) * y1 ;%笛卡尔坐标y
    new_position(i,1)=dikaerX;
    new_position(i,2)=dikaerY;
    fprintf('第%d个点的',i)
    fprintf('坐标是=(%f %f);',new_position(i,1),new_position(i,2))
    fprintf('\n')
end
%经纬度轨迹完全转化为笛卡尔坐标系%

x_goal = real(new_position(:,1))-23615000;
y_goal = real(new_position(:,2))+6941400;
y_goal = y_goal/10;
%%%%%%%%%%%%%%%%%%%%%%%%%%路径加载完毕%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%载入路径%
cx = x_goal;
cy = y_goal;
%载入参数%
k = 1;        % 前视距离增益
Lfc = 1.0;       % 预见性距离
Kp = 1.0 ;       % 速度比例增益
dt = 0.4  ;      % 时间
L = 0.4  ;       %  轮轴距
target_speed = 2; %目标速度
T = 100000;       %仿真时间
qwe = 1200;
%设置预处理过程%
x_process = [cx(1,1)];
y_process = [cy(1,1)];
i = 1;
N = length(x_process);
x =cx(1,1);
y = cy(1,1); 
yaw = 0;
v = 0;
time = 0;
Lf = k * v + Lfc;
f1 = figure(1)
  
while ((T > time)&&(N~=qwe))
    N = length(x_process);
    target_index = calc_target_current_distance(x,y,cx,cy);
   %计算控制量%
    ai = PIDcontrol(target_speed, v,Kp);
    di = pure_pursuit_control(x,y,yaw,v,cx,cy,target_index,k,Lfc,L,Lf);
    
    [x,y,yaw,v] = updatestate(x,y,yaw,v, ai, di,dt,L)
    time = time + dt;
    x_process = [x_process;x];
    y_process = [y_process;y];
    plot(cx,cy,'b',x,y,'r*')
    title('跟踪过程')
    drawnow
    hold on
end
f2 = figure(2)
plot(x_goal,y_goal,'b',x_process,y_process,'r')
title('跟踪路径')
legend('预定轨迹','实际跟踪路径');

%筛选出距离最近的点%
newx = zeros(qwe+1,1);
newy = zeros(qwe+1,1);
locat = 0;
for p = 1:(qwe+1)
    Dis = zeros(3429,1);
    for q = 1:3429
    Dis(q,1) =  (x_goal(q,1)-x_process(p,1))^2 + (y_goal(q,1)-y_process(p,1))^2;
    end
    [~, locat]= min(Dis);
    newx(p,1) = x_goal(locat,1);
    newy(p,1) = y_goal(locat,1);
end
%处理数据，计算均方差%
sum = 0;
for i = 1:(qwe+1)
    sum = sum + sqrt((newx(i,1)-x_process(i,1))^2+(newy(i,1)*10-y_process(i,1)*10)^2);
end
 ave = sum/(qwe+1);
 fprintf('轨迹跟踪的均方差为%f\n',ave);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%函数定义%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x, y, yaw, v] = updatestate(x, y, yaw, v, a, delta,dt,L) %更新当前差速机器人状态%
    x = x + v * cos(yaw) * dt;
    y = y + v * sin(yaw) * dt;
    yaw = yaw + v / L * tan(delta) * dt;
    v = v + a * dt;
end

function [a] = PIDcontrol(target_v, current_v, Kp)
    a = Kp * (target_v - current_v);
end

function [delta] = pure_pursuit_control(x,y,yaw,v,cx,cy,index,k,Lfc,L,Lf)
    tx = cx(index);
    ty = cy(index);
    
    alpha = atan2((ty-y),(tx-x))-yaw;
    
    Lf = k * v + Lfc;
   delta = atan(2*L * sin(alpha)/Lf)  ;
end

function [index] = calc_target_current_distance(x,y, cx,cy)
  N =  length(cx);  %计算目标横轴点的所有点数%
  Distance = zeros(N,1);  %初始化距离矩阵%
for i = 1:N
  Distance(i) =  sqrt((cx(i)-x)^2 + (cy(i)-y)^2); %计算出每个目标点到当前点的距离%
end
  [~, location]= min(Distance);
  index = location;
  index = index+1 ;
end
