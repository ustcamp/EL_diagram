%% 绘制Vb和phi_e(或phi_h)对发光强度的EL相图所需主要参数

clc
% clear all
% 温度
T0 = 1e-2; % 温度, 单位: K
% 能级
Vb= 1.6; % 偏压
% Eh = 1.0; % E_F - E_HOMO
% Ee = 1.3; % E_LUMO - E_F
ES1 = 1.8; % S1能量
ET1 = 1.2; % T1能量
% 速率
% FP = 1e0; % 纳腔Purcell因子
K_S1 = 1e12; % 纳腔中S1→S0的退激发速率
K_T1 = 1e7; % 纳腔中T1→S0的退激发速率
K_ISC = 0; % ISC导致的S1→T1的系间穿越速率
K_tip_IES = 1e10; % 分子<->针尖隧穿速率. 注意, 此值仅在IES区域使用, CI区域, K_tip由电流大小I0确定.
K_sub = 1e9; % 分子<->衬底隧穿速率, 100pA~6.25e8,
% 电流
Ie_target = 30; % DC 电流, 用来限制IES区间电流, pA, 1 pA = 160 ns
pA = 6.25e6;
% IES激发效率
eta_IES = 0.01; % 直接(非自旋交换)IES激发效率
eta_IES_ex = 0.1*eta_IES; % 自旋交换IES激发效率
c_weak_IES = 1e-2; % % 如果发生通过分子轨道的CI, 则认为针尖<->衬底的直接隧穿电流占总电流比例很小
if_IES_deexcite = 0; % 是否考虑IES产生的分子退激发. 0: 不考虑; 1: 考虑.
% 是否保持恒流, 此选项决定是否做迭代计算调整K_tip
if_CC=1;

%% %% EL diagram
% 空穴注入势垒phi_h=EF-E_HOMO, 主要影响负偏压发光特性. 只有满足phi_h>E_T1/E_S1, 才能在负偏压|Vb|>Eh时直接通过CI激发分子到T1/S1.
% 电子注入势垒phi_e=E_LUMO-EF, 主要影响正偏压发光特性. 只有满足phi_e>E_T1/E_S1, 才能在正偏压|Vb|>Ee时直接通过CI激发分子到T1/S1.
% 这里主要考察正偏压发光特性, 所以固定phi_h, 同时改变Vb和phi_e.
%%%%

tic

% 计算相图所用参数

dVb = 0.005; % 绘制相图所用偏压范围和步长
Vbmin = 1.0;
Vbmax = 3.0;
% Vb_list = [-Vbmax:dVb:-Vbmin];
% Vb_list = [Vbmin:dVb:Vbmax];
Vb_list = [-Vbmax:dVb:-Vbmin, Vbmin:dVb:Vbmax];
%
phi_h =1.0; % 空穴注入势垒
% phi_e_list = ES1-phi_h:0.01:3.0; % 绘制相图所用的phi_e范围, 注意存在条件E_HL=phi_h+phi_e=ES1+E_Coul>ES1
phi_e_list = 0.8:dVb:Vbmax; % 绘制相图所用的phi_e范围, 注意存在条件E_HL=phi_h+phi_e=ES1+E_Coul>ES1

% 固定HOMO-LUMO gap
% E_HL = 2.8;
% phi_e_list = 0.4:dVb:E_HL-0.4;

% 初始化临时存储变量

N1 = length(Vb_list);
N2 = length(phi_e_list);
pop0 = zeros(5,N2,N1); % 布居数, pop = [PS0, PT1, PS1, PD0+, PD0-]
iCI0 = zeros(N2,N1);
K_tip0 = zeros(N2,N1);
K_sub0 = zeros(N2,N1);
Ie= zeros(N2,N1); % 电流

% 计算分子布居数

tic
parfor k2 = 1:N2
    for k1 = 1:N1
%         phi_h = E_HL-phi_e_list(k2); % 固定HOMO-LUMO gap, 否则注释掉
        [pop0(:,k2,k1), Ie(k2,k1), iCI0(k2,k1), K_tip0(k2,k1), K_sub0(k2,k1)] = ...
            fun1(Vb_list(k1),Ie_target, phi_e_list(k2),phi_h,ET1,ES1, K_T1,K_S1,K_ISC, K_tip_IES, K_sub, eta_IES, eta_IES_ex, c_weak_IES, if_IES_deexcite, T0);
    end
end
toc

% 绘制二维相图

figure(1)
set(gcf,'color','w','defaultaxesfontsize',20)

% pop = [PS0, PT1, PS1, PD0+, PD0-]
eta_el_to_S1 = abs(squeeze(pop0(3,:,:))*K_S1./(Ie*pA)); % PS1发光强度/电流
z0 = eta_el_to_S1;

% 绘制
surf(Vb_list, phi_e_list, log10(z0))
shading flat; view(2); axis tight equal; grid off; colormap parula(256); box on;

%设置绘图颜色范围
z1min = min(log10(z0(:))); z1max = max(log10(z0(:)));
% z1min = -4.15;
z1min = -10.37;
z1max=-0.57;
% z1max=0;
% 注释掉此部分则绘制全部数据/颜色范围
% caxis([z1min, z1max])
% zlim([z1min, z1max])

% 设置坐标轴标识
xlabel('Bias voltage $V_b$ (V)','interpreter','latex')
ylabel('Hole injection barrier $\phi_h$ (eV)','interpreter','latex')
ylabel('Electron injection barrier $\phi_e$ (eV)','interpreter','latex')
xlim([min(Vb_list),max(Vb_list)])
% set(gca,'Xdir','reverse')

% 设置绘图窗口大小, 坐标轴字体大小
s0=2; % s0可改变绘图窗口大小
width0 = s0*15; % 图形宽度, 请勿修改!
height0 = s0*10.5; % 图形高度
axis_font_size = 22; % 一般请勿修改
set(gcf, 'color', 'w', 'Units', 'centimeters', 'position', [-35,5,width0,height0]) % 图形尺寸
set(gca, 'FontName', 'Helvetica', 'FontSize', axis_font_size, 'color', 'none', 'xcolor','k','ycolor','k') % 坐标轴字体字号
ax=gca; ax.XTick = -5:0.4:5; ax.YTick = ax.XTick; ax.LineWidth = 0.8; ax.TickDir = 'out'; ax.TickLength = [0.010 0.02]; ax.Box = 'on'; axis equal;
xlim([min(Vb_list), max(Vb_list)])
% axis off
colorbar


%% 函数

function [pop0, Ie, iCI0, K_tip0, K_sub0] = ...
    fun1(Vb,Ie_target, Ee,Eh,ET1,ES1, K_T1,K_S1,K_ISC, K_tip_IES, K_sub, eta_IES, eta_IES_ex, c_weak_IES, if_IES_deexcite, T0)

%%%%%%%%
% 一些输入参数的说明
%%%%%%%%

% if_IES_deexcite, 是否考虑IES会退激发分子
% c_weak_IES, 分子轨道进入bias window后, 针尖-衬底直接隧穿电流占总电流的比例

% 是否考虑自旋多重度


j1 = 1;

%%%%%%%%
% 一些常规定义
%%%%%%%%

pA = 6.25e6; % 此因子可以把电流从pA转换为速率, 160ns~1pA, 1ns~160pA
I_elec = Ie_target*pA; % 把电流从pA转换为速率
T0 = T0*0.0862e-3; % 温度, K -> eV
f1 = @(x) 1./(exp(x/T0)+1.0); % Fermi distribution
heaviside0 = @(x) (x>=0)*1+(x<0)*0; % 阶跃函数

%%%%%%%%
% 判断通过HOMO或LUMO能否发生CI, 决定发生IES还是CI
%%%%%%%%

if abs(f1(-Eh) - f1(-Eh-Vb)) > 0.01 || abs(f1(Ee) - f1(Ee-Vb)) > 0.01 % 此式的含义是, 分子HOMO或LUMO在两个电极EF之间, 则发生CI
    iCI = 1;
    K_tip = (3*K_sub*I_elec)/(3*K_sub-I_elec); % CI区域, 通过目标电流估计的分子-针尖隧穿速率
%     K_tip = (1*K_sub*I_elec)/(1*K_sub-I_elec);
%     K_tip = K_tip_IES;
else % IES
    iCI = 0;
    K_tip = K_tip_IES;
end
% 分子能级刚好与EF对齐, 此时Fermi分布布居数只有0.5
if abs(abs(f1(-Eh) - f1(-Eh-Vb)) - 0.5) < 0.001 || abs(abs(f1(Ee) - f1(Ee-Vb))-0.5) <0.001
    iCI = 0.5;
end

% if f1(Vb-Ee) < 0.1
%     K_tip = K_tip_CI_pos;
% elseif f1(-Eh-Vb) < 0.1
%     K_tip = K_tip_CI_neg;
% else 
%     K_tip = K_tip_IES;
% end    

%%%%%%%%
% 跃迁速率, IES
%%%%%%%%

% 如果考虑自旋多重度, 则只需要注意S0→T1和S1→T1前面乘上系数j4=3
j4 = 3;

k_IES = eta_IES*I_elec; % 针尖<->衬底的直接隧穿电流激发非自旋翻转过程的概率
k_IES_ex = eta_IES_ex*I_elec; % 针尖<->衬底的直接隧穿电流激发自旋翻转过程的概率

if iCI == 0 % IES区域, 只有当偏压>=激发阈值, 才能通过IES激发分子
    
    k_S0_T1 = j4 * heaviside0(abs(Vb)-ET1) * k_IES_ex;
    k_T1_S1 = heaviside0(abs(Vb)-(ES1-ET1)) * k_IES_ex;
    k_S0_S1 = heaviside0(abs(Vb)-ES1) * k_IES;
    
    % 除了分子本身退激发过程, IES电子也可以退激发分子, 与IES激发分子不同, IES退激发对偏压无要求
    k_T1_S0 = K_T1 + if_IES_deexcite * k_IES_ex; % T1->S0退激发, 包括T1本身退激发以及IES引起的退激发(可选)
    k_S1_T1 = K_ISC+ j4 * if_IES_deexcite * k_IES_ex; % S1->T1退激发, 包括ISC以及IES引起的退激发(可选)
    k_S1_S0 = K_S1 + if_IES_deexcite * k_IES; % S1->S1退激发, 包括S1本身退激发以及IES引起的退激发(可选)
    
else % CI区域, 假设CI区域针尖<->衬底的直接隧穿电流比例很低 (由c_weak_IES决定), 也考虑激发阈值
    
    k_S0_T1 = j4 * c_weak_IES * heaviside0(abs(Vb)-ET1) * k_IES_ex;
    k_T1_S1 = c_weak_IES * heaviside0(abs(Vb)-(ES1-ET1)) * k_IES_ex;
    k_S0_S1 = c_weak_IES * heaviside0(abs(Vb)-ES1) * k_IES;
    
    % IES电子也可以退激发分子, 与IES激发分子不同, IES退激发对偏压无要求
    k_T1_S0 = K_T1 + if_IES_deexcite * c_weak_IES * k_IES_ex;
    k_S1_T1 = K_ISC + j4 * if_IES_deexcite * c_weak_IES * k_IES_ex;
    k_S1_S0 = K_S1 + if_IES_deexcite * c_weak_IES * k_IES;
    
end

%%%%%%%%
% 跃迁速率, CI
%%%%%%%%

% 考虑自旋多重度
% S0→Dp/Dn, j1 = 2
% Dp/Dn→T1, j2 =1.5
% Dp/Dn→S1, j3 = 0.5
j1 = 2; j2 = 1.5; j3 = 0.5;
% 不考虑自旋多重度, 可以设置: j1 = 1; j2 = 1; j3 = 1;

% D0p -> S0,T1,S1

k_D0p_S0 = K_tip*f1(-Eh-Vb) + K_sub*f1(-Eh); % D0p -> S0
k_D0p_T1 = j2 *( K_tip*f1(-Eh-Vb+ET1) + K_sub*f1(-Eh+ET1) ); % D0p -> T1
k_D0p_S1 = j3 *( K_tip*f1(-Eh-Vb+ES1) + K_sub*f1(-Eh+ES1) ); % D0p -> S1

% S0,T1,S1 -> D0p

k_S0_D0p = j1 * ( K_tip*(1-f1(-Eh-Vb)) + K_sub*(1-f1(-Eh)) ); % S0 -> D0p
k_T1_D0p = K_tip*(1-f1(-Eh-Vb+ET1)) + K_sub*(1-f1(-Eh+ET1)); % T1 -> D0p
k_S1_D0p = K_tip*(1-f1(-Eh-Vb+ES1)) + K_sub*(1-f1(-Eh+ES1)); % S1 -> D0p

% D0n -> S0,T1,S1

k_D0n_S0 = K_tip*(1-f1(Ee-Vb)) + K_sub*(1-f1(Ee)); % D0n -> S0
k_D0n_T1 = j2 *( K_tip*(1-f1(Ee-Vb-ET1)) + K_sub*(1-f1(Ee-ET1)) ); % D0n -> T1
k_D0n_S1 = j3 *( K_tip*(1-f1(Ee-Vb-ES1)) + K_sub*(1-f1(Ee-ES1)) ); % D0n -> S1

% S0,T1,S1 -> D0n

k_S0_D0n = j1 * ( K_tip*f1(Ee-Vb) + K_sub*f1(Ee) ); % S0 -> D0n
k_T1_D0n = K_tip*f1(Ee-Vb-ET1) + K_sub*f1(Ee-ET1); % T1 -> D0n
k_S1_D0n = K_tip*f1(Ee-Vb-ES1) + K_sub*f1(Ee-ES1); % S1 -> D0n

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%
% 速率方程
%%%%

%     syms PS0 PT1 PS1 PD0p PD0n PD1p PD1n
%
%     % d/dt PS0
%     eqn1 = -PS0*(k_S0_D0p+k_S0_D0n+k_S0_T1+k_S0_S1) + PT1*k_T1_S0 + PS1*k_S1_S0 + PD0p*k_D0p_S0 + PD0n*k_D0n_S0  == 0;
%     % d/dt PT1
%     eqn2 = PS0*k_S0_T1 - PT1*(k_T1_S0+k_T1_S1+k_T1_D0p+k_T1_D0n) + PS1*k_S1_T1 + PD0p*k_D0p_T1 + PD0n*k_D0n_T1 == 0;
%     % d/dt PS1
%     eqn3 = PS0*k_S0_S1 +PT1*k_T1_S1 - PS1*(k_S1_S0+k_S1_T1+k_S1_D0p+k_S1_D0n) + PD0p*k_D0p_S1 + PD0n*k_D0n_S1  == 0;
%
%     % d/dt PD0p
%     eqn1 = PS0*k_S0_D0p + PT1*k_T1_D0p + PS1*k_S1_D0p - PD0p*(k_D0p_S0+k_D0p_T1+k_D0p_S1) + PD1p*k_D1p_D0p  == 0; % d/dt PD0p
%     % d/dt PD0n
%     eqn2 = PS0*k_S0_D0n + PT1*k_T1_D0n + PS1*k_S1_D0n - PD0n*(k_D0n_S0+k_D0n_T1+k_D0n_S1)   == 0;
%
%     % sum
%     eqn0 = PS0 + PT1 + PS1 + PD0p + PD0n == 1;
%     eqn7 = eqn0;

% A*X= B, 需要给出系数矩阵, 然后用X=linsolve(A,B)求解, X即布居数
% 注意未知数的顺序是: X = pop = [PS0, PT1, PS1, PD0+, PD0-]
% 因为加上布居数和=1条件后, 方程超定, 故去掉最后一个关于d/dt(PD1-)=0的方程, 代换为布居数和=1

A = ...
    [-(k_S0_D0p+k_S0_D0n+k_S0_T1+k_S0_S1), k_T1_S0, k_S1_S0, k_D0p_S0, k_D0n_S0; % PS0
    k_S0_T1, -(k_T1_S0+k_T1_S1+k_T1_D0p+k_T1_D0n), k_S1_T1, k_D0p_T1, k_D0n_T1 % PT1
    k_S0_S1, k_T1_S1, -(k_S1_S0+k_S1_T1+k_S1_D0p+k_S1_D0n), k_D0p_S1, k_D0n_S1; % PS1
    k_S0_D0p, k_T1_D0p, k_S1_D0p, -(k_D0p_S0+k_D0p_T1+k_D0p_S1), 0; % PD0p
    %     k_S0_D0n, k_T1_D0n, k_S1_D0n, 0, -(k_D0n_S0+k_D0n_T1+k_D0n_S1); % PD0n, 此方程不用
    1,1,1,1,1; % sum_Pn = 1
    ];
B = [0
    0
    0
    0
    1];

% 布居数, pop = [PS0, PT1, PS1, PD0+, PD0-]
pop0 = linsolve(A,B);
[ps0,pt1,ps1,pd0p,pd0n]=deal(pop0(1),pop0(2),pop0(3),pop0(4),pop0(5));

%%%%%%%%
% 计算CI过程电流
%%%%%%%%

Ie0_mol_to_tip = ps0*K_tip*((1-f1(-Eh-Vb))) ... % S0→D0+
    + pt1*K_tip*((1-f1(-Eh-Vb+ET1))) ... % T1→D0
    + ps1*K_tip*((1-f1(-Eh-Vb+ES1))) ... % S1→D0+
    + pd0n*K_tip*((1-f1(Ee-Vb)) + (1-f1(Ee-Vb-ET1)) +(1-f1(Ee-Vb-ES1))); % D0-→S0/T1/S1

Ie0_tip_to_mol = ps0*K_tip*(f1(Ee-Vb)) ... % S0→D0-
    + pt1*K_tip*(f1(Ee-Vb-ET1)) ... % T1→D0-
    + ps1*K_tip*((1-f1(-Eh-Vb+ES1))) ... % S1→D0-
    + pd0p*K_tip*(f1(-Eh-Vb) + f1(-Eh-Vb+ET1) + f1(-Eh-Vb+ES1)); % D0+→S0/T1/S1

Ie0 = Ie0_tip_to_mol - Ie0_mol_to_tip;

if iCI == 0 % IES, 电流恒定, 为给定值
    Ie1 = abs(Ie_target);
elseif iCI == 0.5
    Ie1 = 2*abs(Ie0)/pA;
else % CI, 电流通过分子与电极之间的跃迁速率计算出来
    Ie1 = abs(Ie0)/pA;
end
Ie = Ie1;

%%%%%%%%

K_tip0 = K_tip;
K_sub0 = K_sub;
iCI0 = iCI;

end