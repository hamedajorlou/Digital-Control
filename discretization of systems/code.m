% Hamed Ajorlou 97101167
%% L
A = [0 1;-8/9 -1];
B = [0;10/9];
C = [1 0];
D = 0;
[a, b] = ss2tf(A,B,C,D);
G = tf(a,b);
tau = linspace(0, 250, 1000);
u = ones(1,200);
u1 = [0*u, 0.009 * u, 0* u, -0.009 * u, 0*u];
u2 = [0*u, 0.018 * u, 0* u, -0.018 * u, 0*u];
u3 = [0*u, 0.027 * u, 0* u, -0.027 * u, 0*u];
y1 = lsim(G , u1,tau);
y1 = y1 + 0.1;
y2 = lsim(G,u2,tau);
y2 = y2 + 0.1;
y3 = lsim(G,u3,tau);
y3 = y3 + 0.1;
figure(1)
subplot(2,1,1);
plot(tau , y1,'r', tau , y2,'k', tau , y3,'b');
ylabel('output')
grid on;
legend('10%', '20%', '30%')
subplot(2,1,2); 
u1 = u1 + 0.09;
u2 = u2 + 0.09;
u3 = u3 + 0.09;
plot(tau , u1,'r', tau , u2,'k', tau , u3,'b');
xlabel('t')
ylabel('input')                                         
legend('10%', '20%', '30%')
grid on;
%% NL
u = 0.09*ones(1,200);
tau = linspace(0, 250, 1000); 
u1_NL = [u, 1.1 * u, u, 0.9 * u, u];
u2_NL = [u, 1.2 * u, u, 0.8 * u, u];
u3_NL = [u, 1.3 * u, u, 0.7 * u, u];
ts = [0 250];
y0 = [0.1;0];
[t1,y1_NL] = ode45(@(t,y) foo(t, y, u1_NL,tau), ts, y0);
[t2,y2_NL] = ode45(@(t,y) foo(t, y, u2_NL,tau), ts, y0);
[t3,y3_NL] = ode45(@(t,y) foo(t, y, u3_NL,tau), ts, y0);
figure(2);
subplot(2,1,1)
plot(t1, y1_NL(:,1), 'r', t2, y2_NL(:,1), 'k', t3, y3_NL(:,1), 'b')
grid on;
ylabel('output')
subplot(2,1,2)
plot(tau , u1_NL,'r', tau , u2_NL,'k', tau , u3_NL,'b');
xlabel('t')
ylabel('input')
legend('10%', '20%', '30%')
grid on
%% 3
% results approve what we obtained theoretically
Ts = 0.166;
A_D = eye(2) + Ts * A
B_D = Ts * B
C_D = C
D_D = D
%% 4
% discreting using ZOH 
clc;
Ts = 0.166;
C_system = ss(A, B, C, D);
D_system = c2d(C_system, Ts, 'zoh')
G_D = c2d(G, Ts, 'zoh');

%% 5
clc; 
t = 0:Ts: 250;
pulse1 = pulse(t,0.09, 0.1);
pulse2 = pulse(t,0.09, 0.2);
pulse3 = pulse(t,0.09, 0.3);
x1_i = 0.1; 
x2_i = 0;
y_discrete1 = Numerical(x1_i, x2_i, pulse1, t, Ts);
y_discrete2 = Numerical(x1_i, x2_i, pulse2, t, Ts);
y_discrete3 = Numerical(x1_i, x2_i, pulse3, t, Ts);
figure();
plot(t1, y1_NL(:,1), 'r', t2, y2_NL(:,1), 'k', t3, y3_NL(:,1), 'b')
hold on;
grid on;
xlabel('t');
ylabel('y');
plot(t, y_discrete1, 'r.-', t, y_discrete2, 'k.-', t, y_discrete3, 'b.-')
legend ('10% Continuous', '20% Continuous', '30% Continuous', '10% Discrete', '20% Discrete', '30% Discrete')
%% 6
clc;
pulse_1 = pulse1 - 0.09;
pulse_2 = pulse2 - 0.09;
pulse_3 = pulse3 - 0.09;
y_discrete1 = lsim(D_system, pulse_1, t);
y_discrete2 = lsim(D_system, pulse_2, t);
y_discrete3 = lsim(D_system, pulse_3, t);
y__1 = y1 - 0.1;
y__2 = y2 - 0.1;
y__3 = y3 - 0.1;
figure();
plot(tau , y__1,'r', tau , y__2,'k', tau , y__3,'b');
hold on;
grid on;
xlabel('t');
ylabel('y');
plot(t, y_discrete1, 'r.-', t, y_discrete2, 'k.-', t, y_discrete3, 'b.-')
legend ('10% Continuous', '20% Continuous', '30% Continuous', '10% Discrete', '20% Discrete', '30% Discrete')
%% 7 Nonlinear
for eta = [0.1 0.5 0.9 1 1.1 2 10]
    Ts = 0.166;
    Ts_new = eta * Ts;
    figure();
    plot(t1, y1_NL(:,1),'r', t2, y2_NL(:,1),'k', t3, y3_NL(:,1),'b')
    hold on;
    grid on;
    xlabel('t');
    ylabel('y');
    t = 0:Ts_new:250;
    u = 0.09;
    pulse1 = pulse(t,0.09, 0.1);
    pulse2 = pulse(t,0.09, 0.2);
    pulse3 = pulse(t,0.09, 0.3);
    y_discrete1 = Numerical(x1_i, x2_i, pulse1, t, Ts_new);
    y_discrete2 = Numerical(x1_i, x2_i, pulse2, t, Ts_new);
    y_discrete3 = Numerical(x1_i, x2_i, pulse3, t, Ts_new);
    plot(t, y_discrete1, 'r.-', t, y_discrete2, 'k.-', t, y_discrete3, 'b.-')
    legend ('10% Continuous', '20% Continuous', '30% Continuous', '10% Discrete', '20% Discrete', '30% Discrete')
end
%% 7 Linear
for eta = [0.1 0.5 0.9 1 1.1 2 10]
    Ts = 0.166;
    Ts_new = Ts * eta;
    t = 0:Ts_new: 250;
    pulse1 = pulse(t, u, 0.1);
    pulse2 = pulse(t, u, 0.2);
    pulse3 = pulse(t, u, 0.3);
    pulse_1 = pulse1 - 0.09 ;
    pulse_2 = pulse2 - 0.09;
    pulse_3 = pulse3 - 0.09;
    D_system = c2d(C_system, Ts_new, 'zoh');
    y_discrete1 = lsim(D_system, pulse_1, t);
    y_discrete2 = lsim(D_system, pulse_2, t);
    y_discrete3 = lsim(D_system, pulse_3, t);
    figure;
    plot(tau , y__1,'r', tau , y__2,'k', tau , y__3,'b');
    hold on;
    grid on;
    xlabel('t');
    ylabel('y'); 
    plot(t,y_discrete1, 'r.-', t, y_discrete2, 'k.-', t, y_discrete3, 'b.-')
    legend ('10% Continuous', '20% Continuous', '30% Continuous', '10% Discrete', '20% Discrete', '30% Discrete')
 end
%% 8
clc;
Ts=0.166;
[N,M] = ss2tf(A, B, C, D);
tf_c = tf(N,M)
tf_01 = c2d(tf_c, 0.1 * Ts, 'zoh') 
tf_05 = c2d(tf_c, 0.5 * Ts, 'zoh')
tf_09 = c2d(tf_c, 0.9 * Ts, 'zoh')
tf_1  = c2d(tf_c, 1   * Ts, 'zoh')
tf_11 = c2d(tf_c, 1.1 * Ts, 'zoh')
tf_2  = c2d(tf_c, 2   * Ts, 'zoh')
tf_10 = c2d(tf_c, 10  * Ts, 'zoh') 
%% 9
figure();
bode(tf_c)
hold on;
bode(tf_01);
bode(tf_05);
bode(tf_09);
bode(tf_1);
bode(tf_11);
bode(tf_2);
bode(tf_10); 
legend('Continuous', '0.1T', '0.5T', '0.9T', '1T', '1.1T', '2T', '10T')

%%
function x1_i = Numerical(x1_0, x2_0, u, t, Ts)
x1_i   = zeros(size(t)); 
x2_i = zeros(size(t)); 
x1_i(1)   = x1_0;
x2_i(1) = x2_0;
for n = 1:(length(t) - 1)
    x1_i(n + 1) = x1_i(n) + Ts * x2_i(n) ;
    x2_i(n + 1) = x2_i(n) - Ts * ( x1_i(n) + x2_i(n) ) + Ts * u(n)/(1-x1_i(n));
end
end

%%
function u_out = pulse(t, u, p)
bin1 = find(t < 50 , 1, 'last');
bin2 = find(t < 100, 1, 'last');
bin3 = find(t < 150, 1, 'last');
bin4 = find(t < 200, 1, 'last');
S1 = bin1;
S2 = bin2 - bin1;
S3 = bin3 - bin2;
S4 = bin4 - bin3;
S5 = length(t) - bin4;
u_out1 = u * ones(size(t));
u_out2 = [zeros(1, S1), p * abs(u) * ones(1, S2), zeros(1, S3), -p * abs(u) * ones(1, S4), zeros(1, S5)];
u_out = u_out1 + u_out2;
end
%%
function x_dot = foo(t, y, u, time)
u = interp1(time, u, t);
x_dot = zeros(2,1);
x_dot(1) = y(2);
x_dot(2) = -y(1) - y(2) + u/(1-y(1));
end