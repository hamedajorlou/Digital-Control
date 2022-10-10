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
y2 = lsim(G , u2,tau);
y2 = y2 + 0.1;
y3 = lsim(G , u3,tau);
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
u = 0.09*ones(1,2000);
tau = linspace(0, 250, 10000); 
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
%%
function x_dot = foo(t, y, u, time)
u = interp1(time, u, t);
x_dot = zeros(2,1);
x_dot(1) = y(2);
x_dot(2) = -y(1) - y(2) + u/(1-y(1));
end