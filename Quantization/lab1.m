%% 1) Delayed Sample Function

N = 8; M=4;
d = [1 zeros(1,N-1)];
dm = circshift(d,M); % Shifts by M

subplot(2,1,1); stem(d, 'LineWidth', 1.5); axis([1 10 0 1.5]); 
title('Unit Sample Function, N=8'); ylabel('Unit Sample')
subplot(2,1,2); stem(dm, 'LineWidth', 1.5); axis([1 10 0 1.5]); 
title('Delayed Unit Sample Function, N=8, M=4'); ylabel('Delayed Unit Sample');

%% 2) More Simple Functions

N = 10;
unitSample = [1 zeros(1,N-1)];
unitStep = ones(1,N);
unitRamp = zeros(1,N);
for x = 1:N
    unitRamp(x) = x-1;
end

hold on
stem(unitSample,'r','LineWidth', 1.5)
stem(unitStep,'LineWidth', 1.5)
stem(unitRamp,'LineWidth', 1.5)
hold off
legend('Sample','Step','Ramp')

%% 3) Even More of the Same

A = 8; L = 30; N = 11; Fs = 10000;
t = 0:1:L-1;
x = A*sawtooth(2*pi*Fs/N*t);
y = A*square(2*pi*Fs/N*t, 40);

subplot(2,1,1)
stem(t,x, 'Linewidth',2); axis([0 30 -10 10]);
subplot(2,1,2)
stem(t,y, 'Linewidth',2); axis([0 30 -10 10]);
grid on;
%% 4) Aliasing

t = 0:1/2000:1-1/2000;
x = sin(2*pi*3*t);
y = sin(2*pi*7*t);
z = sin(2*pi*13*t);
subplot(4,1,1);
hold on
plot(t,x,'--r','Linewidth',2)
plot(t,y,':g','Linewidth',2)
plot(t,z,'-b.','Linewidth',2)

t = 0:1/10:1-1/10;
x = sin(2*pi*3*t);
y = sin(2*pi*7*t);
z = sin(2*pi*13*t);
hold on
subplot(4,1,2);
stem(x,'r')
subplot(4,1,3);
stem(y,'g')
subplot(4,1,4);
stem(z,'b')
