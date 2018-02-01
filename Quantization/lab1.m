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

%% 5) Quantization (Rounding)

% This script creates a signal, and then quantizes it to a specified number
% of bits.  It then calculates the quantization error.
% see if you run the script.

fprintf('\nE71 Lab, Sampling and Quantization\n');

% Choose the input type.
choice = questdlg('Choose input','Input',...
    'Sine','Sawtooth','Triangle','Sine');

for b=[2,4,6,8,10]

    N=120;                          % Number of samples in final signal.
    n=0:(N-1);                 %Index

    fprintf('Bits = %g, levels = %g, signal = %s.\n', b, 2^b, choice);

    % Create the  input data sequence.
    switch choice
        case 'Sine'
            x=sin(2*pi*n/N);
        case 'Sawtooth'
            x=sawtooth(2*pi*n/N);
        case 'Triangle'
            x=abs(sawtooth(2*pi*n/N));
        case 'Random'
            x=randn(1,N);       % Random data
            x=x/max(abs(x));    % Scale to +/- 1
    end

    % Signal is restricted to between -1 and +1.
    x(x>=1)=(1-eps);            % Make  signal from -1 to just less than 1.
    x(x<-1)=-1;

    % Quantize a signal to "b" bits.  
    xq=floor((x+1)*2^(b-1));    % Signal is one of 2^n int values (0 to 2^n-1)
    xq=xq/(2^(b-1));            % Signal is from 0 to 2 (quantized)
    xq=xq-(2^(b)-1)/2^(b);      % Shift signal down (rounding)

    xe=x-xq;                    % Quantization error
 
    % Calculate SQNR
    SQNR = 10*log10(sum(x.^2)/sum(xe.^2));
    SQNR2 = 1.76+6.02*b;
    fprintf('SQNR = %g, SQNR2 = %g\n', SQNR, SQNR2);
    
    subplot(3,2,b/2);
    stem(x,'b');
    hold on;
    stem(xq,'r');
    hold on;
    stem(xe,'g');
    legend('exact','quantized','error','Location','Southeast')
    title(sprintf('Signal, Quantized signal and Error for %g bits, %g quantization levels',b,2^b));
    hold off

end

%% 6) Quantization (Truncation)

%% 7) Quantization (a la Proakis et Manolakis)

%% 8) Oversampling


