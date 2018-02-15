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
title('Sawtooth : A=8, L=30, N=11, Fs=10kHz')
xlabel('Time in Units of 0.0001s'); ylabel('Amplitude');
subplot(2,1,2)
stem(t,y, 'Linewidth',2); axis([0 30 -10 10]);
title('Square : A=8, L=30, N=11, Fs=10kHz')
xlabel('Time in Units of 0.0001s'); ylabel('Amplitude');
grid on;


%% 4) Aliasing

t = 0:.0001:.5;
tenHz = 0:.1:.5;

threeHz = sin(2*pi*3*t);
sevenHz = sin(2*pi*7*t);
thirteenHz = sin(2*pi*13*t);

hold on
subplot(4,1,1);
plot(t,threeHz,'--r',t,sevenHz,':g',t,thirteenHz,'-b.','Linewidth',2);
legend('3Hz','7Hz','13Hz');
title('3Hz, 7Hz, 13Hz Sine Waves');

threeHz_sampled = sin(2*pi*3*tenHz);
sevenHz_sampled = sin(2*pi*7*tenHz);
thirteenHz_sampled = sin(2*pi*13*tenHz);

hold on
subplot(4,1,2);
stem(tenHz,threeHz_sampled,'r')
title('3Hz Wave Sampled at 10Hz');
subplot(4,1,3);
stem(tenHz,sevenHz_sampled,'g');
title('7Hz Wave Sampled at 10Hz');
subplot(4,1,4);
stem(tenHz,thirteenHz_sampled,'b')
title('13Hz Wave Sampled at 10Hz');

%% 5A) Quantization (Rounding)
% This script creates a signal, and then quantizes it to a specified number
% of bits.  It then calculates the quantization error.
% see if you run the script.

fprintf('\nE71 Lab, Sampling and Quantization\n');

% Choose the input type.
choice = questdlg('Choose input','Input',...
    'Sine','Sawtooth','Random','Random');

for b=[2,4,6,8,10]

    N=120;                     % Number of samples in final signal.
    n=0:(N-1);                 %Index

    fprintf('Bits = %g, levels = %g, signal = %s.\n', b, 2^b, choice);

    % Create the  input data sequence.
    switch choice
        case 'Sine'
            x=sin(2*pi*n/N);
        case 'Sawtooth'
            x=sawtooth(2*pi*n/N);
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
    figure;
    fprintf('SQNR = %g, SQNR2 = %g\n', SQNR, SQNR2);
    hold on;
    stem(x,'b');
    hold on;
    stem(xq,'r');
    hold on;
    stem(xe,'g');
    legend('exact','quantized','error','Location','Southeast')
    title(sprintf('Signal, Quantized signal and Error for %g bits, %g quantization levels',b,2^b));
    hold off

end


%% 5b) Quantization Cont'd. (Rounding)
% This script creates a signal, and then quantizes it to a specified number
% of bits.  It then calculates the quantization error.
% see if you run the script.

fprintf('\nE71 Lab, Sampling and Quantization\n');

% Choose the input type.
choice = questdlg('Choose input','Input',...
    'Sine','Triangle','Random','Random');

for b=[2,4,6,8,10]

    N=120;                          % Number of samples in final signal.
    n=0:(N-1);                 %Index

    fprintf('Bits = %g, levels = %g, signal = %s.\n', b, 2^b, choice);

    % Create the  input data sequence.
    switch choice
        case 'Sine'
            x=sin(2*pi*n/N);
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
    SQNR2 = 6.02*b;
    
    figure;
    fprintf('SQNR = %g, SQNR2 = %g\n', SQNR, SQNR2);
    stem(x,'b');
    hold on;
    stem(xq,'r');
    hold on;
    stem(xe,'g');
    legend('exact','quantized','error','Location','Southeast')
    title(sprintf('Signal, Quantized signal and Error for %g bits, %g quantization levels',b,2^b));
    hold off

end

%% 6A) Quantization (Truncation)
% This script creates a signal, and then quantizes it to a specified number
% of bits.  It then calculates the quantization error.
% see if you run the script.

fprintf('\nE71 Lab, Sampling and Quantization\n');

% Choose the input type.
choice = questdlg('Choose input','Input',...
    'Sine','Triangle','Random','Random');

for b=[2,4,6,8,10]

    N=120;                          % Number of samples in final signal.
    n=0:(N-1);                 %Index

    fprintf('Bits = %g, levels = %g, signal = %s.\n', b, 2^b, choice);

    % Create the  input data sequence.
    switch choice
        case 'Sine'
            x=sin(2*pi*n/N);
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
    xq=floor((x)*2^(b-1));    % Signal is one of 2^n int values (0 to 2^n-1)
    xq=xq/(2^(b-1));            % Signal is from 0 to 2 (quantized)
    %xq=xq-(2^(b)-1)/2^(b);      % Shift signal down (rounding)

    xe=x-xq;                    % Quantization error
 
    % Calculate SQNR
    SQNR = 10*log10(sum(x.^2)/sum(xe.^2));
    SQNR2 = -4.26+6.02*b;
    
    fprintf('SQNR = %g, SQNR2 = %g\n', SQNR, SQNR2);
    stem(x,'b');
    hold on;
    stem(xq,'r');
    hold on;
    stem(xe,'g');
    legend('exact','quantized','error','Location','Southeast')
    title(sprintf('Signal, Quantized signal and Error for %g bits, %g quantization levels',b,2^b));
    hold off

end

%% 6b) Quantization Cont'd. (Truncation)
% This script creates a signal, and then quantizes it to a specified number
% of bits.  It then calculates the quantization error.
% see if you run the script.

fprintf('\nE71 Lab, Sampling and Quantization\n');

% Choose the input type.
choice = questdlg('Choose input','Input',...
    'Sine','Triangle','Random','Random');

for b=[2,4,6,8,10]

    N=120;                     % Number of samples in final signal.
    n=0:(N-1);                 %Index

    fprintf('Bits = %g, levels = %g, signal = %s.\n', b, 2^b, choice);

    % Create the  input data sequence.
    switch choice
        case 'Sine'
            x=sin(2*pi*n/N);
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
    xq=floor((x)*2^(b-1));    % Signal is one of 2^n int values (0 to 2^n-1)
    xq=xq/(2^(b-1));            % Signal is from 0 to 2 (quantized)
    % xq=xq-(2^(b)-1)/2^(b);      % Shift signal down (rounding)
  
    xe=x-xq;                    % Quantization error
 
    % Calculate SQNR
    SQNR = 10*log10(sum(x.^2)/sum(xe.^2));
    SQNR2 = 6.02*b - 6.02;
    
    figure;
    fprintf('SQNR = %g, SQNR2 = %g\n', SQNR, SQNR2);
    stem(x,'b');
    hold on;
    stem(xq,'r');
    hold on;
    stem(xe,'g');
    legend('exact','quantized','error','Location','Southeast')
    title(sprintf('Signal, Quantized signal and Error for %g bits, %g quantization levels',b,2^b));
    hold off

end

%% 7A) Quantization (à la Proakis et Manolakis)
% This script creates a signal, and then quantizes it to a specified number
% of bits.  It then calculates the quantization error.
% see if you run the script.

fprintf('\nE71 Lab, Sampling and Quantization\n');

% Choose the input type.
choice = questdlg('Choose input','Input',...
    'Sine','Triangle','Random','Random');

for b=[2,4,6,8,10]

    N=120;                          % Number of samples in final signal.
    n=0:(N-1);                 %Index

    fprintf('Bits = %g, levels = %g, signal = %s.\n', b, 2^b, choice);

    % Create the  input data sequence.
    switch choice
        case 'Sine'
            x=sin(2*pi*n/N);
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
    % xq=xq-(2^(b)-1)/2^(b);      % Shift signal down (rounding)
    xq=round(x*2^(b-1))/2^(b-1);
    xe=x-xq;                    % Quantization error
 
    % Calculate SQNR
    SQNR = 10*log10(sum(x.^2)/sum(xe.^2));
       
    figure
    fprintf('SQNR = %g\n', SQNR);
    stem(x,'b');
    hold on;
    stem(xq,'r');
    hold on;
    stem(xe,'g');
    legend('exact','quantized','error','Location','Southeast')
    title(sprintf('Signal, Quantized signal and Error for %g bits, %g quantization levels',b,2^b));
    hold off

end

%% 7b)  Quantization (à la Proakis et Manolakis)
% This script creates a signal, and then quantizes it to a specified number
% of bits.  It then calculates the quantization error.
% see if you run the script.

fprintf('\nE71 Lab, Sampling and Quantization\n');

% Choose the input type.
choice = questdlg('Choose input','Input',...
    'Sine','Triangle','Random','Random');

for b=[2,4,6,8,10]

    N=120;                          % Number of samples in final signal.
    n=0:(N-1);                 %Index

    fprintf('Bits = %g, levels = %g, signal = %s.\n', b, 2^b, choice);

    % Create the  input data sequence.
    switch choice
        case 'Sine'
            x=sin(2*pi*n/N);
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
    %xq=xq-(2^(b)-1)/2^(b);      % Shift signal down (rounding)
    xq=round(x*2^(b-1))/2^(b-1);
    
    xe=x-xq;                    % Quantization error
 
    % Calculate SQNR
    SQNR = 10*log10(sum(x.^2)/sum(xe.^2));
       
    figure;
    fprintf('SQNR = %g\n', SQNR);
    stem(x,'b');
    hold on;
    stem(xq,'r');
    hold on;
    stem(xe,'g');
    legend('exact','quantized','error','Location','Southeast')
    title(sprintf('Signal, Quantized signal and Error for %g bits, %g quantization levels',b,2^b));
    hold off

end

%% 8 Oversampling
% This script creates a random signal, and then quantizes it.  The signal
% is oversampled and then decimated.
% * Oversampling is the process of taking in samples at a faster rate (in
% this script "os" times faster) than you need.
% * Decimation is the process Decimation is the process of taking only one
% of every "os" samples of an oversampled signal to get the final sampling
% rate.
%
% Processing of the oversampled signal can give some benefit, as you will
% see if you run the script.

fprintf('\n\nE71 Lab, Oversampling and Quantization\n');

b=3;                            % Number of bits.
N=100;                          % Number of samples in final signal.
maxOS = 4;                      % Max os rate = 2^maxOS.

% Choose the input type.
choice = questdlg('Choose input','Input',...
    'Sine','Sawtooth','Random','Random');

fprintf('Bits = %g, levels = %g, signal = %s.\n', b, 2^b, choice);
fprintf('SQNR measured in dB.\n\n');

%Generate set of points.
    os=2^maxOS;
    N_os=N*os;                   % Number of samples in oversampled signal
    n=0:(N_os-1);                % Index
    
    % Create the oversampled input data sequence.
    switch choice
        case 'Sine'
            xOrig=sin(2*pi*n/N_os);
        case 'Sawtooth'
             xOrig=sawtooth(2*pi*n/N_os);
        case 'Random'
             xOrig=randn(1,N_os);        % Random data
            % Smooth to begin to remove fast variations.
             xOrig=filter(ones(1,4*os)/4/os,1,xOrig);
             xOrig= xOrig/max(abs(xOrig));        % Scale to +/- 1
    end
    
% This large loop generates and analyzes data at several different
% oversampling rates (all powers of two).
for os_pow=0:maxOS,
    os=2^os_pow;                    % Oversampling rate.
    xds = downsample(xOrig,os);
   % x = filter(ones(1,4*os)/4/os,1,xOrig);
    x = filter(ones(1,4)/4,1,xds);
    x = x/max(abs(x));
    
    %Quantize the oversampled raw signal
    xq=floor((x+1)*2^(b-1));	%Signal is one of 2^b int values (0 to 2^b-1)
    xq=xq/(2^(b-1));            %Signal is from 0 to 2 (quantized)
    xq=xq-(2^(b)-1)/2^(b);      %Shift signal down (rounding)
    
    %Smooth (running average) the quantized oversampled signal
    x_qs=filter(ones(1,os)/os,1,xq);
    
    %Smooth the oversampled signal
    x_s=filter(ones(1,os)/os,1,x);
    
    %Quantize the oversampled smoothed signal
    x_sq=floor((x_s+1)*2^(b-1)); % Signal is one of 2^n int values (0 to 2^n-1)
    x_sq=x_sq/(2^(b-1));         % Signal is from 0 to 2 (quantized)
    x_sq=x_sq-(2^(b)-1)/2^(b);   % Shift signal down (rounding)
       
    xe_qs=x_s-x_qs;             % Error between smoothed and quantized/smoothed
    xe_sq=x_s-x_sq;             % Error between smoothed and smoothed/quantized
    

    subplot(3,1,1);
    stem(x,'m');
    hold on;
    stem(x_s,'b');
    hold off;
    set(gca,'YLim',[-1.1,1.1]);
    legend({'Original Signal','Smoothed Signal'});
    title('Original Signal (downsampled) and smoothed signal');
    
    subplot(3,1,2)
    stem(x_s,'b');
    hold on;
    stem(x_sq,'r');
    hold on;
    stem(xe_sq*2^(b-1),'g');
    legend('exact','quantized','error*2^{(b-1)}','Location','Southeast')
    set(gca,'YLim',[-1.1,1.1]);
    title(sprintf('Smoothed then quantized, %g bits, %g levels, os=%g ',b,2^b, os));
    hold off
    
    subplot(3,1,3)
    stem(x_s,'b');
    hold on;
    stem(x_qs,'r');
    hold on;
    stem(xe_qs*2^(b-1),'g');
    legend('exact','quantized','error*2^{(b-1)}','Location','Southeast')
    set(gca,'YLim',[-1.1,1.1]);
    title(sprintf('Quantized then smoothed, %g bits, %g levels, os=%g ',b,2^b, os));
    hold off
    
    sqnr_sq=10*log10(sum(x_s.^2)/sum(xe_sq.^2));
    sqnr_qs=10*log10(sum(x_s.^2)/sum(xe_qs.^2));
    
    fprintf('Oversampling = %g, ', os);
    fprintf('sqnr_sq = %g, sqnr_qs = %g, sqnr_qs-sqnr_sq = %g\n',...
        sqnr_sq, sqnr_qs, sqnr_qs-sqnr_sq);
    
    pause(1);
end

