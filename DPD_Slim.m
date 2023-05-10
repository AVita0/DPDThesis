%-------------------------------------------------------------------------
%{ 
   Date:   May 2023 
   Author: Alex Vita
   Topic:  Adpative Digital Pre-Distortion
           MSEE Thesis, University of New Haven 

    This is the slim version of the adaptive DPD process. There is only one
    test case and data plot. If you would like a more in depth analysis
    with more test cases please use DPD_Heavy.
    
    Thank you. 

    Test:   Ajacent Channel Interference Reduction
    Input:  1.0 MHz bandwidth 

%}
% -------------------------------------------------------------------------


Fs = 20 * 10^6;
dt = 1/Fs;
StopTime = .05;
t = (0:dt:StopTime-dt)';
tLen = length(t); 

% 1.0 MHz Input
noise = randn(size(t));
noisefilter = lowpass(noise,1*10^6,Fs, Steepness=.999);
s = noisefilter; 


u = .0006;                      %Convergence Factor
N = 19;                         %Filter Taps 

w = zeros(N,1);                
y = zeros(tLen,1); 
e = y;                         
d = y;                          
x = y;
sd = y;
Ts = tLen;

sPA = (4 * (1./(1+exp( -s))-.5));

    for n = (N+1):1:Ts

        sumw = 0; 
        for k = 1:1:N
            sumw = sumw + w(k) * s(n-k);            %Adaptive FIR
        end
        y(n) = sumw;

        sd(n) = s(n)-y(n);                          %Pre-Distorter
        x(n) = (4 * (1./(1+exp( -sd(n)))-.5));      %PA Model      
        e(n) = x(n) - y(n);                         %Error Signal 
        
        for k = 1:1:N
            w(k) = w(k) + 2 * u * e(n) * s(n-k);    %LMS Algorithim 
        end
    
    end

% Power Spectrum Analysis--------------------------------------------------
% Patch function highlights ajacent channel 1.1 MHZ - 2.0 MHz

figure; 
% Orignal Signal 
xTable = timetable(seconds(t),s);
[pxx,f] = pspectrum(xTable,Leakage=.5);
pxx = pow2db(pxx); 
plot((f/10^6),pxx,'LineWidth', 1, LineStyle= '--' ,Color="#ffc61e");
hold on 

% Output with DPD 
xTable = timetable(seconds(t),x);
[pxx,f] = pspectrum(xTable,Leakage=.5);
pxx = pow2db(pxx); 
plot(f/10^6,pxx,'LineWidth', 1, LineStyle= '-' ,Color="#009ade");
hold on

% Output Without DPD 
xTable = timetable(seconds(t),sPA);
[pxx,f] = pspectrum(xTable,Leakage=.5);
pxx = pow2db(pxx); 
plot(f/10^6,pxx,'LineWidth', 1, LineStyle= '-' ,Color="#ff1f5b");
hold on 

% patch([1.1 1.1 2 2], [-120 -60 -60 -120],'g','facealpha',.1 )

grid on;
xlabel('Frequency (MHz)')
ylabel('Gain (dB)')
legend('Original Input', 'Ouput With DPD', 'Output Without DPD');
xlim([0 8])
