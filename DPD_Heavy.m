%-------------------------------------------------------------------------
%{ 
    Date:   May 2023 
    Author: Alex Vita
    Topic:  Adpative Digital Pre-Distortion
            MSEE Thesis, University of New Haven 

    This is the heavy version of the adaptive DPD process. This contains
    all of the test cases and results that were used throughout my thesis
    manuscript. There is a lot of analysis that may subtract from the
    understanding of the DPD process, if you would like a cleaner version
    please use DPD_Slim. 

    I recommend following along with my manuscript specifically Chapter 5
    as it will describe the purpose and results of the test cases below.

    General guidance for tap size (N) and convergence factor(u)
    Input: 1.0 kHZ          N: 19       u: .005
    Input: 1.0 & 1.5 kHz    N: 27       u: .0026     
    Input: 1.0 MHz          N: 19       u: .0006

    Intial setup is for the first test case, Figure 5.1. 

    Thank you. 
%}
% -------------------------------------------------------------------------

minn = 1;
tapn = 1;

% for f2 = 1300:1:1900

Fs = 44.1 * 10^3;                     %For Single and Two-Tone Tests
% Fs = 20 * 10^6;                         %For 1.0 MHz Bandwidth Tests
dt = 1/Fs;
StopTime = .05;
t = (0:dt:StopTime-dt)';
tLen = length(t); 

% noise = randn(size(t));
% noisefilter = lowpass(noise,1*10^6,Fs, Steepness=.999); %1.0 MHz BW
% for m = .1 : .1 : 6
% m = 1;                                  %Mulitplier for 1.0 MHz BW
% s = m * noisefilter; 

f1 = 1000;                            %1.0 KHz
s =  2.5*sin(2*pi*f1*t);              

% f2 = 1500;                            %1.5 KHz
% s2 = 2.5*sin(2*pi*f2*t);  

% s = s+s2;                             %Two-Tone

u = .005;                              %convergence factor
N = 19;                                 %Filter taps 
                         
y = zeros(tLen,1); 
e = y;                               
d = y;
x = y;
sd = y;
Ts = tLen;

sPA = (4 * (1./(1+exp( -s))-.5));       %Various PA Models 
% sPAerf = erf((sqrt(pi)/2)*s);
% sPAtan = tanh(s);
% sPAarc = (2/pi)*atan((pi/2)*s);
% sPAalg = s./(1+abs(s)); 


% xTable = timetable(seconds(t),sPA);   % Original Signal Analysis 
% [pxx,f] = pspectrum(xTable,Leakage=.99,FrequencyResolution=110);
% [pxx,f] = pspectrum(xTable,Leakage=.5); 
% powersPA = pow2db(pxx); 

% for N = 15: 1 : 30                     
%     minn = 1; 
%     for u = .0002 : .0001 : .0007
        w = zeros(N,1);  

    for n = (N+1):1:Ts

        sumw = 0; 
        for k = 1:1:N
            sumw = sumw + w(k) * s(n-k);
        end
        y(n) = sumw;

        sd(n) = s(n)-y(n);
        x(n) = (4 * (1./(1+exp( -sd(n)))-.5));

%         if n < Ts/2+1
%             x(n) = erf((sqrt(pi)/2)*sd(n));
%         else
%             x(n) = tanh((sqrt(pi)/2)*sd(n));
%             x(n) = (2/pi)*atan((pi/2)*sd(n));
%             x(n) = sd(n)./(1+abs(sd(n))); 
%         end
            
        e(n) = x(n) - y(n);
        for k = 1:1:N
            w(k) = w(k) + 2 * u * e(n) * s(n-k); 
        end  
    end

%     xTable = timetable(seconds(t),x);
%     [pxx,f] = pspectrum(xTable,Leakage=.99,FrequencyResolution=110); 
%     powerx = pow2db(pxx);
%     s1f = (round(f1/13.7));
%     s2f = (round(f2/13.7));
%     s1imf = (round((2*f1-f2)/13.7));
%     s2imf = (round((2*f2-f1)/13.7));
%     powerdiffs1(minn, tapn)= max(powersPA(s1f) - powerx(s1f)); 
%     powerdiffs2(minn, tapn)= max(powersPA(s2f) - powerx(s2f)); 
%     powerdiffs1im(minn, tapn)= max(powersPA(s1imf) - powerx(s1imf)); 
%     powerdiffs2im(minn, tapn)= max(powersPA(s2imf) - powerx(s2imf)); 

%     xTable = timetable(seconds(t),x);
%     [pxx,f] = pspectrum(xTable,Leakage=.5); 
%     powerx = pow2db(pxx);
%     powerdiffACI(minn, tapn)= max(powersPA(450:820) - powerx(450:820)); 
   
%     minn = minn+1;
%     end
%     tapn = tapn+1; 
% end

%{
    Test Plots for all figures in Chapter 5 of thesis manuscript. Each test
    requires some commenting or uncommenting of main code. 
%}

% 1 freq input, Figure 5.1-------------------------------------------------
% figure; 
% % Orignal Input
% xTable = timetable(seconds(t),s);
% [pxx,f] = pspectrum(xTable,Leakage=.9,FrequencyResolution=90);
% pxx = pow2db(pxx); 
% plot((f),pxx,'LineWidth', 1, LineStyle= '--' ,Color="#ffc61e");
% hold on 
% % Output with DPD
% xTable = timetable(seconds(t),x);
% [pxx,f] = pspectrum(xTable,Leakage=.9,FrequencyResolution=90);
% pxx = pow2db(pxx); 
% plot(f,pxx,'LineWidth', 1, LineStyle= '-' ,Color="#009ade");
% hold on
% % Output With DPD
% xTable = timetable(seconds(t),sPA);
% [pxx,f] = pspectrum(xTable,Leakage=.9,FrequencyResolution=90);
% pxx = pow2db(pxx); 
% plot(f,pxx,'LineWidth', 1, LineStyle= '-' ,Color="#ff1f5b");
% hold on 
% grid on;
% xlabel('Frequency (Hz)')
% ylabel('Gain (dB)')
% legend('Original Input', 'Ouput With DPD', 'Output Without DPD');
% xlim([500 3500]);

% 1 freq input, AM-AM, Figure 5.2------------------------------------------
% figure;
% % Output without DPD 
% plot(abs(s(100:tLen)),abs(sPA(100:tLen)), LineWidth=2, ...
%     Color="#ff1f5b",Marker=".",LineStyle="none");
% hold on
% % Output with DPD
% plot(abs(s(100:5:tLen)),abs(x(100:5:tLen)),LineWidth=2,...
%     Color="#009ade",Marker=".",LineStyle="none");
% hold on
% % Predistorted PA Input
% plot(abs(s(100:5:tLen)),abs(sd(100:5:tLen)),LineWidth=2,...
%     Color="#af58ba",Marker=".",LineStyle="none");
% hold on 
% %Reference Line
% xlin = linspace(0,2.5,500);
% ylin = xlin; 
% plot(abs(xlin),abs(ylin), LineWidth=2, LineStyle=":",...
%     Color="#ffc61e");
% xlabel('Pin')
% ylabel('Pout')
% xlim([0 2.5]);
% grid on
% legend('PA output without DPD', 'PA output with DPD',...
%     'PA input with DPD', 'Reference Line (x = y)')


% % 2 freq input, Figure 5.3-------------------------------------------------
% uncomment second frequency components 
% figure; 
% %Original Input
% xTable = timetable(seconds(t),s);
% [pxx,f] = pspectrum(xTable,Leakage=.99,FrequencyResolution=110);
% pxx = pow2db(pxx); 
% plot((f),pxx,'LineWidth', 1, LineStyle= '--' ,Color="#ffc61e");
% hold on 
% %Output with DPD
% xTable = timetable(seconds(t),x);
% [pxx,f] = pspectrum(xTable,Leakage=.99,FrequencyResolution=110);
% pxx = pow2db(pxx); 
% plot(f,pxx,'LineWidth', 1, LineStyle= '-' ,Color="#009ade");
% hold on
% %Output without DPD
% xTable = timetable(seconds(t),sPA);
% [pxx,f] = pspectrum(xTable,Leakage=.99,FrequencyResolution=110);
% pxx = pow2db(pxx); 
% plot(f,pxx,'LineWidth', 1, LineStyle= '-' ,Color="#ff1f5b");
% hold on 
% grid on;
% xlabel('Frequency (Hz)')
% ylabel('Gain (dB)')
% legend('Original Input', 'Ouput With DPD', 'Output Without DPD');
% xlim([0 2500]);

% 2 Freq Input, Varying Taps, Figure 5.4-----------------------------------
% Uncomment for N = 1: 1 : 50 , minn = minn+1, and end 
% Powerdiff s2 s1 stuff and xtabel ****
% figure;
% N = 1 : 1 : 50;
% plot(N,powerdiffs1,LineWidth=2, Color="#ff1f5b")
% hold on 
% plot(N,powerdiffs2,LineWidth=2, Color="#009ade", LineStyle=":")
% hold on 
% plot(N,powerdiffs1im,LineWidth=2, Color="#af58ba")
% hold on 
% plot(N,powerdiffs2im,LineWidth=2, Color="#ffc61e", LineStyle=":")
% xline (27, LineWidth=1);
% legend('f1', 'f2', '2f1 - f2', ...
%     '2f2 - f1')
% ylabel('Gain (dB), without DPD - with DPD')
% xlabel('# of Taps')
% grid on 

% Varying convergence, Figure 5.5------------------------------------------
% Set u = .0001, .0025, .0026 and save e for each run to e0025/26/01
% figure;
% plot(e0001, LineWidth=1, Color="#ffc61e", LineStyle="--")
% hold on 
% plot(e0026,LineWidth=1, Color="#ff1f5b", LineStyle="-")
% hold on 
% plot(e0025,LineWidth=2, Color="#009ade")
% hold on 
% grid on 
% box on 
% xline (3000, LineWidth=3,Color="#ffc61e",LineStyle=":");
% xline (400, LineWidth=3,Color="#009ade",LineStyle=":");
% legend('u = .0001', 'u = .0026', 'u = .0025')
% xlabel('Iteration')
% ylabel('Error Signal e(n)')
% xlim ([0 4000])
% ylim ([-1.5 1.5])

% Two Freq, Varying 2nd, Figure 5.6----------------------------------------
% uncomment For loop at top, and end, and min = minn+1
% figure;
% t = tiledlayout('flow'); 
% nexttile
% f2 = 1300:1:1900;
% plot(f2,powerdiffs1,LineWidth=2, Color="#ff1f5b")
% hold on 
% plot(f2,powerdiffs2,LineWidth=2, Color="#009ade", LineStyle=":")
% hold on 
% legend('f1', 'f2')
% ylim([2.1 2.4])
% grid on
% nexttile
% plot(f2,powerdiffs1im,LineWidth=2, Color="#af58ba")
% hold on 
% plot(f2,powerdiffs2im,LineWidth=2, Color="#ffc61e", LineStyle=":")
% legend('2f1 - f2','2f2 - f1')
% ylim ([2 14]);
% ylabel(t,'Gain (dB), without DPD - with DPD')
% xlabel(t,'Frequency (Hz)')
% grid on 

% 1Mhz Input, Figure 5.7---------------------------------------------------
% Switch signal to 1.0 MHz BW, PAtch function displays Ajacent Channel 
% figure; 
% xTable = timetable(seconds(t),s);
% [pxx,f] = pspectrum(xTable,Leakage=.5);
% pxx = pow2db(pxx); 
% plot((f/10^6),pxx,'LineWidth', 1, LineStyle= '--' ,Color="#ffc61e");
% hold on 
% 
% xTable = timetable(seconds(t),x);
% [pxx,f] = pspectrum(xTable,Leakage=.5);
% pxx = pow2db(pxx); 
% plot(f/10^6,pxx,'LineWidth', 1, LineStyle= '-' ,Color="#009ade");
% hold on
% 
% xTable = timetable(seconds(t),sPA);
% [pxx,f] = pspectrum(xTable,Leakage=.5);
% pxx = pow2db(pxx); 
% plot(f/10^6,pxx,'LineWidth', 1, LineStyle= '-' ,Color="#ff1f5b");
% hold on 
% grid on;
% patch([1.1 1.1 2 2], [-120 -60 -60 -120],'g','facealpha',.1 )
% xlabel('Frequency (MHz)')
% ylabel('Gain (dB)')
% legend('Original Input', 'Ouput With DPD', 'Output Without DPD');
% xlim([0 8])

% % Varying mulitplier for 1Mhz, Figure 5.8--------------------------------
% Uncomment m loop, powerdiff ACI, and minn=minn+1, and end 
% figure; 
% m = .1 : .1 : 6;
% m = m';
% plot(m,powerdiffACI, LineWidth=2)
% ylabel('Gain (dB), without DPD - with DPD')
% xlabel('Input Mulitplier')
% legend('1.1 MHz to 2.0 MHz')
% grid on 

% 1.0 MHz, Error Power, Figure 5.9-----------------------------------------
% Set StopTime = .6, m = .5, uncomment th eif else loop in the main code
% run the code for each x(n) in the else loop and save e for all, 
%name etanh, eatan, ealg. in that order. comment out x(n) before if
%satement
% t = tiledlayout(3,1,TileSpacing="tight");
% 
% nexttile
% plot((etanh.^2), LineWidth=1,Color="#009ade"); 
% ylim([0 0.001])
% xlim([4*10^6 10*10^6])
% grid on 
% legend('tanh(sd)')
% 
% nexttile
% plot((eatan.^2), LineWidth=1,Color="#af58ba"); 
% ylim([0 0.002])
% xlim([4*10^6 10*10^6])
% grid on 
% legend('arctan(sd)')
% 
% nexttile
% plot((ealg.^2), LineWidth=1,Color="#00cd6c"); 
% ylim([0 0.006])
% xlim([4*10^6 10*10^6])
% grid on 
% legend('sd/1+|sd|')
% 
% xlabel(t, 'Iteration')
% ylabel(t,'Error Power')

% Averge Error Power, 1.0MHz BW, Figure 5.10-------------------------------
% run after the previous test, this uses the same naming for e(n)
% perf = sum((etanh(4*10^6:5*10^6,1)).^2);
% perf = perf/(1*10^6);
% 
% ptanh = sum((etanh(7*10^6:8*10^6,1)).^2);
% ptanh = ptanh/(1*10^6);
% 
% patan = sum((eatan(7*10^6:8*10^6,1)).^2);
% patan = patan/(1*10^6);
% 
% palg = sum((ealg(7*10^6:8*10^6,1)).^2);
% palg = palg/(1*10^6);
% 
% figure; 
% str = ["erf(sd)" "tanh(sd)" "arctan(sd)" "sd/1+|sd|"];
% plot(1,perf, LineStyle="none", Marker="o", MarkerFaceColor="#ff1f5b", ...
%     MarkerEdgeColor="#ff1f5b")
% hold on 
% 
% plot(2,ptanh, LineStyle="none", Marker="o",MarkerFaceColor="#009ade", ...
%     MarkerEdgeColor="#009ade")
% hold on 
% 
% plot(3,patan, LineStyle="none", Marker="o",MarkerFaceColor="#af58ba", ...
%     MarkerEdgeColor="#af58ba")
% hold on 
% 
% plot(4,palg, LineStyle="none", Marker="o",MarkerFaceColor="#00cd6c", ...
%     MarkerEdgeColor="#00cd6c")
% hold on 
% set(gca,'xtick',1:4,'xticklabel',str)
% grid on 
% xlabel("PA Model function")
% ylabel("Average Error Power")

%Tap to u Correleation, Figure 5.11----------------------------------------
% Set stop time =.05, m = 1, uncomment N and u for loopsm, minn = minn +1,
% tapn = tapn+1, and the end of loops. 
% N = 15 : 1 :30;
% u = .0002 : .0001 : .0007;
% 
% surf(N,(u),powerdiffACI);
% hold on 
% xlabel('# of Taps');
% ylabel('Convergence Factor')
% zlabel('Gain (db) Reduction')






% AM-AM Plot Artifacts-----------------------------------------------------
% figure;
% plot(abs(s(100:tLen)),abs(sPA(100:tLen)), LineWidth=2, ...
%     Color="#ff1f5b",Marker=".",LineStyle="none");
% hold on
% % 
% % Distorted
% plot(abs(s(100:10:tLen)),abs(x(100:10:tLen)),LineWidth=2,...
%     Color="#009ade",Marker=".",LineStyle="none");
% hold on
% 
% % Distorted priot to amplification
% plot(abs(s(100:10:tLen)),abs(sd(100:10:tLen)),LineWidth=2,...
%     Color="#af58ba",Marker=".",LineStyle="none");
% hold on 
% 
% %Reference Line
% x = linspace(0,2.5,500);
% y = x; 
% plot(abs(x),abs(y), LineWidth=2, LineStyle=":",...
%     Color="#ffc61e");
% 
% xlabel('Pin')
% ylabel('Pout')
% xlim([0 2.5]);
% grid on
% legend('PA output without DPD', 'PA output with DPD',...
%     'PA input with DPD', 'Reference Line (x = y)')

% % AM-AM Different Sigmoid Functions--------------------------------------
% Uncomment differnt SPA variables
% figure;
% % Error Function 
% plot(abs(s(100:tLen)),abs(sPAerf(100:tLen)),LineWidth=2,...
%     Color="#ff1f5b",Marker="none",LineStyle="-");
% hold on 
% 
% % Hyperbolic
% plot(abs(s(100:tLen)),abs(sPAtan(100:tLen)),LineWidth=2,...
%     Color="#009ade",Marker="none",LineStyle="-");
% hold on
% 
% %Arctan
% plot(abs(s(100:tLen)),abs(sPAarc(100:tLen)),LineWidth=2,...
%     Color="#af58ba",Marker="none",LineStyle="-");
% hold on 
% 
% plot(abs(s(100:tLen)),abs(sPAalg(100:tLen)),LineWidth=2,...
%     Color="#00cd6c",Marker="none",LineStyle="-");
% hold on 
% 
% 
% %Reference Line
% x = linspace(0,1.5,500);
% y = x; 
% plot(abs(x),abs(y), LineWidth=2, LineStyle=":",...
%     Color="#ffc61e");
% 
% xlabel('Pin')
% ylabel('Pout')
% xlim([0 1.5]);
% ylim ([0 1.5])
% grid on
% legend('erf(sd)', 'tanh(sd)',...
%     'arctan(sd)','sd/1+|sd|', 'Reference Line (x = y)')

