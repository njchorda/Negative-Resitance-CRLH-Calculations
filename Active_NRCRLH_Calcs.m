%% General Setup
clear;clc;close all
format long

f = linspace(0, 6, 1001)*1e9; % Frequency range
w = 2*pi*f;
c0 = 3e8;

%% Input parameters
Z0 = 50;

% Desired frequency points for phase shifts and number of sub-units
f1 = 2e9;
f2 = 3.5e9;
n = 3;

w1 = 2*pi*f1;
w2 = 2*pi*f2;

% Desired phase shifts at w1 and w2
theta1 = pi/2;% the positive number
theta2 = -pi/2;% the negative number

%% Calculate initial CRLH parameters
LR = Z0*(theta1*(w1/w2) - theta2)/(n*w2*(1 - (w1/w2)^2));
CR = (theta1*(w1/w2) - theta2)/(n*w2*Z0*(1 - (w1/w2)^2));
LL = n*Z0*(1 - (w1/w2)^2)/(w1*(theta1 - (w1/w2)*theta2));
CL = n*(1 - (w1/w2)^2)/(w1*Z0*(theta1 - (w1/w2)*theta2));

disp(['LR = ' num2str(LR/1e-9, 4), ' nH'])
disp(['CR = ' num2str(CR/1e-12, 4), ' pF'])
disp(['LL = ' num2str(LL/1e-9, 4), ' nH'])
disp(['CL = ' num2str(CL/1e-12, 4), ' pF'])

Rse = 0*0.625; % Rse ~0.625 from observation of S21 in simulation
Gsh = Rse/Z0^2; % Assume loss-balanced condition (not usually the case)

%% Calculate desired maximum shunt conductance and Rn trend
G_max = 1/(2*w2*LL);
G_desired = -0.8*G_max;

%% Solve for required NR trend with frequency
% Solve eq. (5) at the two frequencies of interest
Rn_req_f1 = (1 + [-1 1]*sqrt(1 - 4*G_desired^2*(w1*LL)^2))/(2*G_desired);
Rn_req_f2 = (1 + [-1 1]*sqrt(1 - 4*G_desired^2*(w2*LL)^2))/(2*G_desired);

Zn_soln1 = (1 + [-1]*sqrt(1 - 4*G_desired.^2*(w.*LL).^2))./(2*G_desired);
Zn_soln2 = (1 + [1]*sqrt(1 - 4*G_desired.^2*(w.*LL).^2))./(2*G_desired);

Rn_req_f1 = max(Rn_req_f1); % Get the value of required Rn with minimum absolute value
Rn_req_f2 = max(Rn_req_f2);

%% Model Rn as -A*exp(-a*w)
a_rn = log(Rn_req_f1/Rn_req_f2)/(w1 - w2);
af_rn = 2*pi*a_rn;
A_rn = -Rn_req_f1/exp(a_rn*w1);

Rn = -A_rn*exp(af_rn*f);

disp(['A = ' num2str(A_rn)])
disp(['alpha = ' num2str(a_rn*1e9) ' (*1e9)'])

if abs(G_desired) == 0
    Rn = zeros(size(f));
end

%% Calculate series and shunt impedances
Zse = 1i.*w.*LR + 1./(1i.*w.*CL) + Rse; % Zse always passive in this case
Ysh = 1./(1i*w*LL + Rn) + 1i*w*CR + Gsh; % Ysh includes effects from Rn
% Ysh = 1./(1i*w*LL) + 1i*w*CR + G_desired; % Perfectly flat G_sh for test
% Ysh = 1./(1i*w*LL + real(Zn_soln1)) + 1i*w*CR + Gsh;
Ysh_passive = 1./(1i*w*LL) + 1i*w*CR + Gsh; % Passive Y_sh for other cells

%% Calculate ABCD parameters (symmetric, active and passive)
A_sym = 1 + Zse.*Ysh/2;
B_sym = Zse.*(1 + Zse.*Ysh/4);
C_sym = Ysh;
D_sym = 1 + Zse.*Ysh/2;

A_sym_passive = 1 + Zse.*Ysh_passive/2;
B_sym_passive = Zse.*(1 + Zse.*Ysh_passive/4);
C_sym_passive = Ysh_passive;
D_sym_passive = 1 + Zse.*Ysh_passive/2;

%% Cascade ABCD parameters for a given value of 'n'
if mod(n, 2) == 0 % For passive-active-passive symmetry
    error('n must be an odd number');
end
At = zeros(size(f));Bt = zeros(size(f));Ct = zeros(size(f));Dt = zeros(size(f)); %Initialize ABCD vectors
for i = 1:length(f)
    if n == 1 %If n = 1, only the active part matters, no passive
        nPassive = 0; %Multiplies the active part by the identity matrix
    else
        nPassive = (n - 1)/2;
    end
    % For passive-active-passive case
    ABCDtemp = [A_sym_passive(i) B_sym_passive(i);C_sym_passive(i) D_sym_passive(i)]^nPassive*[A_sym(i) B_sym(i);C_sym(i) D_sym(i)]*[A_sym_passive(i) B_sym_passive(i);C_sym_passive(i) D_sym_passive(i)]^nPassive;
    % For active-passive case (one active then n-1 passive)
    % ABCDtemp = [A_sym(i) B_sym(i);C_sym(i) D_sym(i)]*[A_sym_passive(i) B_sym_passive(i);C_sym_passive(i) D_sym_passive(i)]^(n-1);
    % For the all-active case (does not have open stopband issue solved, will blow up near transition frequency)
    % ABCDtemp = [A_sym(i) B_sym(i);C_sym(i) D_sym(i)]^n;
    At(i) = ABCDtemp(1, 1);
    Bt(i) = ABCDtemp(1, 2);
    Ct(i) = ABCDtemp(2, 1);
    Dt(i) = ABCDtemp(2, 2);
end

%% Overall cascaded S-parameters
[S11, S12, S21, S22] = SPARAMS.abcd2s(At, Bt, Ct, Dt, Z0);

%% Plot S-parameters
fig1 = figure(1);
% yyaxis left
plot(f/1e9, 20*log10(abs(S11)));
ylim([-20 5])
hold on
% yyaxis right
plot(f/1e9, 20*log10(abs(S21)));
xline([w1 w2]/(2*pi*1e9))
% ylim([-2 4])
xlim([min(f) max(f)]/1e9)
ylim([-40 5])
grid on
xlabel('Freq (GHz)')
ylabel('S-Params (dB)')
fig1.Position = [734 564 560 420];
title('Cascaded S-Params')
legend({'S_{11}', 'S_{21}'}, 'Location','best')

disp(['S21 at w1 = ' num2str(20*log10(abs(S21(closestIdx(f, f1))))) ' dB'])
disp(['S21 at w2 = ' num2str(20*log10(abs(S21(closestIdx(f, f2))))) ' dB'])

%% Calculate and plot dispersion
Betap = (acosd((At + Dt)/2));

fig2 = figure(2);
plot(f/1e9, real(Betap));
ylim([0 180])
xlim([min(f) max(f)]/1e9)
xline([w1 w2]/(2*pi*1e9))
grid on
xlabel('Freq (GHz)')
ylabel('\beta p (deg)')
title('Dispersion')
fig2.Position = [1295 564 560 420];

%% Plot shunt conductance
fig3 = figure(3);
plot(f/1e9, real(Ysh));
xline([w1 w2]/(2*pi*1e9))
yline(G_desired)
ylim([-0.01 0])
xlim([min(f) max(f)]/1e9)
xlabel('Freq (GHz)')
ylabel('Shunt G (S)')
title('Shunt Conductance vs. Frequency')
fig3.Position = [1296 58 560 420];

%% Plot required Rn
fig4 = figure(4);
% Calculate explicit solution across frequency eq. (5)
Rn_exact_1 = (1 + sqrt(1 - 4.*G_desired.^2.*(w.*LL).^2))./(2.*G_desired);
Rn_exact_2 = (1 - sqrt(1 - 4.*G_desired.^2.*(w.*LL).^2))./(2.*G_desired);
plot(f/1e9, Rn)
hold on
plot(f/1e9, real(Rn_exact_1))
plot(f/1e9, real(Rn_exact_2))
ylim([-100 0])
xlim([min(f) max(f)]/1e9)
xline([w1 w2]/(2*pi*1e9))
grid on
xlabel('Freq (GHz)')
ylabel('R_n(f) (\Omega)')
legend({'Exponential Rn', 'Exact Rn soln. 1', 'Exact Rn soln. 2'}, 'Location','best')
title('R_n(\omega) vs. Frequency')
fig4.Position = [734 58 560 420];

%% S-parameters of a single active unit cell
[S11_single, S12_single, S21_single, S22_single] = SPARAMS.abcd2s(A_sym, B_sym, C_sym, D_sym, Z0);
fig5 = figure(5);
plot(f/1e9, 20*log10(abs(S11_single)))
hold on
plot(f/1e9, 20*log10(abs(S21_single)))

Leq = (Rn.^2 + (w.*LL).^2)./(w.^2*LL); % Eq. (2)
A_loss = Rse.*real(Ysh).*(w.*sqrt(Leq*CL)).^2;
B_loss = w.*(CL*Rse + Leq.*Gsh);
wL = 1./sqrt(LL*CL); Zc = sqrt(LL/CL);
alpha_L = 0.5*(Rse./Zc + real(Ysh).*Zc);
% alpha_L = (wL./2./w).*(CL*Rse + LL*real(Ysh));
% IL_approx = 20*log10(abs(exp(-alpha_L)));
IL_approx = -8.686*((Rse./Zc + real(Ysh).*Zc))/2;%-8.686*alpha_L;
disp(['Approx S21 at w1 = ' num2str(IL_approx(closestIdx(f, f1))) ' dB'])
disp(['Approx S21 at w2 = ' num2str(IL_approx(closestIdx(f, f2))) ' dB'])
plot(f/1e9, IL_approx)
% plot(f/1e9, IL_approx1)
ylim([-30 2])
grid on
xlabel('Freq (GHz)')
ylabel('S-Params (dB)')
title('Unit Cell S-Params')
legend({'S_{11}', 'S_{21}', 'S_{21approx}'}, 'Location','best')
fig5.Position = [171   564   560   420];

%% Bloch Impedance
fig6 = figure(6);
Z_B = Bt./sqrt(At.^2 - 1); % From full cascaded
% Z_B = B_sym./sqrt(A_sym.^2 - 1); % From only active
% Z_B = B_sym_passive./sqrt(A_sym_passive.^2 - 1); % From only passive
plot(f/1e9, abs(real(Z_B)))
hold on
plot(f/1e9, abs(imag(Z_B)))
ylim([0 100])
grid on
xlabel('Freq (GHz)')
title('Bloch Impedance')
legend({'re(Z_B)', 'im(Z_B)'}, 'Location','best')
fig6.Position = [173    61   560   420];

%% Insertion Phase
% fig7 = figure(7);
% plot(f/1e9, rad2deg(angle(S21)));
% grid on

%% Rollett (k-delta) stability test
% fig111 = figure(111);
% delt = S11.*S22 - S12.*S21;
% K = (1 - abs(S11).^2 - abs(S22).^2 + abs(delt).^2)./(2*abs(S12.*S21));
% plot(f/1e9, K)
% hold on
% plot(f/1e9, abs(delt))
% grid on
% ylim([0.5 1.5])
% legend({'K', '|\Delta|'})

maxGaindB = -4.343*(Rse./Zc - G_max*Zc);
disp(['Max possible gain = ' num2str(maxGaindB) ' dB'])

%% Plot group velocity
% % vg = (diff(Betap)./diff(2*pi*f)).^(-1);
% vg = -diff(unwrap(angle(S21)))./diff(2*pi*f);
% figure(121)
% plot(f(1:end-1)/1e9, vg)

%% Plot Bsh to check change in transition frequency
% figure(122)
% fsh = 1/(2*pi*sqrt(CR*LL));
% % plot(f/1e9, real(Ysh))
% hold on
% plot(f/1e9, imag(Ysh))
% hold on
% xline([f1 f2 fsh]/1e9)
% grid on
% xlim([f1/1e9 - 1, f2/1e9 + 1])


function [idx, delta] = closestIdx(vec, val)
    [delta, idx] = min(abs(vec - val));
end