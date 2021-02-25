clear; clc; close all;
%% Setup Everything

% Add the submodules to path
addpath(genpath('OFDM-Matlab'))
addpath(genpath('WARPLab-Matlab-Wrapper'))
addpath(genpath('Power-Amplifier-Model'))

%forward distorter (amp simulation)
table = csvread('X:\dev\imc\skippy\gr-predistort\apps\barrett-fwd.csv');
board = LutModel(table(:,end), table(:,1)+1j*table(:,2), 'linear');

%linear LUT amp model
linboard = LutModel(linspace(0,0.13,101)', linspace(0,0.13,101)', 'linear');

%a "not quite right" predistorter LUT. we use a wonky one to prove the
%ILA DPD can pick up the slack.
pdtable = csvread('X:\dev\imc\skippy\gr-predistort\apps\on-14500-0.13.csv');
%pdtable = csvread('X:\dev\imc\skippy\gr-predistort\apps\bypass-14500-0.13.csv');
wonkypd = LutModel(pdtable(:,3), pdtable(:,1)+1j*pdtable(:,2), 'linear');

%the board we're actually transmitting through -- a cascaded LUT
%predistorter (with error) then through the amp in question
casc_board = wonkypd;
casc_board = casc_board.cascade(board);

Fs = 390625 / 4;
ampl = 0.128;

% two-tone test
nsamps = 150000;
f1 = 5e3;
f2 = -6.001e3;
ts = linspace(0, nsamps/Fs, nsamps)';
txd = (ampl/2)*exp(2*pi*1j*ts*f1)+(ampl/2)*exp(2*pi*1j*ts*f2);
%tx_data = Signal(txd, Fs);

%for both qpsk and 16qam
ebw = 0.55;
bw = 7500;
sps = round(Fs / bw);
nsyms = round(nsamps / sps);
filt_sym = 30;
rrc_coeffs = rcosdesign(ebw, filt_sym, sps, 'sqrt');

%qpsk test
bits = randi([0,1],nsyms,2)*2-1;
iqbits = complex(bits(:,1),bits(:,2)) * 2^-0.5;
bbbits = repmat(iqbits,1,sps)';
bbbits = bbbits(:);
bbsyms = conv(rrc_coeffs, bbbits);
bbsyms = bbsyms * ampl / max(abs(bbsyms));
%tx_data = Signal(bbsyms, Fs);

%pi/4 dqpsk test
syms = randi([1,4],round(nsyms/2),1);
phasemap = [pi/4, 3*pi/4, -pi/4, -3*pi/4];
phasechanges = phasemap(syms);
phases = cumsum(phasechanges);
iqsyms = exp(1j*phases)';
bbbits = repmat(iqsyms,1,sps)';
bbbits = bbbits(:);
bbsyms = conv(rrc_coeffs, bbbits);
bbsyms = bbsyms * ampl / max(abs(bbsyms));
%tx_data = Signal(bbsyms, Fs);

%16QAM test
%make matrix with the constellation points in it
ipoints = linspace(-1,1,4);
qpoints = 1j*linspace(-1,1,4)';
ipoints = repmat(ipoints,4,1);
qpoints = repmat(qpoints,1,4);
constellation = ipoints+qpoints;
flat_const = constellation(:);
syms = randi([1,16], nsyms, 1);
bbsyms = flat_const(syms)';
bbrep = repmat(bbsyms, 1, sps)';
filtsyms = conv(rrc_coeffs,bbrep);
filtsyms = filtsyms * ampl / max(abs(filtsyms));
tx_data = Signal(filtsyms, Fs);

casc_board.noise = 1e-17; %TODO should reference to signal power
casc_sig = Signal(casc_board.transmit(tx_data.data), Fs);

PAPR = tx_data.calculate_current_papr();
fprintf('PAPR: %.2fdB\n', PAPR);
Psat = 1.4e3; % saturation power
%these are ONLY valid for a Barrett amp operated with 0.13 peak amplitude.
%also this is approximate, don't count on it
Pavg = Psat * 10^(-PAPR/10);
fprintf('Pavg: %.0fW\n', Pavg);

%coordinates on the plot for the proposed MIC spectral mask
mask_f = [-Fs/2, -30e3, -30e3, -6e3, -6e3, 6e3, 6e3, 30e3, 30e3, Fs/2];
mask_p = [-60,   -60,   -53,   -53,   0,   0,   -53, -53,  -60,  -60]; 

% Setup DPD
dpd_params.order = 5;
dpd_params.memory_depth = 1;
dpd_params.lag_depth = 0;  % 0 is a standard MP. >0 is GMP.
dpd_params.nIterations = 5;
dpd_params.learning_rate = 0.8;
dpd_params.learning_method = 'ema'; % 'newton' or 'ema' for exponential moving average.
dpd_params.use_even = true;
dpd_params.use_conj = 0;    % Conjugate branch. Currently only set up for MP (lag = 0)
dpd_params.use_dc_term = 0; % Adds an additional term for DC
dpd = ILA_DPD(dpd_params);

%% We learn on the CASCADED board! I.e., through the imperfect LUT predistorter,
%% and then the PA itself.
dpd.perform_learning(tx_data, casc_board);
predistorted_out = dpd.predistort(tx_data.data);
w_dpd = Signal(casc_board.transmit(predistorted_out), Fs);
w_lut_dpd = Signal(casc_board.transmit(tx_data.data), Fs);
w_out_dpd = Signal(board.transmit(tx_data.data), Fs);

fprintf('Peak ILA predistorter output value: %.3f\n', max(abs(predistorted_out)));
fprintf('Peak amplifier input: %.3f\n', max(abs(w_dpd.data)));

before = w_out_dpd.measure_all_powers;
after = w_dpd.measure_all_powers;

%% Plot
p1 = w_out_dpd.plot_psd('blue');
p2 = w_dpd.plot_psd('magenta');
p3 = tx_data.plot_psd('red');
p4 = w_lut_dpd.plot_psd('green');
line(mask_f/1e3, mask_p, 'LineWidth', 2, 'Color', 'black');
legend('Without DPD', 'DPD+LUT', 'Ideal', 'With only LUT');
dpd.plot_history;


%% Plot the memoryless amplitude response of the system. This will not work
%% for systems with memory.
figure;

dcterm = 0;
if dpd_params.use_dc_term == 1
    dcterm = dpd.coeffs(end);
    dpd.coeffs = dpd.coeffs(1:end-1);
end

polyvals = flip(dpd.coeffs);
if dpd_params.use_even == 0
    polyvals = polyvals';
    polyvals(2,:) = 0;
    polyvals = polyvals(1:end-1)';
end

polyvals = vertcat(polyvals, dcterm);

x = linspace(0,ampl,301);
learned_curve = polyval(polyvals', x);
ideal_pd = casc_board.invert();
actual_curve = ideal_pd.transmit(x');
initial_pd = wonkypd.invert();
initial_curve = initial_pd.transmit(x');
plot(x, abs(learned_curve), x, abs(actual_curve), x, initial_curve);
legend('Learned', 'Actual', 'Initial', 'Location', 'southeast');