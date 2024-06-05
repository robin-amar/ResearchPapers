% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Copyright (c) 2024, <Robin Amar>, <Mohammad Alaee>, website: https://radarmimo.com/
% All rights reserved.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright
%    notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
% 3. All advertising materials mentioning features or use of this software
%    must display the following acknowledgement:
%    This product includes software developed by the <organization>.
% 4. Neither the name of the UniLU nor the
%    names of authors may be used to endorse or promote products
%    derived from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY <COPYRIGHT HOLDER> ''AS IS'' AND ANY
% EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
% DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
% ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Email: robin.amar@uni.lu, amar.robin.2020@gmail.com
% Code written by : Robin Amar
% Last update : 2024
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% Some typical parameters based on TI-AWR mmWave Sensors
% – TX power: 12 dBm
% – RX noise figure:
% – 14 dB (76 to 77 GHz)
% – 15 dB (77 to 81 GHz)
% -------------------------------------------------------------------------
%                               LRR             MRR             SRR
% -------------------------------------------------------------------------
% Max unambiguous range         225 m           125 m           45 m
% Bandwidth                     300 MHz         540 MHz         750 MHz
% Ramp slope                    10 MHz/us       12 MHz/us       15 MHz/us
% Inter-chirp/Pulse duration    8 us            10 us           12 us
% Number of chirps              256             128             128
% Range resolution              0.5 m           0.28 m          0.2 m
% Chirp/Pulse duration          30 us           45 us           50 us
% Max umambiguous velocity      92.28 km/h      63.75 km/h      56.56 km/h
% Frame time (total)            9.728 ms        7.04 ms         7.94 ms
% Frame time (active)           7.68 ms         5.76 ms         6.4 ms
% Radar data memory required    2048 KB         1024            512
% -------------------------------------------------------------------------
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

clc; clear; close all

rng(2020);
MarkerSize  = 10;
LineWidth   = 2;
FontSize    = 14;

RNDM_SEQUENCE = 2; RANDOM_SLOPE_FMCW = 3; CONSTANT_SLOPE_FMCW = 4; OPTIMAL_SLOPE_MPIR = 5; RANDOM_QUAD_COEFF = 6;
AWGN_TYPE = 1; WGN_TYPE = 2;

DEBUG           = false;
pltSpecRange    = false;
CODE_TYPE       = OPTIMAL_SLOPE_MPIR;
% CODE_TYPE       = RANDOM_QUAD_COEFF;


%% Target Model
target_dist     = 30;
target_speedkmh = 10;
target_rcs      = 10;

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%% Sensor parameters
% Radar System Setup
ant_gain    = 10;                                     % in dB
tx_dB       = 12;                                     % in watts
tx_power    = 10^(tx_dB/10)*1e-3;                     % Tx Power: 12dBm in watts
rx_nf_dB    = 15;                                     % in dB
rx_gain_dB  = 15;

fc          = 79e9;
c           = 3e8;

Nsweep      = 128;

tx_duration = 20e-6;
idle_time   = 0;

pfa         = 1e-2;                                 % Probability of False Alram

k           = 1.38e-23;
T           = 300;
loss        = 0;                                    % dB

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
range_res       = 1;
lambda          = c/fc;
bw              = c / (2 * range_res);
tau_chip        = 1 / bw;                          % Chip Duration. Duration for which a code (say +1 or -1) is applied!
fs              = bw;                                    % Because of IQ Sampling

fc_dummy        = bw;
Ts              = 1/(fc_dummy);

f_adc_FMCW      = 10e6;
Slope_FMCW      = bw/tx_duration;
Rmax_FMCW       = f_adc_FMCW * c /(2*Slope_FMCW);
Rmax_PMCW       = Rmax_FMCW;
% ----------------------------------------------------
% Note: Code Length is deliberately kept lower s.t. it fits the requirement
% of txsig length (i.e. numTxsamples+numIdlesamples)
codeLength      = round(tx_duration*bw);                  % Phase Code
% tx_time         = tx_duration + idle_time;          % Overall Transmit duration

Tp = 60e-6;
tx_time         = Tp;
fp              = 1 / tx_time;                                           % PRF
target_doppler  = 2 * target_speedkmh / 3.6 / lambda;        % Doppler Freq = 2.Target_Velocity/Lambda

speed_maxKmh    = lambda/(4 * tx_time) *3.6;                   % km/s
range_max       = codeLength*range_res;
maxRngSmples    = round(Rmax_PMCW/range_max * codeLength);

numTxsamples = round(tx_duration/tau_chip);     % Phase Code is applied for the whole Tx duration
numIdlesamples = round(idle_time/tau_chip);
totalsamples = numTxsamples+numIdlesamples;
N               = totalsamples;
% ----------------------------------------------------
switch CODE_TYPE
    case RNDM_SEQUENCE
        N           = codeLength;
        code_ref    = mcode(8, N);     
        code        = code_ref;
        for m = 1:Nsweep
            sig(:,m)        = code(:);
        end

%         totalsamples = numTxsamples+numIdlesamples;
    case RANDOM_SLOPE_FMCW
        Tc_min          = 10e-6;
        Tc_max          = 30e-6;
        Tc              = tx_duration;
        Tc_var          = Tc_min + (Tc_max-Tc_min).*rand(Nsweep,1);
        t               = 0:Ts:Tc-Ts;
        sweep_slope     = bw./Tc_var;
        zeta            = randi([-1, 1],Nsweep,1);
        zeta(~zeta)     = 1;
        for m = 1:Nsweep
            theta       = 2*pi*(fc_dummy*t + 0.5*zeta(m)*sweep_slope(m)*t.^2);
            x           = exp(1i*theta); sig(:,m) = x(:);
        end
%         totalsamples    = length(t);
        N               = totalsamples;
    case CONSTANT_SLOPE_FMCW
        Tc              = 20e-6;
        t               = 0:Ts:Tc-Ts;
        sweep_slope     = bw./Tc;
        theta           = 2*pi*(fc_dummy*t + 0.5*sweep_slope*t.^2);
        x               = exp(1i*theta); sig_cnst = x(:);
        totalsamples    = length(t);
%         N               = totalsamples;

        for m = 1:Nsweep
            sig(:,m)        = x(:);
        end

    case RANDOM_QUAD_COEFF
        aQL = 1*randn(3, Nsweep);
        a2_rnd_coeff = 0.001*randn( Nsweep,1);
        a1_rnd_coeff = 0.000001*randn( Nsweep,1);
        a0_rnd_coeff = 1*randn( Nsweep,1);
        N               = totalsamples;
        n               = 1:N;
        for m = 1:Nsweep
            theta_r         = a0_rnd_coeff(m) + a1_rnd_coeff(m)*n + a2_rnd_coeff(m)*n.^2;
            x               = exp(1i*theta_r);
            sig(:,m)        = x(:);
        end

    case OPTIMAL_SLOPE_MPIR
        Tc              = 20e-6;
        t               = 0:Ts:Tc-Ts;
        sweep_slope     = bw./Tc;
        % ------------
        MM = 1; CD = 2; 
        MTHD = MM;
        switch MTHD
            case MM
%                 objct           = load('MPSL_N640k_M2500L256Q2_lp2.mat'); 
                objct           = load('MPSL_M100L128Q2_lp2.mat');
                aQL             = objct.pStruct.aQL;
                a2_quad_coeff   = aQL(3,:);
                a1_quad_coeff   = aQL(2,:);
                a0_quad_coeff   = aQL(1,:);
            case CD
%                 objct           = load('sig_MhmdN2500M100L25Q3.mat');
                objct           = load('sig_MhmdN12800M100L128Q3.mat');
                
                for m = 1:Nsweep
                    a2_quad_coeff(m)   = objct.aseq(m).a(3);
                    a1_quad_coeff(m)   = objct.aseq(m).a(2);
                    a0_quad_coeff(m)   = objct.aseq(m).a(1);
                end
                
        end
        % ------------
        N           = codeLength;
        n               = 1:N;
        for m = 1:Nsweep
            snippet = 2;
            
            switch snippet
                case 1
                    a0              = 0;
                    a1              = 2*pi*bw*Ts;
                    a2              = pi*(bw/Tc)*Ts^2;
                    theta           = a0 + a1*n + a2*n.^2;
                    x               = exp(1i*theta);
                case 2            
                    theta_1          = a0_quad_coeff(m) + a1_quad_coeff(m)*n + a2_quad_coeff(m)*n.^2;
                    x               = exp(1i*theta_1);
                case 3
                    x               = objct.x_PSL;
                case 4
                    bwE = a1_quad_coeff(m)/(2*pi*Ts);
                    TcE = 0.5*(bwE/a2_quad_coeff(m))*Ts^2;
            end
            sig(:,m)        = x(:);
%             % =========== Spectrogram Test ============
%             figure(100);
%             subplot(211); plot(real(x(:)));
%             xlabel('Time (us)'); ylabel('Amplitude (v)');
%             title('Chirp signal'); axis tight;
%             subplot(212); spectrogram(x(:), round(0.2*N),round(0.1*N),'yaxis');
%             % ===========
        end       
end
figure(100);spectrogram(sig(:), round(0.01*N),0,[],100*fs,'yaxis');
txsig = zeros(totalsamples,1);

% ----------------------------------------------------
% min_snr = albersheim(0.9,pfa);                              % pfa -> Probability of False Alarm, 0.9 -> Probability of Detection
min_snr = 11;

%% Maximum detectable range

range = 1 : Rmax_PMCW;
pn = 10*log10(k * T * bw * 10^(rx_nf_dB/10)) + 30; %dBm                % 30 was added to convert the Noise power in dBm
if pltSpecRange
    figure
    pr = 10*log10(tx_power * 10^(ant_gain/10) * 10^(ant_gain/10) ...
        * lambda^2 * target_rcs ./ ( (4*pi)^3 * range.^4 * 10^(loss/10))) + ...
        30 ;
    plot(range, pr, range, pr+ 10*log10(codeLength) , ...               % pr+ 10*log10(codeLength) -> Coherent Integration Gain per chirp
        range, pr + 10*log10(codeLength) + 10 * log10(Nsweep));         % 10 * log10(Nsweep) -> Coherent Integration gain across chirps
    hold on
    plot(range,pn*ones(1,length(range)),range, (pn + min_snr)*ones(1,length(range)));
    xlabel('range (m)')
    ylabel('P_r (dBm)')
    grid on
    title(['RCS = ', num2str(target_rcs)])
    legend('P_r', 'P_r+MF','P_r+MF+FFT','noise','threshold')
end

%% Simulation loop
% Next, run the simulation loop.
nfft_d = 2^(nextpow2(Nsweep));
indexSpeed =  round(speed_maxKmh / (fp*lambda/4*3.6) * nfft_d);
if indexSpeed < 64
    indexSpeed = min(64, nfft_d);
end

indexRange = round(Rmax_PMCW/range_res);
xr = complex(zeros(indexRange,Nsweep));
DSTV = zeros(Nsweep,1);                       % DSTV - Doppler SteeringVector
DSTV(:) = dopsteeringvec(-target_doppler,Nsweep,fp);
vr = sqrt((tx_power * 10^(ant_gain/10) * 10^(ant_gain/10) ...
    * lambda^2 * target_rcs ./ ( (4*pi)^3 * target_dist^4 * 10^(loss/10))) );

for m = 1:Nsweep
    switch CODE_TYPE
        case RNDM_SEQUENCE
            txsig(1:codeLength) = code(1:codeLength);
        case RANDOM_SLOPE_FMCW
            txsig = sig(:,m);
        case CONSTANT_SLOPE_FMCW
            txsig = sig_cnst;
        case OPTIMAL_SLOPE_MPIR
            txsig = sig(:,m);
        case RANDOM_QUAD_COEFF
            txsig = sig(:,m);
    end
    

    % Receive PMCW waveform
    rxsig = zeros(size(txsig));
    rxsig = rxsig + 2 * vr * circshift(txsig,round(target_dist/range_res)) * DSTV(m); % Circshift provides time delay
    
    NOISE_METHOD = WGN_TYPE;
    switch NOISE_METHOD
        case AWGN_TYPE
            y = awgn(rxsig,min_snr,'measured');
            rxsig = y;
        case WGN_TYPE
            % Independent Noise
            noiseSig = wgn(totalsamples,1,pn,'complex','dBm');
            rxsig = rxsig + noiseSig;
    end
    % matched filter
%     xMF = xcorr(txsig,rxsig)/(2*codeLength);
    xMF = xcorr(txsig,rxsig);

    if DEBUG
        rngGrid = linspace(0,Rmax_PMCW,indexRange);
        win     = dsp.Window('WindowFunction', 'Hamming');
        figure;
        plot(10*log10(abs(xMF)));
        figure;
        plot(abs(xMF));
        return
    end
    xr(:,m) = xMF(totalsamples:-1:totalsamples-indexRange+1) ;
end

%% Range and Doppler Estimation
xr = 10^(rx_gain_dB/10) * xr;

rngDop = fft(xr.', nfft_d);
rngDpow = abs(fftshift(rngDop,1));


speedGrid = linspace(- speed_maxKmh,speed_maxKmh, indexSpeed);
rngGrid = linspace(0, Rmax_PMCW, indexRange);


f1 = figure;
hold on
imagesc(rngGrid,speedGrid,10*log10(rngDpow)+30);
set(gca,'YDir','normal')
ylabel('Speed (km/h)');
xlabel('Range (m)');
box on
axis tight
cb = colorbar;
cb.Label.String = 'Power (dBm)';
% caxis([-20 30])
caxis([0 40])
% title('Range- Doppler Map')
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%% SINR Calculation
% --- Track Target Location start ----
var_rng = 1; var_dop = 1;
CUT_speedkmph = target_speedkmh;
dopCUT_idx = round(CUT_speedkmph/speed_maxKmh/2 * indexSpeed +  indexSpeed/2);
CUT_rng = target_dist;
[~, rngCUT_idx] = min(abs(rngGrid - CUT_rng));
% --- Track Target Location end ----
% ------ Target Power Level ---------------
aa = dopCUT_idx - var_dop: dopCUT_idx + var_dop;        % <- Target location Checked
bb = rngCUT_idx - var_rng: rngCUT_idx + var_rng;
[aG2, bG2] = meshgrid(aa,bb);
TgtPwr = 10*log10(max(max(rngDpow(aG2(:), bG2(:)))))+30;
rngDpowTmp = rngDpow;
% rngDpow(aG2(:), bG2(:)) = 1e-10;
rngDpowTmp(aG2(:), bG2(:)) = min(min(rngDpow));
OvrlNoiseThrs = meanabs(rngDpowTmp);
OvrlNoiseThrsdB = 10*log10(OvrlNoiseThrs)+30;
TgtSINR  = TgtPwr - OvrlNoiseThrsdB;

fprintf('Noise Threshold: %f \n', OvrlNoiseThrsdB);
fprintf('Target Power: %f \n', TgtPwr);
fprintf('Target SINR: %f \n', TgtSINR);


% --- Track Noise Location end ----

figure;
imagesc(rngGrid,speedGrid, 10*log10(rngDpowTmp)+30);
set(gca,'YDir','normal')
ylabel('Speed (km/h)');
xlabel('Range (m)');
box on
axis tight
colorbar
% caxis([-20 30])
caxis([0 40])
%% Range and Doppler CUT  %%
% RD_DplrCUT = 10*log10(rngDpow(dopCUT_idx,:));
% RD_rngCUT  = 10*log10(rngDpow(:,rngCUT_idx));
RD_DplrMTX = 10*log10(rngDpow(dopCUT_idx-2:dopCUT_idx+2,:));
RD_rngMTX  = 10*log10(rngDpow(:,rngCUT_idx-2:rngCUT_idx+2));
[MR,IR] = max(RD_DplrMTX,[],2);
[MD,ID] = max(RD_rngMTX,[],1);

[valD,posD]=intersect(MD,MR);
[valR,posR]=intersect(MR, MD);

RD_DplrCUT = RD_DplrMTX(posR,:);
RD_rngCUT  = RD_rngMTX(:,posD);

indxR = find(RD_rngCUT == valD);
indxD = find(RD_DplrCUT == valR);

figure; plot(rngGrid, RD_DplrCUT);
ylabel('Power (dBm)');
xlabel('Range (m)');
box on
axis tight
ylim([-20,15])
% figure; plot(RD_DplrCUT);
figure; plot(speedGrid, RD_rngCUT);
ylabel('Power (dBm)');
xlabel('Doppler (km/h)');
box on
axis tight
ylim([-20,15])
% figure; plot(RD_rngCUT);

%% ISL/ PSL - range/Doppler CUTs
[ dpsl,disl,dPSLRdB,dISLRdB ] = pslisl_RD_CUT( RD_DplrCUT, indxD);
fprintf('----- Range profile --------\n')
fprintf('PSL/PSLR = %f,%f\n', dpsl, dPSLRdB)
fprintf('ISL/ISLR = %f,%f\n', disl, dISLRdB)

[ rpsl,risl,rPSLRdB,rISLRdB ] = pslisl_RD_CUT( RD_rngCUT, indxR);
fprintf('----- Doppler profile --------\n')
fprintf('PSL/PSLR = %f,%f\n', rpsl, rPSLRdB)
fprintf('ISL/ISLR = %f,%f\n', risl, rISLRdB)

function [ psl,isl,PSLRdB,ISLRdB ] = pslisl_RD_CUT(vec,i)
tmp = vec;
vec(i) = min(vec);
r = vec;
N           = length(r); 
psl         = max(r);
isl         = sum(r.^2);

PSLRdB      = 10*log10(psl^2/N^2);
ISLRdB      = 10*log10(isl/N^2);
end