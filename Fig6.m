% This MATLAB script reproduces the results shown in Fig. 6 of 

% O. Martínez Rosabal, O. L. Alcaraz López, V. D. Pegorara Souto, 
% R. D. Souza, S. Montejo-Sánchez, R. Schober, and H. Alves, "Wireless 
% Energy Transfer Beamforming Optimization for Intelligent Transmitting 
% Surface," accepted for publication in IEEE Transactions on Wireless 
% Communications.

% For further information, visit: https://arxiv.org/pdf/2507.06805
%
% This is version 1.00 (Last edited: 2025-10-20)
%
% License: This code is licensed under the MIT license. If you in any way
% use this code for research that results in publications, please cite our
% article as described above.

clear;

%% Simulation settings
M = 10^2;                       % num. ITS elements/antennas
N = 4:10;                       % num. RF chains
fc = 5e9;                       % operating frequency [Hz]
wavelength = 3e8/fc;            % wavelength of the transmitted signal [m]
dxy = wavelength/2;             % antennas' separation [m]
da  = 4*dxy*sqrt(M/pi);         % distance feeder-to-ITS

K = 4;                          % num. IoT devices
p = 1e-3*ones(K,1);             % devices' power requirements [W]

rhoITS = .45;                   % ITS's efficiency [linear]
Pmax = 300;                     % HPA's max. output power [W]
etaMax = .25;                   % HPA's max. efficiency [linear]
l = 2;                          % l-way parameter Doherty HPA
g = 100;                        % HPA's power gain [linear]

bgFeeder = 10;                  % feeder's antennas boresight gain [linear]
bgAnt = 2;                      % ITS elements boresight gain [linear]

Pctrl = 1;                      % control power of the ITS [W]
Pcell = 1e-3;                   % control power of the ITS elements [W]

tol = 1;                        % tolarance of the SCA algorithm [W]
MCIter = 100;                   % number of Monte Carlo runs

RFNetLoss.Ld = 10^(0.5/10);     % insertion losses power splitters [linear]
RFNetLoss.Lc = 10^(0.5/10);     % insertion losses power combiners [linear]
RFNetLoss.Lp = 10^(3.5/10);     % insertion losses phase shifteres [linear]
   
% Positions of the antennas in the array/ITS
posAnt = positionArray(M,dxy);

%% Monte Carlo loop

% memory pre-allocation
HPAPowITS = zeros(numel(N),1);
HPAPowHBFC = zeros(numel(N),1);
HPAPowHBPC = zeros(numel(N),1);
HPAPowFD = 0;

for seed = 1:MCIter
    % devices positions
    rng(seed)
    devPos = rand(3,K);
    
    devPos(1,:) = 3*devPos(1,:) - 1.5;
    devPos(2,:) = 3*devPos(2,:) - 1.5;
    devPos(3,:) = 5;

    % channel coefficients ITS/array -to- device
    h = channelArrayToDevice(posAnt,bgAnt,wavelength,M,K,devPos);

    HPAPowITS_ = zeros(numel(N),1);
    HPAPowHBFC_ = zeros(numel(N),1);
    HPAPowHBPC_ = zeros(numel(N),1);
    for n = 1:numel(N)
        disp(n)
        % channel coefficients feeder -to- ITS
        T = channelFeederToITS(posAnt,N(n),M,bgAnt,bgFeeder,rhoITS,wavelength,da);
    
        % ITS-assisted PB
        [HPAPowITS_(n),~,~] = ITSAssistedPB(h,T,N(n),p,Pmax,l,etaMax,g,K,M,tol);
    
        % PB equipped with a FC hybrid analog-digital beamforming architecture
        [HPAPowHBFC_(n),~,~] = HBFC(h,N(n),p,Pmax,l,etaMax,g,K,M,tol,RFNetLoss);
    
        % PB equipped with a PC hybrid analog-digital beamforming architecture
        [HPAPowHBPC_(n),~,~] = HBPC(h,N(n),p,Pmax,l,etaMax,g,K,M,tol,RFNetLoss);
    end

    % PB equipped with a FD beamforming architecture
    [HPAPowFD_,~] = fullyDigital(h,M,p,Pmax,l,etaMax,g,K,tol);

    % average results
    HPAPowITS = HPAPowITS + 1/MCIter*HPAPowITS_;
    HPAPowHBFC = HPAPowHBFC + 1/MCIter*HPAPowHBFC_;
    HPAPowHBPC = HPAPowHBPC + 1/MCIter*HPAPowHBPC_;
    HPAPowFD = HPAPowFD + 1/MCIter*HPAPowFD_;
end

%% Plot results
fig1 = figure(1);
plotSettings(fig1);

plot(N,10*log10(HPAPowITS),'-o','Color','#0072BD','LineWidth',2); hold on
plot(N,10*log10(HPAPowFD)*ones(numel(N),1),'-.s','Color','#7E2F8E','LineWidth',2);
plot(N,10*log10(HPAPowHBFC),'--^','Color','#A2142F','LineWidth',2);
plot(N,10*log10(HPAPowHBPC),':p','Color','#77AC30','LineWidth',2); hold off

grid on
box on
ax = gca;
ax.FontSize = tckFontSize; 
ax.TickLabelInterpreter = 'latex';
xlabel('$N$','FontSize',15,'Interpreter','latex')
ylabel('$P_T$ (dBW)','FontSize',15,'Interpreter','latex')
legend('ITS-equipped PB','FD-equipped PB','HBFC-equipped PB','HBPC-equipped PB',...
    'Position',[0.135063539696936 0.23064834482447 0.301403559457635 0.209933769781858],...
    'FontSize',14,'Interpreter','latex')