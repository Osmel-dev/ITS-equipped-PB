% This MATLAB script reproduces the results shown in Fig. 9 of 

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
M = 20^2;                       % num. ITS elements/antennas
N = 1;                          % num. RF chains
fc = 5e9;                       % operating frequency [Hz]
wavelength = 3e8/fc;            % wavelength of the transmitted signal [m]
dxy = wavelength/2;             % antennas' separation [m]
da  = [1.3541 0.2];               % distance feeder-to-ITS [m]

K = 1;                          % num. IoT devices
p = 1e-3*ones(K,1);             % devices' power requirements [W]
devPos = [0 0 1.5]';             % device position [m]

rhoITS = .45;                   % ITS's efficiency [linear]
Pmax = 300;                     % HPA's max. output power [W]
etaMax = .25;                   % HPA's max. efficiency [linear]
l = 2;                        % l-way parameter Doherty HPA
g = 100;                        % HPA's power gain [linear]

bgFeeder = 10;                  % feeder's antennas boresight gain [linear]
bgAnt = 2;                      % ITS elements boresight gain [linear]

Pctrl = 1;                      % control power of the ITS [W]
Pcell = 1e-3;                   % control power of the ITS elements [W]

tol = 100;                        % tolarance of the SCA algorithm [W]
MCIter = 1;                   % number of Monte Carlo runs

RFNetLoss.Ld = 10^(0.5/10);     % insertion losses power splitters [linear]
RFNetLoss.Lc = 10^(0.5/10);     % insertion losses power combiners [linear]
RFNetLoss.Lp = 10^(3.5/10);     % insertion losses phase shifteres [linear]

%% Main simulation

% Positions of the antennas in the array/ITS
posAnt = positionArray(M,dxy);
    
% channel coefficients ITS/array -to- device
h = channelArrayToDevice(posAnt,bgAnt,wavelength,M,K,devPos);

% probes grid (device side, 200 per dimension)
probesDev = 200;
[probesPosYDev,probesPosZDev] = meshgrid(linspace(-5,5,probesDev),linspace(.5,12,probesDev));

% probes grid (ITS side, 20 per dimension)
probesITS = 20;
[probesPosXITS,probesPosYITS] = meshgrid(linspace(-0.285,0.285,probesITS),linspace(-0.285,0.285,probesITS));

% memory pre-allocation
rPowDev = zeros(probesDev, probesDev, numel(da));
rPowITS = zeros(probesITS, probesITS, probesITS, numel(da));

for ii = 1:numel(da)
    % channel coefficients feeder -to- ITS
    T = channelFeederToITS(posAnt,N,M,bgAnt,bgFeeder,rhoITS,wavelength,da(ii));
    
    phiOpt = exp(-1i*sum(-angle(h) + angle(T),2));
    bOpt = sqrt(p/(g*rhoITS*abs(h'*diag(phiOpt)*T)^2));

    % ==================================================================
    % normalized received power (device side)        
    rPowDevSlice = zeros(probesDev,probesDev);
    for index = 1:numel(rPowDevSlice)
        probePos = [0 probesPosYDev(index) probesPosZDev(index)]';
    
        % Positions of the antennas in the array/ITS
        posAnt = positionArray(M,dxy);
    
        % channel coefficients ITS/array -to- device
        [~, hNorm] = channelArrayToDevice(posAnt,bgAnt,wavelength,M,K,probePos);
    
        % normalized received power
        rPowDevSlice(index) = norm(hNorm'*diag(phiOpt)*T*bOpt)^2;
    end

    rPowDev(:,:,ii) = rPowDevSlice;
      
    % convert the positions of the probes into a 3 x num. probes.
    probesPos = [probesPosXITS(:)'; probesPosYITS(:)'; zeros(1,numel(probesPosXITS))];
    
    % channel coefficients feeder -to- ITS
    T = channelFeederToITS(probesPos,N,M,bgAnt,bgFeeder,1,wavelength,da(ii));
    
    % normalized received power at the ITS
    rPowITS(:,:,ii) = reshape(abs(T).^2,[probesITS,probesITS]);
end

% global normalization
gmax = max( [max(rPowDev(:)), max(rPowITS(:))] ); 

rPowDev = rPowDev/gmax;
rPowITS = rPowITS/gmax;

%% Plot results
fig1 = figure(1);
plotSettings(fig1)

t = tiledlayout(2,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

for ii = 1:numel(da)
    % normalized received power (ITS side)
    nexttile(2*ii-1)
    pcolor(probesPosXITS,probesPosYITS,rPowITS(:,:,ii));
    box on
    shading interp;
    colormap jet;
    colorbar
    
    xlabel('y (m)')
    ylabel('x (m)')
    axis equal
    axis tight
    
    % normalized received power (device side)
    nexttile(2*ii)
    pcolor(probesPosZDev,probesPosYDev,rPowDev(:,:,ii)); hold on
    shading interp;
    colormap jet;
    colorbar
    xlabel('z (m)')
    ylabel('y (m)')
    
    % device position
    devPos = [0 0 1.5]';
    plot(devPos(3,:),devPos(2,:),'sw','MarkerSize',10,'LineWidth',1.5);
    pbaspect([1 1 1]);
end