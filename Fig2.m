% This MATLAB script reproduces the results shown in Fig. 2 of 

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

clear; clc;

% Settings
etaMax = pi/4;                % Doherty maximum efficiency
x = linspace(0,1,100);        % normalized input power
l = 1:4;                      % number of HPAs in the Doherty architecture

% memory pre-allocation & drain efficiency computation
eta = zeros(numel(3),numel(x));
for jj = 1:numel(l)
    for ii = 1:numel(x)
        if l(jj) == 1
            % Class B amplifier
            eta(jj,ii) = l(jj)*etaMax*sqrt(x(ii));
        else
            % l-way Doherty
            if x(ii) <= 1/l(jj)^2
                eta(jj,ii) = l(jj)*etaMax*sqrt(x(ii));
            else
                eta(jj,ii) = l(jj)*etaMax*x(ii)/((l(jj)+1)*sqrt(x(ii))-1);
            end
        end
    end
end

fig1 = figure(1);
plotSettings(fig1)

step = 5;
plot(10*log10(x),eta(1,:),'-o','LineWidth',1.5,'Color','#0072BD',...
    'MarkerIndices',1:step:length(x)); hold on
plot(10*log10(x),eta(2,:),'-^','LineWidth',1.5,'Color','#7E2F8E',...
    'MarkerIndices',1:step:length(x))
plot(10*log10(x),eta(3,:),'-p','LineWidth',1.5,'Color','#A2142F', ...
    'MarkerIndices',1:step:length(x)); 
plot(10*log10(x),eta(4,:),'-s','LineWidth',1.5,'Color','#77AC30', ...
    'MarkerIndices',1:step:length(x)); hold off

ax = gca;
ax.FontSize = 12; 

grid on
xlabel('normalized input power (dB)','FontSize',15,'Interpreter','latex')
ylabel('drain efficiency','FontSize',15,'Interpreter','latex')
legend('class B','$2-$way','$3-$way','$4-$way','FontSize',14,...
    'location','southeast','Interpreter','latex') 