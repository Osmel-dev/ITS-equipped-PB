function [h,hNorm] = channelArrayToDevice(posAnt,bgAnt,wavelength,M,K,devPos)
    % This function computes the channel coefficients between the
    % ITS/antenna array (for the )

    %% IRS to user channel (OK)
    h = zeros(M,K);
    hNorm = zeros(M,K);
    for k = 1:K
        for m = 1:M
            % distance k^th device to m^th antenna in the array
            r = norm(devPos(:,k) - posAnt(:,m));
            
            % gain 
            theta = acos(devPos(3,k)/r);
            F = (theta>=0 & theta<=pi/2)*2*(bgAnt+1)*cos(theta)^bgAnt;
            
            % channel coefficient array -to- device
            h(m,k) = sqrt(F)*wavelength/(4*pi*r)*exp(-1i*2*pi*r/wavelength);

            hNorm(m,k) = sqrt(F)*exp(-1i*2*pi*r/wavelength);
        end
    end
end

