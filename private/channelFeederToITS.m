function T = channelFeederToITS(posAnt,N,M,bgAnt,bgFeeder,rhoITS,wavelength,da)
    % This function computes the channel coefficients between the feeder
    % and the ITS.

    delta = wavelength;
    dxy = wavelength/2;

    %% radius of the feeder
    if N == 1
        ra = 0;
    else
        ra = delta/(2*sin(pi/N));
    end  

    %% position of the feeders 
    posFeeder = zeros(3,N);
    for n = 1:N
        posFeeder(1,n) = ra*cos(pi/2+(n-1)*2*pi/N);
        posFeeder(2,n) = ra*sin(pi/2+(n-1)*2*pi/N);
        posFeeder(3,n) = -da;
    end
    
    %% active antennas -to- ITS channel coefficients
    T = zeros(M,N);
    for n = 1:N
        for m = 1:M
            % distance n^th feeder antenna to m^th passive element
            dFeedITS = norm(posFeeder(:,n) - posAnt(:,m));
    
            % ITS-to-feeder (passive antenna gain)
            theta1 = acos(da/dFeedITS);
            FAnt = (theta1>=0 & theta1<=pi/2)*2*(bgAnt+1)*cos(theta1)^bgAnt;
    
            % feeder-to-ITS (active antenna gain)
            dFeed = norm(posFeeder(:,n));
            dITS = norm(posAnt(:,m));
            cos_theta2 = (dFeed^2 + dFeedITS^2 - dITS^2)/(2*dFeed*dFeedITS);
            FFeeder = 2*(bgFeeder+1)*cos_theta2^bgFeeder;
    
            % channel coefficient feeder -to- ITS
            T(m,n) = wavelength*sqrt(rhoITS*FAnt*FFeeder)*exp(-1i*2*pi*dFeedITS/wavelength)/(4*pi*dFeedITS);
        end
    end
end