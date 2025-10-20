function [HPAPow,phiOpt,bOpt] = ITSAssistedPB(h,T,N,p,Pmax,l,etaMax,g,K,M,tol)
    % This function minimizes the power consumption of a PB equipped with
    % an ITS.

    % INPUT:
    % h           => channel coeffs. ITS-to-devs. [MxK].
    % T           => channel coeffs. feeder-to-ITS [MxN].
    % N           => num. RF chains [scalar].
    % p           => received pow. requirement [Kx1], [W].
    % Pmax        => HPA max. output power, scalar, [W].
    % l           => l-way Doherty HPA parameter, [scalar].
    % etaMax      => HPA max. efficiency, [scalar], [linear].
    % g           => HPA's power gain, [scalar], [linear].
    % K           => number of deployed IoT devices, [scalar].
    % M           => number of transmit antennas, [scalar].
    % tol           => tolerance SCA [W].

    % OUTPUT:
    % HPAPow      => HPA's minimum output power, [scalar], [W].
    % phiOpt      => phase shifter configuration, [Mx1].
    % bOpt        => digital beamforming, [NxQ].

    % For further information, visit: https://arxiv.org/pdf/2507.06805
    %
    % This is version 1.00 (Last edited: 2025-10-20)
    %
    % License: This code is licensed under the MIT license. If you in any way
    % use this code for research that results in publications, please cite our
    % article as described above.
    

    %% Initialization procedure 
    
    % Assuming that N >= K, divide the ITS elements into N
    % non-overlapping single-RF chain PBs.
    idxITSElementsPerAntenna = zeros(M,1);
    for m = 1:M
        [~,idx] = max(abs(T(m,:)).^2);
        idxITSElementsPerAntenna(m) = idx;
    end

    % Compute the number of RF chains assigned to each device.   
    pathLoss = vecnorm(h)'.^2;
    numRFChainsPerDev = ones(K,1);
    for n = 1:N-K
        metric = zeros(K,1);
        for k = 1:K
            metric(k) = numRFChainsPerDev(k)^2 + sum(numRFChainsPerDev) - numRFChainsPerDev(k);
        end

        [~,idx] = min(pathLoss.*metric);
        numRFChainsPerDev(idx) = numRFChainsPerDev(idx) + 1;
    end

    devID = 1:K;
    idxRFChainsPerDeviceVec = [];
    for k = 1:K
        idxRFChainsPerDeviceVec = [idxRFChainsPerDeviceVec repmat(devID(k), 1, numRFChainsPerDev(k))];
    end

    % Compute the phase shifters configuration
    idxRFChainsPerDeviceMx = unique(perms(idxRFChainsPerDeviceVec),"rows");
    totalPermutations = size(idxRFChainsPerDeviceMx,1);

    % memory preallocation
    HPAPowVec = zeros(totalPermutations,1);
    bCell = cell(totalPermutations,1);
    phiCell = cell(totalPermutations,1);

    for ii = 1:totalPermutations
        % select the i-th RF chains to devices assignment
        idxDevice = idxRFChainsPerDeviceMx(ii,:);

        % feeder-to-ITS channel slice
        ITSElementsMaskFeeder = (idxITSElementsPerAntenna == (1:N));
        TSliced = T.*ITSElementsMaskFeeder;
        TSliced = sum(TSliced,2);

        % ITS-to-devs channel slice
        ITSElementsMaskDevices = arrayfun(@(x) sum(ITSElementsMaskFeeder(:,idxDevice == x),2),...
            unique(idxDevice), 'UniformOutput', false);
        ITSElementsMaskDevices = horzcat(ITSElementsMaskDevices{:});
        hSliced = h.*ITSElementsMaskDevices; 
        
        % ITS phase shift MRT
        phiTemp = exp(-1i*sum(-angle(hSliced) + angle(TSliced),2));

        % Compute the digital beamforming SDP formulation    
        [bTemp, HPAPow] = ITSAssistedPBInit(K,g,h,T,p,N,Pmax,l,etaMax,phiTemp);

        HPAPowVec(ii) = HPAPow;
        bCell{ii} = bTemp;
        phiCell{ii} = phiTemp;
    end

    % Select the init solution
    [~,idx] = min(HPAPowVec);
    bInit = bCell{idx};
    phiInit = phiCell{idx};   

    %% main optimization loop (SCA)

    % init. optimization variables
    Q = N; 

    z = zeros(K,Q);
    nu = cell(K,Q);
    for k = 1:K
        hTilde = (h(:,k)'*diag(phiInit)*T)';
        for q = 1:Q        
            z(k,q) = hTilde'*bInit(:,q);
            nu{k,q} = z(k,q)*hTilde + bInit(:,q);
        end
    end

    % SCA loop
    HPAPowPrev = 0;
    HPAPowNext = Inf;
    iter = 0;
    while (abs(HPAPowNext - HPAPowPrev) > tol)
        iter = iter + 1;
        disp(iter)

        cvx_begin quiet
            variable t(N) nonnegative
            variable phi(M) complex
            variable b(N,Q) complex
            minimize ( sum(t) )
            subject to 

                for n = 1:N                        
                    % epigraph objective
                    if g*norm(bInit(n,:))^2 <= Pmax/l^2
                        sqrt(g*Pmax)/(l*etaMax)*norm(b(n,:)) <= t(n)

                        g*square_pos(norm(b(n,:))) <= Pmax/l^2
                    else
                        (l+1)*sqrt(g*Pmax)/(l*etaMax)*norm(b(n,:)) - Pmax/(l*etaMax) <= t(n)

                        g*square_pos(norm(b(n,:))) <= Pmax
                        2*g*real(conj(bInit(n,:))*b(n,:).') + g*norm(bInit(n,:))^2 >= Pmax/l^2
                    end

                    for k = 1:K
                        % effective channel
                        hTilde = (h(:,k)'*diag(phi)*T)';

                        temp2 = 0;
                        for q = 1:Q
                            temp2 = temp2 + real(nu{k,q}'*(z(k,q)*hTilde + b(:,q))) - ...
                                1/2*norm(nu{k,q})^2 - 1/2*square_pos(norm(z(k,q)*hTilde - b(:,q))) - abs(z(k,q))^2;
                        end

                        g*temp2/p(k) >= 1
                    end
                end

                for m = 1:M
                    abs(phi(m)) <= 1
                end
        cvx_end

        % update optimization variables

        % return the previous value before the optimization returned
        % NaN/Inf
        if isnan(cvx_optval) || isinf(cvx_optval)
            break;
        end

        bInit = b;
        phiInit = phi;

        for k = 1:K
            hTilde = (h(:,k)'*diag(phi)*T)';
            for q = 1:Q        
                z(k,q) = hTilde'*b(:,q);
                nu{k,q} = z(k,q)*hTilde + b(:,q);
            end
        end

        HPAPowPrev = HPAPowNext;

        HPAPowNext = 0;
        for n = 1:N
            %
            if g*norm(b(n,:))^2 <= Pmax/l^2
                HPAPowNext = HPAPowNext + sqrt(g*norm(b(n,:))^2*Pmax)/(l*etaMax);
            else
                HPAPowNext = HPAPowNext + ((l+1)*sqrt(g*norm(b(n,:))^2*Pmax) - Pmax)/(l*etaMax);
            end
        end
    end

    %% results
    fprintf('ITSAssistedPB algorithm successfully converged to optimal value %.4e after %d iterations. Number of antennas and ITS elements is %d and %d, respectively.\n', ...
            HPAPowNext, iter, N, M);

    HPAPow = HPAPowNext;
    phiOpt = phiInit;
    bOpt = bInit;
end

function [bTemp, HPAPow] = ITSAssistedPBInit(K,g,h,T,p,N,Pmax,l,etaMax,phiTemp)
    
    % optimization problem: minimum the maximum power consumption  
    % channel matrix (SDP formulation) 
    HTilde = zeros(N,N,K);
    for k = 1:K
        hTilde = (h(:,k)'*diag(phiTemp)*T)';
        HTilde(:,:,k) = hTilde*hTilde';
    end

    % basis vectors in rank-1 matrix form (SDP formulation)
    E = zeros(N,N,K);
    e_ = eye(N);
    for n = 1:N
        E(:,:,n) = e_(:,n)*e_(:,n)';
    end

    cvx_begin sdp quiet
    variable B(N,N) hermitian semidefinite
    variable t nonnegative
    minimize ( t )
    subject to 
    for k = 1:K
        g*real(trace(HTilde(:,:,k)*B))/p(k) >= 1
    end
    
    for n = 1:N
        g*trace(E(:,:,n)*B) <= t
    end
    cvx_end
    
    % solution recovery via eigenvalue decomposition
    r = rank(B);
    [U,S] = eig(B);
    bTemp = zeros(N,r);
    for rr = 1:r
        bTemp(:,rr) = sqrt(S(rr,rr))*U(:,rr);
    end
    
    HPAPow = 0;
    for n = 1:N
        %
        if g*norm(bTemp(n,:))^2 <= Pmax/l^2
            HPAPow = HPAPow + sqrt(g*norm(bTemp(n,:))^2*Pmax)/(l*etaMax);
        else
            HPAPow = HPAPow + ((l+1)*sqrt(g*norm(bTemp(n,:))^2*Pmax) - Pmax)/(l*etaMax);
        end
    end
end