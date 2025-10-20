function [HPAPow,bOpt,aOpt] = HBPC(h,N,p,Pmax,l,etaMax,g,K,M,tol,RFNetLoss)
    % This function minimizes the power consumption of a PB equipped with a
    % partially-connected hybrid analog-digital beamforming architecture

      % INPUT
    % h             => channel coefficients, [MxK].
    % N             => number of RF chains, [scalar].
    % Q             => number of streams, [scalar].
    % p             => received power requirements, [Kx1], [W].
    % Pmax          => HPA max. output power, [scalar], [W].
    % l             => l-way Doherty HPA parameter, [scalar].
    % etaMax        => HPA max. efficiency, [scalar], [linear].
    % g             => HPA's power gain, [scalar], [linear].
    % K             => number of deployed IoT devices, [scalar].
    % M             => number of transmit antennas, [scalar].
    % RFNetLoss     => loss in RF network, [struct].
    % tol           => tolerance SCA [W].

    % OUTPUT:
    % HPAPow        => HPA's minimum output power, [scalar], [W].
    % bOpt          => digital beamforming, [NxQ].
    % aOpt          => phase shifter configuration, [NxM].

    % For further information, visit: https://arxiv.org/pdf/2507.06805
    %
    % This is version 1.00 (Last edited: 2025-10-20)
    %
    % License: This code is licensed under the MIT license. If you in any way
    % use this code for research that results in publications, please cite our
    % article as described above.

    %% number of antennas per RF chain (revise!!)
    L = floor(M/N);
    MTilde = N*L;
    h = h(1:MTilde,:);

    %% loss in the passive RF network
    LossPowDividers = ceil(log2(L))*RFNetLoss.Ld;
    LossPowCombiners = 1; % RFNetLoss.Lc;
    LossPhaseShifters = RFNetLoss.Lp;

    TLoss = LossPowDividers*LossPowCombiners*LossPhaseShifters;

    Q = N;

    %% Initialization procedure 
    % 3. Compute the number of RF chains per user
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

    % 4. Compute the phase shifters configuration
    idxRFChainsPerDeviceMx = unique(perms(idxRFChainsPerDeviceVec),"rows");
    totalPermutations = size(idxRFChainsPerDeviceMx,1);
    HPAPowVec = zeros(totalPermutations,1);
    bCell = cell(totalPermutations,1);
    aCell = cell(totalPermutations,1);

    for ii = 1:totalPermutations
        % phase shifter matrix construction (block diagonal)
        aTemp = zeros(MTilde,N);
        idxInit = 1;
        idxEnd = L;
        for n = 1:N
            phase = angle(h(idxInit:idxEnd,idxRFChainsPerDeviceMx(ii,n)));
            aTemp(idxInit:idxEnd,n) = exp(1i*phase);

            idxInit = idxInit + L;
            idxEnd = idxEnd + L;
        end

        [bTemp, HPAPow] = HBPCInit(N,K,L,TLoss,h,aTemp,g,p,Pmax,l,etaMax);

        HPAPowVec(ii) = HPAPow;
        bCell{ii} = bTemp;
        aCell{ii} = aTemp;
    end

    %
    [~,idx] = min(HPAPowVec);
    bInit = bCell{idx};
    aInit = aCell{idx};

    %% Solution of the main problem via SCA optimization

    % init. optimization variables
    z = zeros(K,Q);
    nu = cell(K,Q);
    for kk = 1:K
        % effective channel 1/sqrt(L*TLoss)
        hTilde = 1/sqrt(L*TLoss)*(h(:,kk)'*aInit)';
        for qq = 1:Q       
            % auxiliary variables
            z(kk,qq) = hTilde'*bInit(:,qq);
            nu{kk,qq} = z(kk,qq)*hTilde + bInit(:,qq);
        end
    end

    % SCA loop
    HPAPowPrev = 0;
    HPAPowNext = Inf;
    iter = 0;
    while (abs(HPAPowNext - HPAPowPrev) > tol)
        % display the state of the process
        iter = iter + 1;

        % cvx block
        cvx_begin quiet
        variable t(N) nonnegative
        variable b(N,Q) complex 
        variable aVec(MTilde) complex

        minimize ( sum(t) )
        subject to
            % epigraph constraints from piecewise objective
            for nn = 1:N
                if g*norm(bInit(nn,:))^2 <= Pmax/l^2
                    sqrt(g*Pmax)/(l*etaMax)*norm(b(nn,:)) <= t(nn)
                    g*square_pos(norm(b(nn,:))) <= Pmax/l^2
                else 
                    % epigraph constraint
                    (l+1)*sqrt(g*Pmax)/(l*etaMax)*norm(b(nn,:)) - Pmax/(l*etaMax) <= t(nn)

                    % box constraint
                    g*square_pos(norm(b(nn,:))) <= Pmax
                    2*g*real(conj(bInit(nn,:))*b(nn,:).') + g*norm(bInit(nn,:))^2 >= Pmax/l^2
                end
            end

            % phase shift matrix construction
            expression aMx(MTilde,N)

            idxInit = 1;
            idxEnd = L;
            for nn = 1:N
                aMx(idxInit:idxEnd,nn) = aVec(idxInit:idxEnd);  
                if idxInit == 1 
                    aMx(idxEnd+1:end,nn) = 0;
                elseif idxEnd == MTilde
                    aMx(1:idxInit-1,nn) = 0;
                else
                    aMx(1:idxInit-1,nn) = 0;
                    aMx(idxEnd+1:end,nn) = 0;
                end
                idxInit = idxInit + L;
                idxEnd = idxEnd + L;
            end

            % power requirement constraint
            for kk = 1:K
                % effective channel
                hTilde = 1/sqrt(L*TLoss)*(h(:,kk)'*aMx)';

                temp = 0;
                for qq = 1:Q
                    temp = temp + real(nu{kk,qq}'*(z(kk,qq)*hTilde + b(:,qq))) - ...
                        1/2*norm(nu{kk,qq})^2 - 1/2*square_pos(norm(z(kk,qq)*hTilde - b(:,qq))) ...
                        - abs(z(kk,qq))^2;
                end

                g*temp/p(kk) >= 1
            end

            % relaxed constant modulus constraint (analog beamforming)
            for mm = 1:MTilde
                abs(aVec(mm)) <= 1
            end

        cvx_end

        % update optimization variables

        % return the previous value before the optimization returned
        % NaN/Inf
        if isnan(cvx_optval) || isinf(cvx_optval)
            break;
        end

        bInit = b;
        aInit = aVec;

        z = zeros(K,Q);
        nu = cell(K,Q);
        for kk = 1:K
            hTilde = 1/sqrt(L*TLoss)*(h(:,kk)'*aMx)';
            for qq = 1:Q        
                z(kk,qq) = hTilde'*bInit(:,qq);
                nu{kk,qq} = z(kk,qq)*hTilde + bInit(:,qq);
            end
        end

        % update ref. previous value
        HPAPowPrev = HPAPowNext;

        % Compute the HPA's power consumption
        HPAPowNext = HPAPowConsumptFnc(N,g,b,Pmax,l,etaMax);
    end

    %% output
    fprintf('HBPC algorithm successfully converged to optimal value %.4e after %d iterations. Number of RF chains and antennas are %d and %d, respectively.\n', ...
            HPAPowNext, iter, N, M);

    HPAPow = HPAPowNext;
    bOpt = bInit;

    % phase shift matrix construction
    aOpt = zeros(MTilde,N);
    idxInit = 1;
    idxEnd = L;
    for nn = 1:N
        aOpt(idxInit:idxEnd,nn) = aInit(idxInit:idxEnd);  

        idxInit = idxInit + L;
        idxEnd = idxEnd + L;
    end
end

function [b, HPAPow] = HBPCInit(N,K,L,TLoss,h,aTemp,g,p,Pmax,l,etaMax)

    % effective channel matrix 
    HTilde = zeros(N,N,K);
    for k = 1:K
        hTilde = 1/sqrt(L*TLoss)*(h(:,k)'*aTemp)';
        HTilde(:,:,k) = hTilde*hTilde';
    end
    
    % basis vectors in rank-1 matrix form
    E = zeros(N,N,K);
    e = eye(N);
    for n = 1:N
        E(:,:,n) = e(:,n)*e(:,n)';
    end
    
    % SDP problem: minimize the HPA's input power subject to a
    % minimum devices' received power requirement.
    cvx_begin sdp quiet
    variable B(N,N) hermitian semidefinite
    variable t nonnegative
    minimize ( t )
    subject to 
    for k = 1:K
        g*real(trace(HTilde(:,:,k)*B)/p(k)) >= 1
    end
    
    for n = 1:N
        g*trace(E(:,:,n)*B) <= t
    end
    cvx_end
    
    % solution recovery via eigenvalue decomposition
    r = rank(B);
    [U,S] = eig(B);
    b = zeros(N,r);
    for rr = 1:r
        b(:,rr) = sqrt(S(rr,rr))*U(:,rr);
    end

    HPAPow = HPAPowConsumptFnc(N,g,b,Pmax,l,etaMax);
end

function HPAPow = HPAPowConsumptFnc(N,g,b,Pmax,l,etaMax)
    HPAPow = 0;
    for n = 1:N
        %
        if g*norm(b(n,:))^2 <= Pmax/l^2
            HPAPow = HPAPow + sqrt(g*norm(b(n,:))^2*Pmax)/(l*etaMax);
        else
            HPAPow = HPAPow + ((l+1)*sqrt(g*norm(b(n,:))^2*Pmax) - Pmax)/(l*etaMax);
        end
    end
end