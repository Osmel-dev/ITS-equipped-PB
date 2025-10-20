function [HPAPow,bOpt,aOpt] = HBFC(h,N,p,Pmax,l,etaMax,g,K,M,tol,RFNetLoss)
    % This function minimizes the power consumption of a PB equipped with a
    % fully-connected hybrid analog-digital beamforming architecture. 

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

    % loss in the passive RF network
    LossPowDividers = ceil(log2(M))*RFNetLoss.Ld;
    LossPowCombiners = ceil(log2(N))*RFNetLoss.Lc;
    LossPhaseShifters = RFNetLoss.Lp;

    TLoss = LossPowDividers*LossPowCombiners*LossPhaseShifters;

    Q = N;

     %% Initialization procedure 

    % Compute the number of RF chains per user
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
    HPAPowVec = zeros(totalPermutations,1);
    bCell = cell(totalPermutations,1);
    aCell = cell(totalPermutations,1);

    for ii = 1:totalPermutations
        aTemp = zeros(M,N);
        for n = 1:numel(idxRFChainsPerDeviceVec)
            aTemp(:,n) = exp(1i*angle(h(:,idxRFChainsPerDeviceMx(ii,n))));
        end 

        [bTemp, HPAPow] = HBFCInit(h,N,p,Pmax,l,etaMax,g,K,M,TLoss,aTemp);

        % HPAPow
        HPAPowVec(ii) = HPAPow;
        bCell{ii} = bTemp;
        aCell{ii} = aTemp;
    end

    %
    [~,idx] = min(HPAPowVec);
    bInit = bCell{idx};
    aInit = aCell{idx};

    %% Solution via SCA optimization

    % init. optimization variables
    z = zeros(K,Q);
    nu = cell(K,Q);
    for k = 1:K
        % previosuly 1/sqrt(M*TLoss)
        hTilde = 1/sqrt(M*TLoss)*(h(:,k)'*aInit)';
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
        variable b(N,Q) complex 
        variable a(M,N) complex

        minimize ( sum(t) ) % - lambda*(2*real(trace(aTemp'*a)) + norm(aTemp,'fro')^2)
        subject to
            % epigraph constraints from piecewise objective
            for n = 1:N
                if g*norm(bInit(n,:))^2 <= Pmax/l^2
                    sqrt(g*Pmax)/(l*etaMax)*norm(b(n,:)) <= t(n)
                    g*square_pos(norm(b(n,:))) <= Pmax/l^2
                else 
                    % epigraph constraint
                    (l+1)*sqrt(g*Pmax)/(l*etaMax)*norm(b(n,:)) - Pmax/(l*etaMax) <= t(n)

                    % box constraint
                    g*square_pos(norm(b(n,:))) <= Pmax
                    2*g*real(conj(bInit(n,:))*b(n,:).') + g*norm(bInit(n,:))^2 >= Pmax/l^2
                end
            end

            % power requirement constraint
            for k = 1:K
                % "effective channel"
                hTilde = 1/sqrt(M*TLoss)*(h(:,k)'*a)';

                temp = 0;
                for q = 1:Q
                    temp = temp + real(nu{k,q}'*(z(k,q)*hTilde + b(:,q))) ...
                        - 1/2*norm(nu{k,q})^2 - 1/2*square_pos(norm(z(k,q)*hTilde - b(:,q))) ...
                        - abs(z(k,q))^2;
                end

                g*temp/p(k) >= 1
            end

            % relaxed constant modulus constraint (analog beamforming)
            for m = 1:M
                for n = 1:N
                    abs(a(m,n)) <= 1;
                end
            end

        cvx_end

        % update optimization variables

        % return the previous value before the optimization returned
        % NaN/Inf
        if isnan(cvx_optval) || isinf(cvx_optval)
            break;
        end

        bInit = b;
        aInit = a;

        for k = 1:K
            hTilde = 1/sqrt(M*TLoss)*(h(:,k)'*a)';
            for q = 1:Q        
                z(k,q) = hTilde'*bInit(:,q);
                nu{k,q} = z(k,q)*hTilde + bInit(:,q);
            end
        end

        % update ref. previous value
        HPAPowPrev = HPAPowNext;

        % Compute the HPA's power consumption
        HPAPowNext = HPAPowConsumptFnc(N,Pmax,l,etaMax,g,b);
    end

    %% output
    fprintf('HBFC algorithm successfully converged to optimal value %.4e after %d iterations. Number of RF chains and antennas are %d and %d, respectively.\n', ...
            HPAPowNext, iter, N, M);

    HPAPow = HPAPowNext;
    aOpt = aInit;
    bOpt = bInit;
end

function [b, HPAPow] = HBFCInit(h,N,p,Pmax,l,etaMax,g,K,M,TLoss,aTemp)

    % channel matrix (SDP formulation)
    HTilde = zeros(N,N,K);
    for k = 1:K
        hTilde = 1/sqrt(M*TLoss)*(h(:,k)'*aTemp)';
        HTilde(:,:,k) = hTilde*hTilde';
    end
    
    % Basis vectors in rank-1 matrix form (SDP formulation)
    E = zeros(N,N,K);
    e = eye(N);
    for n = 1:N
        E(:,:,n) = e(:,n)*e(:,n)';
    end
    
    % SDP problem 
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
    
    % Solution recovery via eigenvalue decomposition
    r = rank(B);
    [U,S] = eig(B);
    b = zeros(N,r);
    for rr = 1:r
        b(:,rr) = sqrt(S(rr,rr))*U(:,rr);
    end

    HPAPow = HPAPowConsumptFnc(N,Pmax,l,etaMax,g,b);
end

function HPAPow = HPAPowConsumptFnc(N,Pmax,l,etaMax,g,b)
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