function [HPAPow,bOpt] = fullyDigital(h,M,p,Pmax,l,etaMax,g,K,tol)
    % This function minimizes the power consumption of a PB equipped with a
    % hybrid analog-digital beamforming architecture. 

    % INPUT
    % h             => channel coefficients, [NxK].
    % N             => number of RF chains, [scalar].
    % Q             => number of streams, [scalar].
    % p             => received power requirements, [Kx1], [W].
    % Pmax          => HPA max. output power, [scalar], [W].
    % l             => l-way Doherty HPA parameter, [scalar].
    % etaMax        => HPA max. efficiency, [scalar], [linear].
    % g             => HPA's power gain, [scalar], [linear].
    % K             => number of deployed IoT devices, [scalar].
    % tol           => tolerance SCA [W].

    % OUTPUT:
    % HPAPow        => HPA's minimum output power, [scalar], [W].
    % bOpt          => digital beamforming, [MxQ].

    % For further information, visit: https://arxiv.org/pdf/2507.06805
    %
    % This is version 1.00 (Last edited: 2025-10-20)
    %
    % License: This code is licensed under the MIT license. If you in any way
    % use this code for research that results in publications, please cite our
    % article as described above.
    
    Q = M;

    %% Initialization via SDP formulation

    % channel matrix (SDP formulation)
    H = zeros(M,M,K);
    for k = 1:K
        H(:,:,k) = h(:,k)*h(:,k)';
    end

    % Basis vectors in rank-1 matrix form (SDP formulation)
    E = zeros(M,M,K);
    e = eye(M);
    for m = 1:M
        E(:,:,m) = e(:,m)*e(:,m)';
    end

    % SDP problem
    cvx_begin sdp quiet
        variable B(M,M) hermitian semidefinite
        variable t nonnegative
        minimize ( t )
        subject to 
            for k = 1:K
                g*real(trace(H(:,:,k)*B)/p(k)) >= 1
            end

            for m = 1:M
                g*trace(E(:,:,m)*B) <= t
            end
    cvx_end

    % Solution recovery via eigenvalue decomposition
    r = rank(B);    
    [U,S] = eig(B);
    bTemp = zeros(M,r);
    for rr = 1:r
        bTemp(:,rr) = sqrt(S(rr,rr))*U(:,rr);
    end

    %% Solution via SCA optimization

    % init. optimization variables
    z = zeros(K,Q);
    for k = 1:K
        for q = 1:Q
            z(k,q) = h(:,k)'*bTemp(:,q);
        end
    end

    % SCA loop 
    HPAPowPrev = 0;
    HPAPowNext = Inf;
    iter = 0;
    while (abs(HPAPowNext - HPAPowPrev) > tol)
        % update ref. previous value
        HPAPowPrev = HPAPowNext; 
        iter = iter + 1;

        cvx_begin quiet
        variable t(M) nonnegative
        variable b(M,Q) complex
        minimize ( sum(t) )
        subject to
            % epigraph constraints from piecewise objective
            for m = 1:M
                if g*norm(bTemp(m,:))^2 <= Pmax/l^2
                    sqrt(g*Pmax)/(l*etaMax)*norm(b(m,:)) <= t(m)
                    g*norm(b(m,:)) <= Pmax/l^2
                else 
                    (l+1)*sqrt(g*Pmax)/(l*etaMax)*norm(b(m,:)) - Pmax/(l*etaMax) <= t(m)

                    g*norm(b(m,:)) <= Pmax
                    2*g*real(conj(bTemp(m,:))*b(m,:).') + g*norm(bTemp(m,:))^2 >= Pmax/l^2
                end
            end

            % power requirement constraint
            for k = 1:K                
                temp = 0;
                for q = 1:Q
                    temp = temp + 2*real(z(k,q)'*h(:,k)'*b(:,q)) - abs(z(k,q))^2;
                end

                g*temp/p(k) >= 1
            end

        cvx_end

        % update optimization variables
        bTemp = b;

        for k = 1:K
            for q = 1:Q
                z(k,q) = h(:,k)'*b(:,q);
            end
        end

        % Compute the HPA's power consumption
        HPAPowNext = 0;
        for m = 1:M
            %
            if g*norm(b(m,:))^2 <= Pmax/l^2
                HPAPowNext = HPAPowNext + sqrt(g*norm(b(m,:))^2*Pmax)/(l*etaMax);
            else
                HPAPowNext = HPAPowNext + ((l+1)*sqrt(g*norm(b(m,:))^2*Pmax) - Pmax)/(l*etaMax);
            end
        end
    end

    %% output
    fprintf('FD algorithm successfully converged to optimal value %.4e after %d iterations. Number antennas is %d. \n', ...
            cvx_optval, iter, M);

    HPAPow = HPAPowNext;
    bOpt = b;
end