function posAnt = positionArray(M,dxy)
    % This function computes the position of the antenna/passive elements
    % within the array. 

    % INPUT:
    % M             => num. antennas/passive elements [scalar]
    % dxy           => inter-element separation [scalar]

    % OUTPUT:
    % posAnt        => positions of the antennas/passive elements [3xM]

    % For further information, visit: https://arxiv.org/pdf/2507.06805
    %
    % This is version 1.00 (Last edited: 2025-10-20)
    %
    % License: This code is licensed under the MIT license. If you in any way
    % use this code for research that results in publications, please cite our
    % article as described above.

    %% Compute the side size of the antenna array/ITS
    l = (sqrt(M) - 1)*dxy;

    %% position of the IRS' elements (centered in the x-y plane)
    posAnt = zeros(3,M);
    for m = 1:M
        posAnt(1,m) = -l/2 + mod(m-1,sqrt(M))*dxy;          % x coordinate
        posAnt(2,m) = l/2 - floor((m-1)/sqrt(M))*dxy;       % y coordinate
        posAnt(3,m) = 0;                                     % z coordinate
    end
    
end

