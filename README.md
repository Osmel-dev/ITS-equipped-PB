Wireless Energy Transfer Beamforming Optimization for Intelligent Transmitting Surface
==================

This code package is related to the following scientific article:

O. Martínez Rosabal, O. L. Alcaraz López, V. D. Pegorara Souto, R. D. Souza, S. Montejo-Sánchez, R. Schober, and H. Alves, "Wireless Energy Transfer Beamforming Optimization for Intelligent Transmitting Surface," accepted for publication in IEEE Transactions on Wireless Communications.

Available at: https://arxiv.org/pdf/2507.06805

## Abstract of Article

Radio frequency (RF) wireless energy transfer (WET) is a promising technology for powering the growing ecosystem of Internet of Things (IoT) using power beacons (PBs). Recent research focuses on efficient PB architectures that can support numerous antennas. In this context, PBs equipped with intelligent surfaces present a promising approach, enabling physically large, reconfigurable arrays. Motivated by these advantages, this work aims to minimize the power consumption of a PB equipped with a passive intelligent transmitting surface (ITS) and a collocated digital beamforming-based feeder to charge multiple single-antenna devices. To model the PB's power consumption accurately, we consider power amplifiers nonlinearities, ITS control power, and feeder-to-ITS air interface losses. The resulting optimization problem is highly nonlinear and nonconvex due to the high-power amplifier (HPA), the received power constraints at the devices, and the unit-modulus constraint imposed by the phase shifter configuration of the ITS. To tackle this issue, we apply successive convex approximation (SCA) to iteratively solve convex subproblems that jointly optimize the digital precoder and phase configuration. Given SCA's sensitivity to initialization, we propose an algorithm that ensures initialization feasibility while balancing convergence speed and solution quality. We compare the proposed ITS-equipped PB's power consumption against benchmark architectures featuring digital and hybrid analog-digital beamforming. Results demonstrate that the proposed architecture efficiently scales with the number of RF chains and ITS elements. We also show that nonuniform ITS power distribution influences beamforming and can shift a device between near- and far-field regions, even with a constant aperture.

## Content of Code Package

The repository contains Matlab scripts and user-defined functions required to reproduce the numerical results in the article. To run the code, you need to install the latest version of the modeling framework CVX available at https://cvxr.com/cvx/download/. 

See each file for further documentation.

## Acknowledgements

This work is partially supported by UPRISING  (Grant Number 348515) and 6G Flagship (Grant Number 369116) funded by the Research Council of Finland; European Commission through the Horizon Europe/JU SNS projects Hexa-X-II (Grant no. 101095759) and AMBIENT-6G (Grant 101192113); CNPq (305021/2021-4), RNP/MCTI Brazil 6G (01245.020548/2021-07); FAPEMIG under projects No APQ-05305-23, APQ-04523-23, PPE-00124-23, RED-00194-23, and APQ-02782-25; CNPQ under project No 443974/2024-1; ANID FONDECYT Regular No. 1241977; and by Deutsche Forschungsgemeinschaft (DFG, German Research Foundation) under projects SCHO 831/19-1 (Project-ID 556167353) and SFB 1483 (Project-ID 442419336, EmpkinS). The authors also wish to acknowledge CSC - IT Center for Science, Finland, for computational resources.

## License and Referencing

This code package is licensed under the MIT license. If you in any way use this code for research that results in publications, please cite our original article listed above.
