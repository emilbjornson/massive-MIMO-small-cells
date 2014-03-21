function [Wsolution,transmitpower,BSallocation] = function_QoSproblem_relaxation(H,Nantennas,SINRconstraints,Q,q)
%Solves problem "minimize the total power consumption while satisfying
%the QoS constraints and the power constraints" in Eq. (7) of the article:
%
%Emil Björnson, Marios Kountouris, Merouane Debbah, "Massive MIMO and
%Small Cells: Improving Energy Efficiency by Optimal Soft-Cell
%Coordination," Proceedings International Conference on Telecommunications
%(ICT'13), Casablanca, Morocco, May 2013.
%
%This is version 1.3 (Last edited: 2014-03-21)
%
%%The implementation utilizes and requires CVX: http://cvxr.com/
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%
%This file considers the general case with a macro base station and
%multiple small cells. The power minimization under QoS requirements and
%power constraints is
%
%minimize   total transmit power
%subject to SINR_k >= SINRconstraints(k) for all users k,
%           Power constraints.
%
%This optimization problem is convex. The computational complexity is
%therefore polynomial in the number of users, antennas, and power
%constraints. The implementation can, at least, handle 10 users, 100
%antennas, and 100 power constraints.
%
%
%INPUT:
%H               = Kr x N matrix with row index for receiver and column
%                  index transmit antennas. N is the total number of
%                  antennas and Kr is the number of users.
%Nantennas       = Kt x 1 vector with number of antennas at each
%                  transmitter
%SINRconstraints = Kr x 1 vector with SINR constraints for all users.
%Q               = N x N x L matrix with the L weighting matrices for the
%                  L power constraints
%q               = Limits of the L power constraints
%
%OUTPUT:
%Wsolution     = N x Kr matrix with beamforming that solves the QoS problem. 
%                This matrix is empty if this problem isinfeasible.
%transmitpower = Total transmit power at the final solution
%BSallocation  = N x Kr matrix with ones for all antennas that serve a
%                given user and zeros otherwise.



Kt = length(Nantennas); %Number of transmitters
Kr = size(H,1);  %Number of users
N = size(H,2); %Number of transmit antennas (in total)
L = length(q); %Number of power constraints


%Vector with indicies where the antenna indices of each transmitter starts and ends.
antennaInds=[0; cumsum(Nantennas)];

%Pre-compute the products of channel matrices for each transmitter
HH = zeros(N,N,Kr);
for j = 1:length(Nantennas)
    for k = 1:Kr
        HH(1+antennaInds(j):antennaInds(j+1),1+antennaInds(j):antennaInds(j+1),k) = H(k,1+antennaInds(j):antennaInds(j+1))' * H(k,1+antennaInds(j):antennaInds(j+1));
    end
end

%Pre-compute a block-diagonal matrix with the channels from each
%transmitter.
Rh = zeros(Kt,N,Kr);
for j = 1:Kt
    for k = 1:Kr
        Rh(j,1+antennaInds(j):antennaInds(j+1),k) =  H(k,1+antennaInds(j):antennaInds(j+1));
    end
end


%Solve the power minimization under QoS requirements problem using CVX
cvx_begin
cvx_quiet(true); % This suppresses screen output from the solver
variable W(N,N,Kr) complex; %Variable for N x N x Kr beamforming matrix

minimize real(trace(sum(W,3))) %Minimize the total power

subject to

for k = 1:Kr
    %Semi-definite formulation for the SINR constraint of user k
    real(trace(W(:,:,k)*HH(:,:,k)))*(1+SINRconstraints(k))  >= SINRconstraints(k)*(1 + real(trace(sum(W,3)*HH(:,:,k))));
    
    %Define the N x N beamforming matrix as Hermitian pos. semi-definite
    W(:,:,k) == hermitian_semidefinite(N);
end

%Power constraints (L constraints)
for l = 1:L
    real(trace(Q(:,:,l)*sum(W,3))) <= q(l);
end

cvx_end


%Store solutions from the optimization problem
Wsolution = W; %Final beamforming solution
transmitpower = trace(sum(Wsolution,3)); %Total transmit power


%Find which transmitter that allocates non-zero power (more than 1e-4) to
%the each user. This is the one that serves the user.
BSallocation = zeros(length(Nantennas),Kr);
for k = 1:Kr
    signalpower = zeros(Kt,1);
    for j = 1:Kt
        signalpower(j) = trace(HH(1+antennaInds(j):antennaInds(j+1),1+antennaInds(j):antennaInds(j+1),k)*W(1+antennaInds(j):antennaInds(j+1),1+antennaInds(j):antennaInds(j+1),k));
    end
    BSallocation(:,k) = signalpower > 1e-4;
end

