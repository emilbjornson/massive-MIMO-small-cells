function [Wsolution,totalpower,BSallocation] = function_QoSproblem_multiflowRZF(H,Nantennas,SINRconstraints,Q,q)
%Implementation of the Multiflow-RZF algorithm from the article:
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
%multiple small cells. Each transmitter applies regularized zero-forcing
%locally and perform joint power allocation to solve the following problem:
%
%minimize   total transmit power
%subject to SINR_k >= SINRconstraints(k) for all users k,
%           Power constraints.
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

g = zeros(Kt,Kr,Kr); %Placeholder for effective channel gain using RZF
QQ = zeros(Kt,Kr,L); %Placeholder for effective power constraints using RZF: QQ = u^H Q u



Wsolution=zeros(N,Kr); %Placeholder for beamforming solution

%Vector with indicies where the antenna indices of each transmitter starts and ends.
antennaInds=[0; cumsum(Nantennas)];

%Step 1: Go through all transmitters and all users. Transform the original
%problem into a simplified version using RZF
for j=1:Kt
    
    for i=1:Kr
        %Compute beamforming direction using RZF
        A = H([1:i-1 i+1:Kr],1+antennaInds(j):antennaInds(j+1))' * H([1:i-1 i+1:Kr],1+antennaInds(j):antennaInds(j+1))+diag(Kr./q(1+antennaInds(j):antennaInds(j+1))/max(SINRconstraints));
        w = A\H(i,1+antennaInds(j):antennaInds(j+1))';
        w = w/norm(w);
        Wsolution(1+antennaInds(j):antennaInds(j+1),i) = w;
        
        %Compute the corresponding channel gain
        for k = 1:Kr
            g(j,i,k) = abs(H(k,1+antennaInds(j):antennaInds(j+1))*w).^2;
        end
        
        %Transform the power constraints into the equivalent form in (9)
        for l = 1:L
            QQ(j,i,l) = real(w'*Q(1+antennaInds(j):antennaInds(j+1),1+antennaInds(j):antennaInds(j+1),l)*w);
        end
    end
end


%Pre-compute a block-diagonal matrix with the channels from each
%transmitter.
Rh = zeros(Kt,N,Kr);
for j = 1:Kt
    for k = 1:Kr
        Rh(j,1+antennaInds(j):antennaInds(j+1),k) = H(k,1+antennaInds(j):antennaInds(j+1));
    end
end


%Solve the power minimization under QoS requirements problem using CVX
cvx_begin
cvx_quiet(true); % this suppresses screen output from the solver
variable P(Kt,Kr); %Variable for Kt x Kr power allocation matrix

minimize real(sum(P(:))) %Minimize the total power

subject to

P >= 0; %Power must be non-negative

%SINR constraints (Kr constraints)
for k = 1:Kr
    g(:,k,k)'*P(:,k)*(1+SINRconstraints(k)) >= SINRconstraints(k)*(1 + sum(sum(g(:,:,k).*P)) );
end

%Power constraints (L constraints)
for l = 1:L
    sum(sum(QQ(:,:,l).*P)) <= q(l);
end

cvx_end


%Store solutions from the optimization problem and compute user allocation
totalpower = real(sum(P(:)));

BSallocation = zeros(Kt,Kr);

for j = 1:Kt
    for i = 1:Kr
        Wsolution(1+antennaInds(j):antennaInds(j+1),i) = Wsolution(1+antennaInds(j):antennaInds(j+1),i)*sqrt(P(j,i));
        
        BSallocation(j,i) = g(j,i,i)'*P(j,i)>1e-4;
    end
end



