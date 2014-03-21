function [Wsolution,transmitpower] = function_QoSproblem_singleBS(H,SINRconstraints,Qsqrt,q)
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
%This file considers the special case of only a macro base station.
%The power minimization under QoS requirements and power constraints is
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
%                  index transmit antennas
%SINRconstraints = Kr x 1 vector with SINR constraints for all users.
%Qsqrt           = N x N x L matrix with matrix-square roots of the L 
%                  weighting matrices for the L power constraints
%q               = Limits of the L power constraints
%
%OUTPUT:
%Wsolution     = N x Kr matrix with beamforming that solves the QoS problem. 
%                This matrix is empty if this problem isinfeasible.
%transmitpower = Total transmit power at the final solution


Kr = size(H,1); %Number of users
N = size(H,2); %Number of transmit antennas
L = length(q); %Number of power constraints


%Solve the power minimization under QoS requirements problem using CVX
cvx_begin
cvx_quiet(true); % This suppresses screen output from the solver

variable W(N,Kr) complex; %Variable for N x Kr beamforming matrix

minimize norm(W,'fro') %Minimize the square root of the total power

subject to

%SINR constraints (Kr constraints)
for k = 1:Kr
    %SOCP formulation for the SINR constraint of user k
    sqrt(1+SINRconstraints(k))*real(H(k,:)*W(:,k))  >= sqrt(SINRconstraints(k)) * norm([1 zeros(1,Kr-1); H(k,:)*W],'fro');
    
    imag(H(k,:)*W(:,k)) == 0;  %Useful link is assumed to be real-valued
end

%Power constraints (L constraints)
for l = 1:L
    norm(Qsqrt(:,:,l)*W,'fro') <= sqrt(q(l));
end
    
cvx_end


%Store solutions from the optimization problem
Wsolution = W;
transmitpower = norm(Wsolution,'fro').^2;

