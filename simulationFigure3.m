%This Matlab script can be used to generate Figure 3 in the article:
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
%Please note that the channels are generated randomly, thus the results
%will not be exactly the same as in the paper.

%Initialization
close all;
clear all;


%%Simulation parameters

nbrOfRealizations = 10; %Number of channel realizations

%There are extreme cases when the problem becomes infeasible due to
%insufficient power at macro BS and a poor user location realization. This
%is the result of bad scheduling (random scheduling in this code) and
%should either be removed (manually) or the power constraints at the macro
%BS can be relaxed to make the problem feasible (when this option is true).
relaxPowerConstraintsIfNeeded = true; 

nbrOfSubCarriers = 600; %Number of subcarriers in the LTE-like system
subcarrierBandwidth = 15e3; %Hz  in the LTE-like system

noiseFigure = 5; %dB
noiseFloordBm = -174+10*log10(subcarrierBandwidth)+noiseFigure; %Noise floor in dBm

Kt = 5; %Number of transmitters: One BS + Kt-1 Small Cells
Kr = 10; %Number of users

QoSconstraints = 2; %QoS constraint in bits/s/Hz per user
SINRconstraints = (2^QoSconstraints-1)*ones(Kr,1); %Corresponding SINR constraints per user

rng('shuffle'); %Initiate the random number generators with a random seed

%%If rng('shuffle'); is not supported by your Matlab version, you can use
%%the following commands instead:
%randn('state',sum(100*clock));
%rand('state',sum(100*clock));



%%Parameters for Macro Base Station (BS)

P_BS = 66; %Maximal transmit power per subcarrier at BS (in mW) ---Change to something very large to guarantee a feasible solution!
rho0_inv = 0.388; %Efficiency of power amplifiers at BS
circuitPowerBS = 189/nbrOfSubCarriers; %Circuit power per antenna (in mW)
minimalUserDistanceBS = 0.035; %kilometers

%Different alternatives of number of antennas at macro BS
NBScases = 20:10:100;

%Propagation Parameters
pentrationlossdB = 20; %Penetration loss (indoor users)
shadowFadingBS = 10^(7/10); %Standard deviation of shadow fading for all users to the macro BS



%%Parameters for Small Cell Access points (SCAs)

P_SCA = 0.08; %%Maximal transmit power per subcarrier at SCA (in mW)
rhoj_inv = 0.052; %Efficiency of power amplifiers at SCA
circuitPowerSCA = 5.6/nbrOfSubCarriers; %Circuit power per antenna (in mW)
minimalUserDistanceSCA = 0.003; %kilometers

%Different alternatives of number of antennas at SCAs
NSCAcases = [0 1 2 3];

%Propagation Parameters
clusterSize = 0.04; %Users within this distance have better propagation properties
pentrationlossdB_incluster = 0; %Pentration loss for users in cluser (dB)
pentrationlossdB_outcluster = 20; %Pentration loss for users not in cluser (dB)

%Standard deviation of shadow fading for all users to the SCA
shadowFadingSCA_incluster = 10^(7/10);
shadowFadingSCA_outcluster = 10^(7/10);



%%Preliminaries for generating simulation environment

radius = 0.5; %SCAs are located uniformly on a circle with this radius (in km)
BSlocations = [0 0; radius/2 radius/2; radius/2 -radius/2; -radius/2 -radius/2; -radius/2 radius/2]; %Position of BS and SCAs in km.
minimalUserDistances = [minimalUserDistanceBS; minimalUserDistanceSCA*ones(Kt-1,1)]; %Minimal user distance for each transmitter in km.


%Matrix for saving the total power in the optimal solutions
totalpower = zeros(length(NBScases),length(NSCAcases),nbrOfRealizations);

%Matrix for saving which antennas that serve each user at the optimal
%solutions.
optimalUserAllocation = zeros(Kt,Kr,length(NBScases),length(NSCAcases),nbrOfRealizations);


%Generate channel realization (before correlation and pathlosses)
maxNbrOfAntennas = max(NBScases)+(Kt-1)*max(NSCAcases);
Huncorr = (randn(Kr,maxNbrOfAntennas,nbrOfRealizations)+1i*randn(Kr,maxNbrOfAntennas,nbrOfRealizations))/sqrt(2);

%Generate realizations for shadow fading
shadowFadingRealizations = randn(Kt,Kr,nbrOfRealizations);


%For-loop that goes through each system realization
for iter = 1:nbrOfRealizations
    [iter nbrOfRealizations] %Simple output to keep track of progress
    
    %Step 1: Generate user locations
    
    %Uniformly distributed users in a cell with radius defined above
    distance_normalized = sqrt(rand(Kr,1)); %Uniformly distributed distances from cell center
    angles = rand(Kr,1)*2*pi; %Uniformly distributed angles
    
    MSlocations = zeros(Kr,2); %Vector to store user loaction
    
    %Go through all the transmitters
    for j = 1:Kt
        if j==1 %Kr-Kt+1 users are placed uniformly in the whole cell
            MSlocations(1:Kr-Kt+j,:) = radius*repmat(distance_normalized(1:Kr-Kt+j),[1 2]).*[cos(angles(1:Kr-Kt+j)) sin(angles(1:Kr-Kt+j))]; %Store user location
        else %1 user per SCA is placed within its cluster (defined above)
            position = clusterSize*repmat(distance_normalized(Kr-Kt+j),[1 2]) .*[cos(angles(Kr-Kt+j)) sin(angles(Kr-Kt+j))];
            
            %Check that the user is not too close to the SCA. By
            %design, this user will satisfy the minimal distance
            %requirements of the other transmitters as well.
            while norm(position) < minimalUserDistances(j)
                position = clusterSize*repmat(sqrt(rand(1)),[1 2]) .*[cos(angles(Kr-Kt+j)) sin(angles(Kr-Kt+j))]; %Create a new user position if minimal distance is not satisfied
            end
            
            MSlocations(Kr-Kt+j,:) = BSlocations(j,:) + position;  %Store user location around the SCA
        end
    end
    
    %Go through all users and check that they satisfy the minimal
    %distance requirements for all transmitters. If not, generate a new
    %random distance and angle from the origin and try again.
    for k = 1:Kr
        while min(sqrt(sum(abs(BSlocations-repmat(MSlocations(k,:),[Kt 1])).^2,2))-minimalUserDistances)<0
            distance = radius*sqrt(rand(1,1)); %Uniformly distributed distance from cell center
            angles = rand(1,1)*2*pi;  %Uniformly distributed angle
            MSlocations(k,:) = repmat(distance,[1 2]).*[cos(angles) sin(angles)]; %Store new user location
        end
    end
    


        
    
    %Go through all the different numbers of BS antennas to generate
    %channels and compute the results.
    for scenario = 1:length(NBScases)
        %[scenario length(NBScases)] %Simple output to keep track of progress
        
        Nantennas = [NBScases(scenario); max(NSCAcases)*ones(Kt-1,1)]; %Vector with number of antennas for each transmitter
        N = sum(Nantennas); %Total number of transmit antennas
        H = Huncorr(:,[1:NBScases(scenario) max(NBScases)+1:maxNbrOfAntennas],iter);
        
        antennaInds = [0; cumsum(Nantennas)]; %Vector with indicies where the antenna indices of each transmitter starts and ends.
        
        %Generate channel covariance structure for macro BS using the physical channel model in [2, Eq. (34)]
        nbrOfScatteringClusters = Nantennas(1)/2;
        uniformangles = -pi/2 + (1:nbrOfScatteringClusters)*pi/nbrOfScatteringClusters;
        omega = 0.3;
        A = [exp( - 1i*2*pi*omega* (0:Nantennas(1)-1)' *  sin(uniformangles))/sqrt(nbrOfScatteringClusters) zeros(Nantennas(1),nbrOfScatteringClusters)];
        


        
        %Compute the average channel attentuation (pathloss and shadow
        %fading) divided by the noise power and renormalize the channel
        %realization.
        lossovernoise = zeros(Kt,Kr);
        
        %Go through all users
        for k = 1:Kr
            
            dista = sqrt(sum(abs(BSlocations-repmat(MSlocations(k,:),[Kt 1])).^2,2)); %Compute distance from each transmitter
            
            %Go through all transmitters
            for j = 1:Kt
                if j==1
                    %Compute channel attenuation over noise power for the macro BS
                    lossovernoise(j,k) = 10^( (128.1+37.6*log10(dista(j)) + shadowFadingBS*shadowFadingRealizations(j,k,iter) + noiseFloordBm + pentrationlossdB)/10);
                    
                    %Apply this to the current channel realization and
                    %multiply with A that gives the spatial correlation.
                    %(Note that we increase the channel gain by rho0_inv,
                    %which will be neutralized later by reducing the power
                    %by 1/rho0_inv.)
                    H(k,1+antennaInds(j):antennaInds(j+1)) = (H(k,1+antennaInds(j):antennaInds(j+1))*A')*sqrt(rho0_inv/lossovernoise(j,k));
                else
                    %Compute channel attenuation over noise power for the
                    %SCAs, depending on whether the user is within 40 m or
                    %not.
                    if dista(j)<=clusterSize
                        lossovernoise(j,k) = 10^( (127+30*log10(dista(j)) + shadowFadingSCA_incluster*shadowFadingRealizations(j,k,iter) + noiseFloordBm + pentrationlossdB_incluster)/10);
                    else
                        lossovernoise(j,k) = 10^( (128.1+37.6*log10(dista(j)) + shadowFadingSCA_outcluster*shadowFadingRealizations(j,k,iter) + noiseFloordBm + pentrationlossdB_outcluster)/10);
                    end
                    
                    %Apply this to the current channel realization.
                    %(Note that we increase the channel gain by rhoj_inv,
                    %which will be neutralized later by reducing the power
                    %by 1/rhoj_inv.)
                    H(k,1+antennaInds(j):antennaInds(j+1)) = H(k,1+antennaInds(j):antennaInds(j+1))*sqrt(rhoj_inv/lossovernoise(j,k));
                end
            end
            
        end
        
        
        %Create the power constraints in vector using the notation from the
        %book "Optimal Resource Allocation in Coordinated Multi-Cell Systems"
        %by Emil Björnson and Eduard Jorswieck.
        
        L = N; %Number of constraints (one per antenna)

        Q = zeros(N,N,L); %Weighting matrices
        Qsqrt = zeros(N,N,L); %Matrix square root of weighting matrices
  
        %Create per-antenna constraints
        for j = 1:L
            Q(j,j,j) = 1;
            Qsqrt(j,j,j) = 1;
        end
        
        %Store the maximum power that correspond to each power constraint
        %(Note that the power normalized by divided by rho0_inv and
        %rhoj_inv, to avoid feeding this parameter into the optimization
        %solvers. The inverse normalization was done for the channels
        %earlier in the code.)
        q = [(P_BS/rho0_inv)*ones(Nantennas(1),1); (P_SCA/rhoj_inv)*ones(max(NSCAcases)*(Kt-1),1)];
        
        
        %Compute the minimal transmit power for each of the different
        %number of antennas at the SCAs.
        for m = 1:length(NSCAcases)
            
            %This is the case with only the macro BS
            if m==1
                %Compute the optimal solution by solving a convex
                %optimization problem
                [Wsolution_BSonly,transmitpower_BSonly] = function_QoSproblem_singleBS(H(:,1:Nantennas(1)),SINRconstraints,Qsqrt(1:Nantennas(1),1:Nantennas(1),1:Nantennas(1)),q(1:Nantennas(1)));
                
                %Takes care of extreme cases when the problem is infeasible
                if relaxPowerConstraintsIfNeeded && isnan(transmitpower_BSonly)
                    [Wsolution_BSonly,transmitpower_BSonly] = function_QoSproblem_singleBS(H(:,1:Nantennas(1)),SINRconstraints,[],[]);
                end
                
                %Compute the total power using the definition in Eq. (4)
                %and Eq. (5)
                totalpower(scenario,m,iter) = Nantennas(1)*circuitPowerBS + transmitpower_BSonly;
                
                optimalUserAllocation(1,:,scenario,m,iter) = 1; %Store which antennas that serve the user in the final solution
            else
                
                %Find the current channel matrix by simply removing the
                %columns for SCA-antennas that are not in use.
                Hcurrent = zeros(Kr,NBScases(scenario)+NSCAcases(m)*(Kt-1));
                Hcurrent(:,1:NBScases(scenario)) = H(:,1:NBScases(scenario));
                for j = 2:Kt
                    Hcurrent(:,NBScases(scenario)+(j-2)*NSCAcases(m)+1:NBScases(scenario)+(j-1)*NSCAcases(m)) = H(:,NBScases(scenario)+(j-2)*max(NSCAcases)+1:NBScases(scenario)+(j-2)*max(NSCAcases)+NSCAcases(m));
                end
                
                %Update vectors with current number of antennas
                Nantennas_current = [NBScases(scenario); NSCAcases(m)*ones(Kt-1,1)];
                Ncurrent = sum(Nantennas_current);
                
                %Compute the optimal solution by solving a convex
                %optimization problem
                [Wsolution,transmitpower_BS_SCA,BSallocation] = function_QoSproblem_relaxation(Hcurrent,Nantennas_current,SINRconstraints,Q(1:Ncurrent,1:Ncurrent,1:Ncurrent),q(1:Ncurrent));
                
                %Takes care of extreme cases when the problem is infeasible
                if relaxPowerConstraintsIfNeeded && isnan(transmitpower_BS_SCA)
                    [Wsolution,transmitpower_BS_SCA,BSallocation] = function_QoSproblem_relaxation(Hcurrent,Nantennas_current,SINRconstraints,Q(1:Ncurrent,1:Ncurrent,NBScases(scenario)+1:Ncurrent),q(NBScases(scenario)+1:Ncurrent));
                end
                
                %Compute the total power using the definition in Eq. (4)
                %and Eq. (5)
                totalpower(scenario,m,iter) = Nantennas(1)*circuitPowerBS + (Kt-1)*NSCAcases(m)*circuitPowerSCA + transmitpower_BS_SCA(1);
                
                optimalUserAllocation(:,:,scenario,m,iter) = BSallocation; %Store which antennas that serve the user in the final solution
            end
            
        end
        
        
    end
    
end



%Compute the average over channel realizations (and remove infeasiblities and outliers)
feasible = sum(~isnan(totalpower) & totalpower<1e6,3);
totalpower(isnan(totalpower) | totalpower>=1e6) = 0;
totalpower_average = sum(totalpower,3)./feasible;


%Plot the numerical results
figure; hold on; box on;

plot(NBScases,10*log10(totalpower_average(:,1)),'kd:','LineWidth',1);
plot(NBScases,10*log10(totalpower_average(:,2)),'bo-.','LineWidth',1);
plot(NBScases,10*log10(totalpower_average(:,3)),'kd--','LineWidth',1);
plot(NBScases,10*log10(totalpower_average(:,4)),'rh-','LineWidth',1);

legend('No SCAs','1 antenna/SCA','2 antenna/SCA','3 antenna/SCA','Location','SouthWest');

xlabel('Antennas at the BS, N_{BS}')
ylabel('Total Power per Subcarrier [dBm]')

