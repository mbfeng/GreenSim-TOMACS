% File2_SimPaths.m
% Simulate sample paths and calculate true prices

%% Simulate sample paths
Paths = zeros(NPeriods*NPaths, 2);      % placeholder for sample paths
PriceTrue = zeros(NPeriods*NPaths, 1);  % placeholder for true prices

for i = 1:NPeriods*NPaths
    if(mod(i, NPeriods) == 1)
        % every paths starts from stationary distribution
        Paths(i,:) = normrnd(station_mean, station_std);
    else
        % sample paths evolves as an AR(1) process
        Paths(i,:) = mu + eta.*Paths(i-1,:) + normrnd([0, 0], sigma);
    end
end

% transformation from state space to parameter space
Paths(:,1) = sigmoid(lambda_lb, lambda_ub, Paths(:,1));
Paths(:,2) = sigmoid(theta_lb, theta_ub, Paths(:,2));

%% Calculate true prices in each period in each sample path
NTerms = 40; % number of approximation terms in calculating true prices
for i = 1:NPeriods*NPaths
    PriceTrue(i) = CatBondPrice(Paths(i,1), Paths(i,2), K, p, NTerms);
end

PriceTrue = FaceVal.*PriceTrue;
Price_lo = FaceVal.*CatBondPrice(lambda_lo, theta_lo, K, p, NTerms);
Price_mi = FaceVal.*CatBondPrice(lambda_mi, theta_mi, K, p, NTerms);
Price_hi = FaceVal.*CatBondPrice(lambda_hi, theta_hi, K, p, NTerms);