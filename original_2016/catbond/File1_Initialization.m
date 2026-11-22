% File1_Initialization.m
% Initilizes workspace and parameters

%% Initialize Workspace
clear;      % clear variables
close all;  % close all figures
clc;        % clear screen

rng(0);     % setup random number generator

%% Initialize Parameters 
% ----- CAT bond parameters
K = 16;     % triggering point of Cat bond 
p = 0.5;    % recovery rate in Cat bond    
FaceVal = 100; % face value of Cat bond

% ----- Experiment parameters
NPaths = 1e4;   % number of sample paths
NPeriods = 1e2;   % number of periods per sample path
NReps = 1e2;    % number of replications per period
SKStart = 5;    % starting period of stochastic kriging

% ----- State process parameters
% xk+1 = mu + eta * xk + eps, var[eps] = sigma^2
% AR(1) process <===> Euler discretization of OU-process
mu = [0.0, 0.0]; % constant long-run mean
eta = [0.6, 0.5]; % autoregressive parameter
sigma = [0.8, 0.5]; % standard deviation of noise term
station_mean = mu./(1 - eta);   % stationary mean
station_std = sqrt((sigma.^2)./(1 - (eta).^2));   % stationary std

% ----- State transformation
% bounds for Poisson mean (lambda) & exponential mean (theta)
lambda_ub = 4; lambda_lb = 2;
theta_ub = 6; theta_lb = 4;
% sigmoid function with specified bounds
sigmoid = @(lb, ub, x) lb + (ub-lb)./(1 + exp(x));

% ----- Corner and certer parameters in the numerical experiment
lambda_lo = lambda_lb; theta_lo = theta_lb;    % bottom-left corner of parameter space
lambda_hi = lambda_ub; theta_hi = theta_ub;    % top-right corner of parameter space

% center of parameter space
lambda_mi = 0.5 * (lambda_lo + lambda_hi); 
theta_mi = 0.5 * (theta_lo + theta_hi);    
