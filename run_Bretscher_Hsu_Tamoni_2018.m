%This small script illustrates
% Bretscher, Gsu, Tamoni (2018)
clearvars; clearvars -global; close all; clc;
addpath('utils')

%% common settings
order_app = 3;
randn('seed',1)
num_sim = 2000;
shocks = randn(2,num_sim);
MODELS = ["rbc_sv_direct", "rbc_sv_lagged_direct", "rbc_sv_FGR", "rbc_sv_FL"];
Data.label = ["ln(c)" "ln(k)" "ln(a)" "ln(siga)"];

%% "direct" specification
if ismember("rbc_sv_direct",MODELS)
    [X_sim.rbc_sv_direct, Y_sim.rbc_sv_direct] = sgu_run('rbc_sv_direct',order_app,shocks);
    Data.rbc_sv_direct = zeros(4,num_sim);
    Data.rbc_sv_direct(1,1:num_sim) = Y_sim.rbc_sv_direct(1,1:num_sim); % ln_c
    Data.rbc_sv_direct(2,1:num_sim) = X_sim.rbc_sv_direct(1,1:num_sim); % ln_k
    Data.rbc_sv_direct(3,1:num_sim) = Y_sim.rbc_sv_direct(2,1:num_sim); % ln_a
    Data.rbc_sv_direct(4,1:num_sim) = X_sim.rbc_sv_direct(3,1:num_sim); % ln_siga
end

%% "Lagged Direct" specification
if ismember("rbc_sv_lagged_direct",MODELS)
    [X_sim.rbc_sv_lagged_direct, Y_sim.rbc_sv_lagged_direct] = sgu_run('rbc_sv_lagged_direct',order_app,shocks);
    Data.rbc_sv_lagged_direct = zeros(4,num_sim);
    Data.rbc_sv_lagged_direct(1,1:num_sim) = Y_sim.rbc_sv_lagged_direct(1,1:num_sim); % ln_c
    Data.rbc_sv_lagged_direct(2,1:num_sim) = X_sim.rbc_sv_lagged_direct(1,1:num_sim); % ln_k
    Data.rbc_sv_lagged_direct(3,1:num_sim) = Y_sim.rbc_sv_lagged_direct(2,1:num_sim); % ln_a
    Data.rbc_sv_lagged_direct(4,1:num_sim) = Y_sim.rbc_sv_lagged_direct(4,1:num_sim); % ln_siga
end

%% "FGR" specification
if ismember("rbc_sv_FGR",MODELS)
    [X_sim.rbc_sv_FGR, Y_sim.rbc_sv_FGR] = sgu_run('rbc_sv_FGR',order_app,shocks);
    Data.rbc_sv_FGR = zeros(4,num_sim);
    Data.rbc_sv_FGR(1,1:num_sim) = Y_sim.rbc_sv_FGR(1,1:num_sim); % ln_c
    Data.rbc_sv_FGR(2,1:num_sim) = X_sim.rbc_sv_FGR(1,1:num_sim); % ln_k
    Data.rbc_sv_FGR(3,1:num_sim) = Y_sim.rbc_sv_FGR(2,1:num_sim); % ln_a
    Data.rbc_sv_FGR(4,1:num_sim) = Y_sim.rbc_sv_FGR(3,1:num_sim); % ln_siga
end

%% "FL" specification
if ismember("rbc_sv_FL",MODELS)
    [X_sim.rbc_sv_FL, Y_sim.rbc_sv_FL] = sgu_run('rbc_sv_FL',order_app,shocks);
    Data.rbc_sv_FL = zeros(4,num_sim);
    Data.rbc_sv_FL(1,1:num_sim) = Y_sim.rbc_sv_FL(1,1:num_sim); % ln_c
    Data.rbc_sv_FL(2,1:num_sim) = X_sim.rbc_sv_FL(1,1:num_sim); % ln_k
    Data.rbc_sv_FL(3,1:num_sim) = Y_sim.rbc_sv_FL(2,1:num_sim); % ln_a
    Data.rbc_sv_FL(4,1:num_sim) = X_sim.rbc_sv_FL(3,1:num_sim); % ln_siga
end

%% Comparison
clc
for m = MODELS
    TBL.(m) = zeros(size(Data.label,2),4);
    for j = 1:size(Data.label,2)
        TBL.(m)(j,1) = mean(Data.(m)(j,:)');
        TBL.(m)(j,2) = std(Data.(m)(j,:)');
        TBL.(m)(j,3) = skewness(Data.(m)(j,:)');
        TBL.(m)(j,4) = kurtosis(Data.(m)(j,:)');
    end
    fprintf('SPECIFICATION: %s\n',upper(m));
    disp(array2table(TBL.(m),'RowNames',Data.label,'VariableNames',{'MEAN','STD','SKEWNESS','KURTOSIS'}));
end

rmpath('utils')