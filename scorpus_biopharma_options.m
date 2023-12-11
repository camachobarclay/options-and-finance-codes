clear;
clc;
close all;

diary BSCex

%% Defining the parameters of the warrant setup

A = 10^9;           % Total company assets, an assumed quantity subject to change
n = 10^6;           % Number of total shares outstanding, an assumed quantity subject to change
P = 10^5;           % Principle amount invested
un_pr = 3;          % Price of unit = 1 share common stock, 2 Series A warrants
K = 6;              % Series A Warrant strike price
m = 100000/un_pr;   % Number of Series A units bought
SA_m = m*2;         % Number of Series A warrants available to exercise
                    % Determining the corrective factor cf as required by BSC
cf = n/(SA_m+n);    % Note that corrective factor can be modified to something other than the share dilution n/(m+n)

format bank;

fprintf('Assumptions:\n\n')
fprintf('Number of shares n outstanding post sale of %f units: n = %d.\n\n',m,n);
fprintf('(N.B. We assume this value. This is important to know in order to get a sense\nof how much to dilute the n shares after m warrants are exercised.)\n\n');
fprintf('Series A warrant strike price K: K = $%.2f.\n',K);
format;
fprintf('Number of series A warrants SA_m assumed exercised: SA_m = %f.\n',SA_m);

%% Parameters for the blsprice function

sigma = 0:0.01:2;                       % Volatility rate for the company uniformly distributed
sigma_SP = [23.37450975, 41.28723142];
sigma_SP = sigma_SP/100;                % Volatility rate for company based on SP500 companies
R = 0.29; R = R/100;                    % The risk free rate of a 5 yr Treasury note per annum
T = 5;                                  % Time to expiry
prices = [5,6];                         % Simulated stock prices of interest
delta = 0;                              % Dividend (which is always 0 for this company)

% Use these for loops and indexing

M = length(prices);
N = length(sigma);
L = length(sigma_SP);

fprintf('SP500 average, median implied volatility rates sigma_SP: sigma_SP = %.2f%%, %.2f%%.\n', sigma_SP(1)*100, sigma_SP(2)*100);
fprintf('Simulated implied volatility rates sigma: sigma = %d%%, %d%%,..., %d%%, %d%%.\n',sigma(1)*100,sigma(2)*100,sigma(N-1)*100,sigma(N)*100);
fprintf('Risk free premium interest rate R of a 5 year treasury note: R = %.2f%%.\n', R*100);
fprintf('Time T until expiration of the warrants: T = %d years.\n',T);
fprintf('Test prices for the underlying shares S: S = ');

for k=1:M
    fprintf('$%.2f',prices(k));
    if k<M
        fprintf(', ');
    end
end

fprintf('.\n\n');

format bank

% Vectors used for graphing

War_pr_vals = [zeros(size(sigma)), zeros(size(sigma_SP))];
Opt_pr_vals = War_pr_vals;

War_pr_vals_cm = zeros(size(sigma));
Opt_pr_vals_cm = War_pr_vals_cm;

comb_sigma = [sigma sigma_SP];
[comb_sigma, ind] =sort(comb_sigma);

%% Running Black Scholes

for j = 1:M
    
    pr = prices(j);
    

    fprintf('\n-------------------------\n\nComputed using MATLAB\n');
    fprintf('\nStock Price\t\t\tVolatility\t\tWarrant Price\t\tUndiluted Option Price\n\n');
    
    % Calculating Black Scholes for the simulated implied volatilitis listed in sigma, sigma_SP
    
    for k = 1:N
        sd = sigma(k);
        [War_pr,~] = blsprice(pr, K, R, T, sd);
        Opt_pr = War_pr;                % Saving the base option price
        War_pr = War_pr*cf;             % Applying the corrective for dilution to get the warrant price
        Opt_pr_vals(k) = Opt_pr;        % Saving values
        War_pr_vals(k) = War_pr;        
        fprintf('$%.2f\t\t\t\t%f\t\t$%f\t\t\t$%f\n',pr,sd,War_pr,Opt_pr);
    end
        
    fprintf('\n');
    
    % Calculating Black Scholes for the SP500 based implied  volatilities
    
    for l = 1:L
        sd = sigma_SP(l);
        [War_pr,~] = blsprice(pr, K, R, T, sd);
        Opt_pr = War_pr;                % Saving the base option price
        War_pr = War_pr*cf;             % Applying the corrective for dilution to get the warrant price
        Opt_pr_vals(l+N) = Opt_pr;      % Saving values
        War_pr_vals(l+N) = War_pr;
        fprintf('$%.2f\t\t\t\t%f\t\t$%f\t\t\t$%f\n',pr,sd,War_pr,Opt_pr);
    end

    fprintf('\n-------------------------\n\nComputed Manually\n');
    fprintf('\nStock Price\t\t\tVolatility\t\tWarrant Price\t\tUndiluted Option Price\n\n');
    
    %% Repeat the simulation using manual calculation

    for k = 1:N
        sd = sigma(k);
        d1 = log(pr/K) +(R - delta + 0.5*sd^2)*T;               % Calculating d1 according to the BSC formula
        d1 = d1/(sd*sqrt(T));                                   % Calculating d1 according to the BSC formula
        d2 = d1 - sd*sqrt(T);                                   % Calculating d2 according to the BSC formala
        War_pr = pr*normcdf(d1) - K*exp(-R*T)*normcdf(d2);      % Calculating the option price manually
        Opt_pr = War_pr;                                        % Saving the base option price
        War_pr = War_pr*cf;                                     % Applying the corrective for dilution to get the warrant price
        Opt_pr_vals_cm(k) = Opt_pr;                             % Saving values
        War_pr_vals_cm(k) = War_pr;
        fprintf('$%.2f\t\t\t\t%f\t\t$%f\t\t\t$%f\n',pr,sd,War_pr,Opt_pr);
    end
    
    fprintf('\n');
    
    for l = 1:L
        sd = sigma_SP(l);
        d1 = log(pr/K) +(R - delta + 0.5*sd^2)*T;               % Calculating d1 according to the BSC formula.
        d1 = d1/(sd*sqrt(T));                                   % Calculating d1 according to the BSC formula.
        d2 = d1 - sd*sqrt(T);                                   % Calculating d2 according to the BSC formala
        War_pr = pr*normcdf(d1) - K*exp(-R*T)*normcdf(d2);      % Calculating the option price manually.
        Opt_pr = War_pr;                                        % Saving the base option price
        War_pr = War_pr*cf;                                     % Applying the corrective for dilution to get the warrant price
        Opt_pr_vals_cm(l+N) = Opt_pr;                           % Saving values
        War_pr_vals_cm(l+N) = War_pr;
        fprintf('$%.2f\t\t\t\t%f\t\t$%f\t\t\t$%f\n',pr,sd,War_pr,Opt_pr);
    end
    
 %% Plot Results   
    figure;
    hold on
    plot(comb_sigma,War_pr_vals(ind),'LineWidth',5)
    plot(comb_sigma, War_pr_vals_cm(ind), '--g', 'LineWidth',2.5)
    plot(comb_sigma, Opt_pr_vals(ind), 'magenta','LineWidth',5)
    plot(comb_sigma, Opt_pr_vals_cm(ind),'--c','LineWidth',2.5);
    
    xlabel('Implied Volatility Values');
    ylabel('Price');
    legend('Warrant Price','Warrant Price Manually Computed','Option Price','Option Price Manually Computed') 
    title(['Derivative vs. Volatility Graph when Stock Price = ', num2str(pr)]);

    grid on 
    grid minor
    set(gca,'xtick',linspace(0, sigma(end), 5))
    set(gca,'ytick',linspace(0,pr,pr+1)) 
    
    fprintf('\n');
    
end

%% Produce 3D Plots of the the option prices vs. two variables

% When sigma and the stock price are allowed to vary

figure;

Prices = 5:0.5:20;
sigma = 0:.1:2;
[SD, PR] = meshgrid(sigma, Prices);

BSCpr = @(x,y) blsprice(y, K, R, T, x);

OPT_PR = BSCpr(SD, PR);
surf(SD,PR,OPT_PR);

xlabel('Implied Volatility Rates');
ylabel('Stock Price');
zlabel('Option Price');
title('Option Prices vs. Variable Volatility, Stock Price')

% When sigma and the time until expiration are allowed to vary (when the
% stock price is 5, 6 dollars.

Tvals = 0:0.25:5;
[SD, TVALS] = meshgrid(sigma, Tvals);

for l = 1:L
    
    figure;

    pr = prices(l);
    BSCt = @(x,y) blsprice(pr, K, R, y, x);
    OPT_PR = BSCt(SD, TVALS);
    surf(SD,TVALS,OPT_PR);
    
    xlabel('Implied Volatility Rates');
    ylabel('Time in Years');
    zlabel('Option Price');
    title(['Option Prices vs. Variable Volatility, Time when the Stock Price = ', num2str(pr)]);
    
end

diary off
