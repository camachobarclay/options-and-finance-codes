clear;
close all;
clc;

%% Setting up initial parameters (taken from Bloomberg) 

S = 9.67;                           % SLCR stock price
AssetVal = 330.7*10^6;              % Sum total of SLCR's assets
n = 30.5*10^6;                      % Number of SLCR shares outstanding
div = 0;                            % SLCR dividend yield rate

% Picking a time period to work with 

security = 'KRNL US Equity';
start_date = '09/16/2020';
%start_date = '01/02/2020';
end_date = '09/16/2021';

%% Retrieving historical end of day price data for SLCR based on above variables

javaaddpath('C:\blp\BloombergWindowsSDK\JavaAPI\v3.15.1.1\lib\blpapi3.jar')
c = blp;
[price, sec] = history(c,security, 'LAST_PRICE',...
    start_date, end_date);
close(c);

price = price(:,2);     % Getting rid of extraneous outputs

%% Computing annualized volatility for SLCR based on 

% First compute daily returns percentages.
daily_returns = 100*((price(2:end) - price(1:end-1))./price(1:end-1));
sigma = std(daily_returns)*sqrt(252);  % Scale by sqrt(252) to annualize sigma
sigma = sigma/100;                     % Bl-Scholes function requires decimal form for sigma 

%% Setting up other parameters

K = 11.5;       % Warrant strike price

T = 1962/365.25;     % Years until warrant expiration

%r = 3
r = 1.50;       % Annualized risk free interest rate
r = r/100;      % Bl-Scholes function requires decimal form for r 

M = [10^3*(1:9), 10^4*(1:9), 100000:10000:290000, 10^5*[3:9], 10^6]'; 

%% Running the Black-Scholes Algorithm for multiple values of m = # of warrants issued

% Saving output in a file and printing column/table headings
delete 'KRNL Warrant Pricing';
diary 'KRNL Warrant Pricing';
fprintf('KRNL US Equity Warrant Pricing %s\n\n',date);
fprintf('m = # of Warrants Exercised\tComputed (Call) Warrant Price\n');

for k = 1:length(M)
    m = M(k);
    [call,put] = blsprice(AssetVal/n,K,r,T,sigma);  % Running Black-Scholes
   
    % Normalizing the prices to reflect dilution after warrants are
    % exercised.
    q = m/n;
    call = call/(1+q);                     
    put = put/(1+q);
    
    % printing more results
    fprintf('%i\t\t\t\t',m);
    fprintf('$%.4f\n',call);
end

diary off;


%%