%% Tidying up the workspace

clear;
close all;
clc;

%% Setting up initial parameters using a call Bloomber API Call
warning off;
SPAC = readcell('C:\blp\data\W_11_SPAC_Monitor_Copy');
warning on;

[m, n] = size(SPAC);

%%

ColHeads = cell(n,1);

for k = 1:n
    ColHeads{k} = SPAC{1,k};
end

%%
tickerind = zeros(m,1);
ticker = '';
war_ticker = ticker;
for j = 1:m
    ticker = SPAC{j,1};
    if all(ticker(1,end-1:end)=='US')
        SPAC{j,1} = strcat(ticker,' EQUITY');
        tickerind(j) = j;
    end
end

[tickerind,~] = find(tickerind);         
SPAC = SPAC(tickerind,:);
[m, n] = size(SPAC);

tickerind = zeros(m,1);
ticker = '';
kval = 0;
%war_exist = 0;
price_exist = 0;
cashstring = '';
cashoff = 0;
existenceRecord = zeros(m,2);
for j = 1:m
%     ticker = SPAC{j,13};
%     if ischar(ticker)
%         SPAC{j,13} = strcat(ticker,' EQUITY');
%         j
%         war_exist = 1
%     else
%         continue;
%     end
    
    kval = SPAC{j,14};
    if isnumeric(kval)
        price_exist = 1;
    end
    
    
    cashoffstr = SPAC{j,5};
    cashoff = str2num(cashoffstr(1:end-3));
    if cashoffstr(end-2) =='B'
        cashoff = cashoff*10^3;
    end
    SPAC{j,5} = cashoff;
    tickerind(j) = j;
    
end

[tickerind,~] = find(tickerind);         
SPAC = SPAC(tickerind,:);
[m, n] = size(SPAC);

testint = 1;
indlist = randi(m,testint,1);
loopind = 0;


start_date = '09/23/2020';
end_date = '09/23/2021';

for i = 1:testint
    loopind = indlist(i)
    security = SPAC{loopind,1}
    fprintf('Underlying Security: %s\n', security)

    A = SPAC{loopind,5}
    K = SPAC{loopind,14}



stock_data_req = {'LAST_PRICE';  'CURRENT_SHARES_OUTSTANDING_RT';...
    'WRT_EXPIRE_DT';'SECOND_UNDERLYING_TICKER'};

    %'FUND_TOTAL_ASSETS'; 'FUND_TOTAL_ASSETS';...
    %'FUND_CRNCY_ADJ_TOTAL_ASSETS';'FUND_NET_ASSET_VAL';...
% 'BS_TOT_ASSET'; 'NET_ASSETS';...
%     'BS_CASH_NEAR_CASH_ITEM'; 'CURR_ENTP_VAL'; 'ENTERPRISE_VALUE';...
%     'SHORT_AND_LONG_TERM_DEBT'; 'CUR_MKT_CAP'; 'EQY_DVD_YLD_12M';...
%     'PREFERRED_EQUITY_&_MINORITY_INT';
war_data_req = {'OPT_STRIKE_PX'; 'STRIKE_PRICE_REALTIME'; 'WRT_DAYS_EXPIRE';...
     'OPT_DAYS_EXPIRE';'LAST_PRICE';'WRT_IMPLIED_VOLATILITY_LAST';...
     'DERIVATIVE_STRIKE_PRICE'};

javaaddpath('C:\blp\BloombergWindowsSDK\JavaAPI\v3.15.1.1\lib\blpapi3.jar')

% starting blp
c = blp;

% Getting bulk data on the chosen security, warrant
[stock_data,~] = getdata(c,security,stock_data_req);
warrant = stock_data.SECOND_UNDERLYING_TICKER{1};

[warrant_data,~] = getdata(c,warrant,war_data_req);
% Retrieving historical end of day price data for SLCR based on above variables
[price, sec] = history(c,security, 'LAST_PRICE',...
    start_date, end_date);

% closing blp
close(c);

%% Printing values from the stock/warrant structs

% Printing Stock and Warrant Struct Values

format bank
fprintf('Stock Struct Parameters:\n')
stock_data_req = fieldnames(stock_data);
L = length(stock_data_req);
fprintf('\n');
for j = 1:L
    c = eval(['stock_data.',stock_data_req{j}]);
    if isnumeric(c)
        fprintf([stock_data_req{j},': %.2f\n'],c);
    end
end

fprintf('\n');

end
% Importing Stock Struct Values


fprintf('Warrant Struct Parameters:\n')
fprintf('\n');
war_data_req = fieldnames(warrant_data);
L = length(war_data_req);
for j = 1:L
    c = eval(['warrant_data.',war_data_req{j}]);
    if isnumeric(c)
        fprintf([war_data_req{j},': %.2f\n'],c);
    end
end

fprintf('\n');

%% Assigning struct values from the stock/warrant structs to variables

S = stock_data.LAST_PRICE;                           % Stock price assignment
N = stock_data.CURRENT_SHARES_OUTSTANDING_RT;        % Number of SLCR shares outstanding
%EQ_DIV = stock_data.EQY_DVD_YLD_12M;                 % Purported dividend yield rate

TotAssets = stock_data.BS_TOT_ASSET;                 % Potential value for A = Total Assets
NetAssets = stock_data.NET_ASSETS;                   % ""


% Download basic financial data to compute company cash value using the
% formula EV = MCAP - Cash + Pref + Debt, or rather Cash = MCAP - EV + Pref
% + Debt.

MCAP = stock_data.CUR_MKT_CAP/10^6;                   
CEV = stock_data.CURR_ENTP_VAL;
EV = stock_data.ENTERPRISE_VALUE;
PE_MINT = stock_data.PREFERRED_EQUITY_a_MINORITY_INT;
CCE = stock_data.BS_CASH_NEAR_CASH_ITEM;
DEBT = stock_data.SHORT_AND_LONG_TERM_DEBT;

% Importing Warrant Struct Vals

OPT_K  = warrant_data.OPT_STRIKE_PX; 
OPT_KRT = warrant_data.STRIKE_PRICE_REALTIME;
DER_K = warrant_data.DERIVATIVE_STRIKE_PRICE;
WDE = warrant_data.WRT_DAYS_EXPIRE;
ODE = warrant_data.OPT_DAYS_EXPIRE;
WRT_lpx = warrant_data.LAST_PRICE;
WRT_imvol = warrant_data.WRT_IMPLIED_VOLATILITY_LAST;


%% Setting up time until expiration

T = 0;

if isnumeric(WDE)&&WDE>=0
    T = WDE;       % Warrant strike price
    WDEsucc = 1;
elseif isnumeric(ODE)&&ODE>=0
    warning('Non-valid value for WRT_DAYS_EXPIRE. Using OPT_DAYS_EXPIRE instead.')
    T = ODE;
    ODEsucc = 1;
end

T = T/365.25;     % Years until warrant expiration

if ~(WDEsucc||ODEsucc)
    T = 3;
    warning('Warrant expiration data missing, using default value (3 years) for maturity.');
end
% 
% %% Setting up the total asset value A for the company
% 
% AssetVal = 0;
% CompCashAssets = 0;
% AssetValsucc = 0;
% CompCashsucc = 0;
% EVsucc = 0;
% CEVsucc = 0;
% PMIDsucc = 0;
% NEVValsucc = 0;
% 
% if isnumeric(TotAssets)&&(TotAssets>=0)
%     AssetVal = TotAssets;                       % Prefer to use total asset val
%     AssetValsucc = 1;
%     % Computing company cash value using the formula CASH = MCAP - EV +
%     % PREF + DEBT.
% elseif isnumeric(MCAP)&&(MCAP>=0)
%     if isnumeric(CEV)&&(CEV>=0)
%         CompCashAssets = MCAP - CEV;
%         NEV = CEV;
%         CEVsucc = 1;
%     elseif isnumeric(EV)&&(EV>=0)
%         CompCashAssets = MCAP - EV;
%         NEV = EV;
%         EVsucc = 1;
%     end
%     if CEVsucc||EVsucc
%         CompCashsucc = 1;
%         warning('Using computed value for total company value, i.e., computed cash.');
%     end
%     if (isnumeric(DEBT)&&(DEBT>=0))&&(isnumeric(PE_MINT)&&(PE_MINT>=0))
%         CompCashAssets = CompCashAssets + DEBT + PE_MINT;
%         NEV = NEV - DEBT;
%         PMIDsucc = 1;
%     end
%     
% end
% 
% clist = ["AssetVal";
%          "CompCashAssets";
%          "BBGCashVal"];
% 
% cnum = length(clist);
% 
% CashVals = zeros(2,cnum);
% CashVals(2,:) = 1:cnum;
% 
% if AssetValsucc
%     CashVals(1,1) = AssetVal;
% end
% 
% if CompCashsucc
%     CashVals(1,2) = CompCashAssets;
% end
% 
% if isnumeric(CCE)&&(CCE>=0)
%     CashVals(1,3) = CCE;
% end
% 
% newind = find(CashVals(1,:));
% cnum = length(newind);
% 
% CashVals = CashVals(:,newind);
% clist = clist(newind);
% Asucc = 0;
% if any(CashVals(1,:))
%     A = CashVals(1,1);
%     powerdiff = zeros(1,nchoosek(cnum,2));
%     counter = 0;
%     Asucc = CashVals(2,1);
%     for j = 1:cnum
%         for k = j+1:cnum
%             counter = counter + 1;
%             powerdiff(counter) = abs(floor(log10(CashVals(j))) - floor(log10(CashVals(k))));
%             fprintf('Difference in signficant digits between ');
%             fprintf(clist(j));
%             fprintf(' and ');
%             fprintf(clist(k));
%             fprintf('.\nDifference in power: %d.\n',powerdiff(counter));
%         end
%     end
% 
%     if counter
%         fprintf('\n');
%     end
% else
%     A = S*N;
%     warning('Failed to find firm asset value, using computed stock price S instead.')
% end
% 
% %%
% K = OPT_KRT;            % Strike Price
% 
% r = [1.50];       % Annualized risk free interest rate
% r = r/100;              % Bl-Scholes function requires decimal form for r 
% 
% %m = [10^3*(1:9), 10^4*(1:9), 100000:10000:290000, 10^5*[3:9], 10^6]'; 
% m = 10^5;
% 
% %% Computing annualized volatility for SLCR based on 
% mystdcalc = @(u,ubar) sqrt((sum((u-ubar).^2)/(length(u) - 1)));
% sq252 = sqrt(252);
% 
% % First compute daily returns percentages.
% daily_returns = 100*((price(2:end) - price(1:end-1))./price(1:end-1));
% log_daily_returns = log(price(2:end)./price(1:end-1));
% 
% dr_mean = mean (daily_returns);
% ldr_mean = mean(log_daily_returns);
% 
% dr_sigma = mystdcalc(daily_returns,dr_mean)*sq252;
% ldr_sigma = mystdcalc(log_daily_returns,ldr_mean)*sq252*100;
% 
% %sigma = std(daily_returns)*sq252;  % Scale by sqrt(252) to annualize sigma
% sigma = ldr_sigma;
% sigma = sigma/100;                     % Bl-Scholes function requires decimal form for sigma 
% 
% 
% 
% 
% %% Running the Black-Scholes Algorithm for multiple values of m = # of warrants issued
% 
% % Saving output in a file and printing column/table headings
% delete 'SLCR Warrant Pricing';
% diary 'SLCR Warrant Pricing';
% fprintf('SLCR US Equity Warrant Pricing %s\n\n',date);
% fprintf('m = # of Warrants Exercised\tComputed (Call) Warrant Price\n')
% 
% for k = 1:length(m)
%     M = m(k);
%     [call,put] = blsprice(A/N,K,r,T,sigma);  % Running Black-Scholes
%    
%     % Normalizing the prices to reflect dilution after warrants are
%     % exercised.
%     q = M/N;
%     call = call/(1+q);                     
%     put = put/(1+q);
%     
%     % printing more results
%     fprintf('%i\t\t\t\t',m);
%     fprintf('$%.4f\n',call);
% end
% 
% % 
% diary off;
% 
% 
