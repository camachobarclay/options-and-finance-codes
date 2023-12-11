%% Tidying up the workspace, adding paths

clear;
close all;
clc;

javaaddpath('C:\blp\BloombergWindowsSDK\JavaAPI\v3.15.1.1\lib\blpapi3.jar');

%% Uploading SPAX Monitor from local drive (must be locally saved to the below address).

warning off;
SPAC = readcell('C:\blp\data\W_SPAC_Monitor.xlsx');
warning on;



%% Massaging SPAC data for later use

% Filtering out first column of the SPAX monitor so that it contains only
% equities

[m, n] = size(SPAC);
ColHeads = cell(n,1);
Mspaceval = floor(log10(n));

for k = 1:n
    ColHeads{k} = SPAC{1,k};
    for kk = 1:(Mspaceval-floor(log10(k)))
        fprintf(' ');
    end
    fprintf('%i. %s\n',k,ColHeads{k});
end
tickerind = zeros(m,1);
ticker = '';

for j = 1:m
    ticker = SPAC{j,1};
    if all(ticker(1,end-1:end)=='US')
        SPAC{j,1} = strcat(ticker,' EQUITY');
        tickerind(j) = j;
    end
end

tickerind = nonzeros(tickerind);
SPAC = SPAC(tickerind,:);       
[m, ~] = size(SPAC);

% Finding indices Eliminating SPACs that don't have a (nonzero) strike
% price and cash offer.

tickerind = zeros(m,1);
ticker = ''; cashstring = '';
kval = 0; colwrtpx = 1; price_exist = 0; cashoff = 0; 
billcount = 0; billcheck = 0;
cashcond1 = 0; cashcond2 = 0; cashcond3 = 0;
filterm = 150; filterwrt = 1;

for j = 1:m
    % Making sure the strike price exists and is valid.
    kval = SPAC{j,14};
    if ~isnumeric(kval)||~(kval>0)
        continue;
    end
    
    % Filtering out warrants that are not cheap enough.
    colwrtpx = SPAC{j,16};
    if ~isnumeric(colwrtpx)
        if (colwrtpx > filterwrt)
            fprintf('Warrant price for SPAC %s is greater than $%.2f limit.',...
                SPAC{j,1}, filterwrt);
            continue;
        end
    end
    
    % Check to see if cash offer exists, skip SPAC if non-existent.
    cashoffstr = SPAC{j,5};
    cashoff = str2double(cashoffstr(1:end-3));
    
    % Making sure we operate on the same scale (convert to million if
    % billion).
    if (cashoffstr(end-2) =='B')||(cashoffstr(end-2)=='b')
        cashoff = cashoff*10^3;
        billcheck = 1;
    end
    
    % Skip cash size offers that are not large enough
    cashcond1 = isnumeric(cashoff);
    cashcond2 = (~isempty(cashoff));
    cashcond3 = (cashoff>=filterm);

    if cashcond1&&cashcond2&&cashcond3
        billcount = billcount + billcheck;
        billcheck = 0;
        SPAC{j,5} = cashoff;
    else
        if ~cashcond1
            fprintf('Non-numeric value for SPAC %s cash offer',SPAC{j,1});
        elseif ~cashcond2
            fprintf('Empty value for SPAC%s cash offer',SPAC{j,1});
        else
            fprintf('Cash offer value for SPAC %s is less than %fM threshold',...
                SPAC{j,1}, filterm);
        end
        fprintf('\n');
        continue;
    end

    tickerind(j) = j;
end

fprintf('\nThe total number of IPO cash offers that exceeded a billion dollars is %i.\n',...
    billcount);

tickerind = nonzeros(tickerind);         
SPAC = SPAC(tickerind,:);
[m, ~] = size(SPAC);

%% Setting up parameters to run the driver loop 

% Setting up parameters to run through the loop
testint = 5;
index = randi(m,testint,1);
%index = 1:m;
SPACind = 0;
p = 0;

%% Running the driver loop
for p = 1:testint
    
    SPACind = index(p);
    
    % Assigning the BSC parameters
    security = SPAC{SPACind,1};     % Security Ticker
    A = SPAC{SPACind,5};            % Total cash offer in millions
    K = SPAC{SPACind,14};           % Warrant Strike Price
    
    % Assigning FLDS for stock, warrant import from BBG

    stock_data_req = {'LAST_PRICE';  'CURRENT_SHARES_OUTSTANDING_RT';...
        'WRT_EXPIRE_DT';'SECOND_UNDERLYING_TICKER'; 'EQY_INIT_PO_DT';...
     'MERGERS_AND_ACQUISITIONS'};
    war_data_req = {'OPT_STRIKE_PX'; 'STRIKE_PRICE_REALTIME'; 'WRT_DAYS_EXPIRE';...
     'OPT_DAYS_EXPIRE';'LAST_PRICE';'WRT_IMPLIED_VOLATILITY_LAST';...
     'DERIVATIVE_STRIKE_PRICE'; 'WRT_ISSUE_AMT'; 'WRT_OUTSTANDING'};
    
    %% Importing values from BBG through API call

    c = blp;

    % Getting bulk data on the chosen security, warrant
    [stock_data,~] = getdata(c,security,stock_data_req);
    
    warrant = strcat(stock_data.SECOND_UNDERLYING_TICKER{1},' Equity');    
    [warrant_data,~] = getdata(c,warrant,war_data_req);
    
    % Skipping warrants that are more that $1.00.
    WRT_lpxsucc = 0;
    WRT_lpx = warrant_data.LAST_PRICE;
    if isnumeric(WRT_lpx)
        if (WRT_lpx>filterwrt)
            warning('Warrant price for SPAC %s is greater than $%.2f limit.',...
                security, filterwrt);
            continue;
        elseif (WRT_lpx>0)
            WRT_lpxsucc = 1;
        end
    end

    MA = stock_data.MERGERS_AND_ACQUISITIONS;
    MA = MA{:,:};
    MAm = 0;
    MAn = 0;
    [MAm, MAn] = size(MA);
    MA_DATES = zeros(MAm,1);
    
    if MAm
        for k = 1:MAm
            MA_DATES(k) = MA{1,3};
        end
        MA_DATES = sort(MA_DATES);
        dates_cell = cell(MAm,1);
        for k = 1:MAm
            dates_cell{k} = datestr(datetime(MA_DATES(k),'ConvertFrom','datenum'),'mm/dd/yyyy');
        end
        % dates_cell{:}
    end
    %MA_DATES = datetime(MA_DATES);
    

    %% Setting up start and end dates for historical price collection
    
    % Counters to track success.
    
    % SD_Lsucc = 0;
    SD_DTsucc = 0; 
    SD_succ = 0; 
    API_SDDTsucc = 0;
    API_SDsucc = 0;

    % Loading the start date
    start_date = SPAC{SPACind,7};
    OSD = start_date;

    % Trying to load the start date using datetime, making sure that
    % doesn't give trouble.

    try
        start_date = datetime(start_date,'InputFormat','MM/dd/yyyy');
        SD_DTsucc = 1;
        start_date = datestr(start_date,'mm/dd/yyyy');
        SD_succ = 1;
    catch
        if ~SD_DTsucc
            warning('Could not perform datetime operations on %s EQY_INIT_PO_DT value.',...
                security);
        elseif ~SD_succ
            warning('Could not perform datestr operations on %s start_date value.',...
                security);
        end
    end    

    % Trying to load the start date using datetime, making sure that
    % doesn't give trouble.
    try
        API_start_date = datetime(stock_data.EQY_INIT_PO_DT,...
            'ConvertFrom','datenum');
        API_SDDTsucc = 1;
        API_start_date = datestr(API_start_date,'mm/dd/yyyy');
        API_SDsucc = 1;
    catch
        if ~API_SDDTsucc
            warning('Could not perform datestr operations on %s EQY_INIT_PO_DT value.',...
                    security);
        elseif ~API_SDsucc
            warning('Could not perform datetime operations on %s EQY_INIT_PO_DT value.',...
                security);
        end
    end
    
    % Running through the different start_date, API_date combinations.
    if SD_succ||API_SDsucc      
    elseif (~SD_succ)&&API_SDsucc
        start_date = API_start_date;
        fprintf('Program is assigning API_start_date to pricing start_date variable for SPAC %s',...
            security);
    elseif (~SD_succ)&&(~API_SDsucc)
        warning('No dates were loaded for SPAC %s.',security);
        continue;
    else
        warning('An error occured in initial pricing date assignment section for SPAC %s.',...
            security);
        continue;
    end
    
    % Finalizing states for price history
    start_date = datestr(start_date,'mm/dd/yyyy');
    end_date = datestr(datetime(date)-1,'mm/dd/yyyy');

    % Not including SPACs that don't have enough trading history.
    if daysact(start_date,end_date)<2
        warning('Not enough stock price history for %s to proceed with trial.',...
            security); 
        continue;
    end
    
  %%  
    % Retrieving historical end of day price data for SLCR based on above variables
    [Price, sec] = history(c,security, 'LAST_PRICE',...
        start_date, end_date);
    
    % closing blp
    close(c);

    price = Price(:,2);
    dates = datetime(Price(:,1),'ConvertFrom','datenum');

    %% Assigning struct values from the stock/warrant structs to variables
    S_succ = 1;
    S = stock_data.LAST_PRICE;      % Stock price assignment
    if ~(isnumeric(S))
        S = price(end);
        warning('Struct parameter S for stock_data is non numerical. Assigning last end of day trading price value.')
        S_succ = 0;
    end

    if S_succ
        price = [price; S];
        end_date = datestr(date,'mm/dd/yyyy');
        dates = [dates;date];
    end

    if MAm
        dind = zeros(MAm,1);
        for k = 1:MAm
            dind(k) = find(dates == datetime(dates_cell{k}));
        end
    end
    
    % Number of shares outstanding
    N = stock_data.CURRENT_SHARES_OUTSTANDING_RT;       

    % Importing (remaining) Warrant Struct Vals
    
    OPT_K  = warrant_data.OPT_STRIKE_PX; 
    OPT_KRT = warrant_data.STRIKE_PRICE_REALTIME;
    DER_K = warrant_data.DERIVATIVE_STRIKE_PRICE;
    WDE = warrant_data.WRT_DAYS_EXPIRE;
    ODE = warrant_data.OPT_DAYS_EXPIRE;
    WRT_imvol = warrant_data.WRT_IMPLIED_VOLATILITY_LAST;
    WRT_count = warrant_data.WRT_OUTSTANDING;
    WRT_amount = warrant_data.WRT_ISSUE_AMT;
    %% Setting up other BSC parameters

    T = 0;
    
    if isnumeric(WDE)&&WDE>=0
        T = WDE;       % Warrant strike price
        WDEsucc = 1;
    elseif isnumeric(ODE)&&ODE>=0
        warning('Non-valid value for WRT_DAYS_EXPIRE. Using OPT_DAYS_EXPIRE instead.')
        T = ODE;
        ODEsucc = 1;
    else
        warning('Time to expiration does not exist. Skipping SPAC %s.',security);
    end
    
    T = T/365.25;     % Years until warrant expiration
    
    r = [1.50];             % Annualized risk free interest rate
    r = r/100;              % Bl-Scholes function requires decimal form for r 
    tol = 10^-7;            % Tolerance for the implied volatility calcilation
    MV = [100000:10000:300000]'; 
    %MV = 10^5;
    
    %% Computing annualized volatility for SLCR based on 

    [DR, MEAN, SIGMA] = sigmacalc(price);
    [daily_returns, dr_mean, dr_sigma] = drsigmacalc(price);
    [log_daily_returns, ldr_mean, ldr_sigma] = ldrsigmacalc(price);
    sigma = SIGMA;
    
  %% Printing values from the stock/warrant structs, printing results, graphs.
    
    fprintf('Security Name: %s\n',SPAC{SPACind,2});
    figure
    hold on;
    plot(dates,price);
    xlabel([security,'''s Entire Trading History']);
    ylabel('Share Price');
    title([security,' Closing Price/Entire Trading History']);
    hold off;
    fprintf('Security Ticker: %s\n', security);
    fprintf('Warrant Ticker: %s\n', warrant);
    if WRT_lpxsucc
        fprintf('Warrant Price (in dollars): %.2f\n', WRT_lpx);
    end
    fprintf('Stock Price: $%.2f\n',S);
    fprintf('Cash Offer (in millions): $%i\n',A);
    fprintf('Number of shares outstanding (in millions): %i\n', N);
    fprintf('SPAC Asset Value/share: $%.4f\n',A/N);
    fprintf('Warrant Strike Price: $%.2f\n',K);
    if isnumeric(WRT_count)
        fprintf('Number of warrants outstanding: %i\n',WRT_count);
    end
    if isnumeric(WRT_amount)
        fprintf('Number of warrants issued (in thousands): %i\n',WRT_amount);
    end  
    fprintf('Time until maturity expiration: %i days | %.4f years\n',...
        T*365.25,T);
    fprintf('Annualized Volatility: %%%.4f\n', sigma*100);
    figure
    hold on;
    bar(dates(2:end), DR*100);
    xlabel([security,'''s Entire Trading History']);
    ylabel('100*(S_{close} - S_{prior})/S_{prior}');
    title([security,' Daily Returns (% form)']);
    hold off;
    if WRT_lpxsucc&&S_succ
        mybls = @(z) blsprice(S,K,r,T,z);
        imvol = IMVOL(mybls,WRT_lpx,tol);
        fprintf('MATLAB calculated implied volatility: %%%.4f\n',imvol*100);
    else
        fprintf('BBG current stock and warrant price do no exist simultaneously. Cannot calculate current implied volatility.\n');
    end
    
    fprintf('Number of trading day prices used to compute volatility: %i\n',...
        length(price));
    fprintf('Annualized risk free rate: %%%.4f\n',r*100);
    fprintf('SPAC IPO Date (SPAC Monitor): %s\n',start_date);
    fprintf('Last close date: %s\n',end_date);
    if ischar(OSD)
        fprintf('Original SPAX date column value: %s\n',OSD);
    end
    if ischar(API_start_date)
        fprintf('SPAC IPO Date (API Call): %s\n',API_start_date);
    end
    
    % Printing Stock Struct Values
    format bank
    fprintf('\nStock Struct Parameters:\n');
    stock_data_req = fieldnames(stock_data);
    L = length(stock_data_req);
    
    for j = 1:L
        c = eval(['stock_data.',stock_data_req{j}]);
        if isnumeric(c)
            fprintf([stock_data_req{j},': %.4f\n'],c);
        end
    end
    
    % Printing Warrant Struct Values
    fprintf('\nWarrant Struct Parameters:\n')
    war_data_req = fieldnames(warrant_data);
    L = length(war_data_req);
    for j = 1:L
        c = eval(['warrant_data.',war_data_req{j}]);
        if isnumeric(c)
            fprintf([war_data_req{j},': %.4f\n'],c);
        end
    end

    fprintf('\n');
    
    %% Computing the bare call price for the warrant using Black Scholes.

    [call,~] = blsprice(S,K,r,T,sigma);
    
    fprintf('For inputs:\n');
    fprintf('Stock price S = $%.2f\n',S);
    fprintf('Strike price K = $%.2f\n',K);
    fprintf('Annualized risk free rate R = %%%.4f\n',r*100);
    fprintf('Time T = %.4f years\n',T);
    fprintf('Annualized volality sigma = %%%.4f\n\n',sigma*100);
    fprintf('Black-Scholes MATLAB calculated call warrant price is $%.2f.\n',call);
    if WRT_lpxsucc
        percalculate = @(x,y) (y./x-1)*100;
        perupdownside = percalculate(WRT_lpx, call);
        fprintf(['The upside/downside in the computed warrant prices\nvs. ',...
            'current market prices is %%%.4f.\n'],perupdownside);
    end
    WRT_intrval = max(0,S-K); TimeVal = call - WRT_intrval;
    fprintf('The intrinsic/time value of the warrant is $%.4f/$%.4f.\n',...
        WRT_intrval,TimeVal);
    
    %% Computing the associated (European) Greek variables.

    delta = 0; theta = 0; gamma = 0;  vega = 0; rho = 0; 
    
    d1 = D1(S,K,r,sigma,T);
    d2 = D2(d1,sigma,T);
    
    delta = DELTA(d1);
    theta = THETA(S,K,r,sigma,T,d1,d2);
    gamma = GAMMA(S,sigma,T,d1);
    vega = VEGA(S,T,d1);
    rho = RHO(K,r,T,d2);

    fprintf('Greek value MATLAB calculations for when options are European style:\n');
    fprintf('DELTA = %.4f\n',delta);
    fprintf('THETA = %.4f\n',theta);
    fprintf('GAMMA = %.4f\n',gamma);
    fprintf('VEGA = %.4f\n',vega);
    fprintf('RHO = %.4f\n',rho);

    fprintf('\n');
    
    %% Computing the diluted warrant prices for different number of warrants m.
    
    fprintf('M = # of Warrants | Computed Warrant Price | Time Value');
    
    if WRT_lpxsucc
        fprintf('| Up/downside vs BBG WRT PX\n')
    else
        fprintf('\n');
    end
    
    for k = 1:length(MV)
        M = MV(k);
        [call,~] = blsprice(A/N,K,r,T,sigma);  % Running Black-Scholes
       
        % Normalizing the prices to reflect dilution after warrants are
        % exercised.
        q = M/(N*10^6);
        call = call/(1+q);                     
        TimeVal = call - WRT_intrval;
    
        % printing more results
        fprintf('%i\t |',M);
        fprintf('\t$%.4f\t|',call);
        fprintf('\t$%.4f',TimeVal);
        if WRT_lpxsucc
            perupdownside = percalculate(WRT_lpx, call);
            fprintf('\t|\t%%%.4f\n',perupdownside);
        else
            fprintf('\n');
        end
    end
    

   fprintf('\n--------------\n');

    if MAm
        SIGMAmat = zeros(MAm+1,1);
        [DRs, MEANs, SIGMAs] = sigmacalc(price(1:dind(1)-1));
        [daily_returnss, dr_means, dr_sigmas] = drsigmacalc(price(1:dind(1)-1));
        [log_daily_returnss, ldr_means, ldr_sigmas] = ldrsigmacalc(price(1:dind(1)-1));
        sigma = SIGMAs;
        SIGMAmat(1) = sigma;
        for k = 2:(MAm - 1)
            dindd = dind(k);
            dinddd = dind(k+1);
            [DR, MEAN, SIGMA] = sigmacalc(price(dindd:dindd-1));
            [daily_returns, dr_mean, dr_sigma] = drsigmacalc(price(dindd:dindd-1));
            [log_daily_returns, ldr_mean, ldr_sigma] = ldrsigmacalc(price(dindd:dindd-1));
            sigma = SIGMA;
            SIGMAmat(k) = sigma;
        end
        [DRe, MEANe, SIGMAe] = sigmacalc(price(dind(end):end));
        [daily_returnse, dr_meane, dr_sigmae] = drsigmacalc(price(dind(end):end));
        [log_daily_returnse, ldr_meane, ldr_sigmae] = ldrsigmacalc(price(dind(end):end));
        sigma = SIGMAe;
        SIGMAmat(end) = sigma; 

        for k = 1:(MAm+1)
            fprintf('SIGMAmat(%i) is %%%.4f\n',k,SIGMAmat(k)*100);
        end
        fprintf('\nComputing warrant prices using new volatilities:\n\n');

        %% Computing the bare call price for the warrant using Black Scholes.
        for k = 1:(MAm+1)
            sigma = SIGMAmat(k);
            [call,~] = blsprice(S,K,r,T,sigma);
            
            fprintf('For inputs:\n');
            fprintf('Stock price S = $%.2f\n',S);
            fprintf('Strike price K = $%.2f\n',K);
            fprintf('Annualized risk free rate R = %%%.4f\n',r*100);
            fprintf('Time T = %.4f years\n',T);
            fprintf('Annualized volality sigma = %%%.4f\n\n',sigma*100);
            fprintf('Black-Scholes MATLAB calculated call warrant price is $%.2f.\n',call);
            if WRT_lpxsucc
                percalculate = @(x,y) (y./x-1)*100;
                perupdownside = percalculate(WRT_lpx, call);
                fprintf(['The upside/downside in the computed warrant prices\nvs. ',...
                    'current market prices is %%%.4f.\n'],perupdownside);
            end
            TimeVal = call - WRT_intrval;
            fprintf('The intrinsic/time value of the warrant is $%.4f/$%.4f.\n',...
                WRT_intrval,TimeVal);
            
            %% Computing the associated (European) Greek variables.
        
            delta = 0; theta = 0; gamma = 0;  vega = 0; rho = 0; 
            
            d1 = D1(S,K,r,sigma,T);
            d2 = D2(d1,sigma,T);
            
            delta = DELTA(d1);
            theta = THETA(S,K,r,sigma,T,d1,d2);
            gamma = GAMMA(S,sigma,T,d1);
            vega = VEGA(S,T,d1);
            rho = RHO(K,r,T,d2);
        
            fprintf('Greek value MATLAB calculations for when options are European style:\n');
            fprintf('DELTA = %.4f\n',delta);
            fprintf('THETA = %.4f\n',theta);
            fprintf('GAMMA = %.4f\n',gamma);
            fprintf('VEGA = %.4f\n',vega);
            fprintf('RHO = %.4f\n',rho);
        
            fprintf('\n');
        end



    end 
end

% stocks = {'CRHC US Equity'; 'NGAB US Equity';'FRON US Equity';...
%         'HUGS US Equity';'IIAC US Equity';'THMA US Equity'};
% warrants = {'CRHC/WS US Equity'; 'NGAB/WS US Equity';'FRONW US Equity';...
%         'HUGS/WS US Equity';'IIAC/WS US Equity';'THMAW US Equity'};
% AssetA = [828.9,414.5,231.2,289.1,403.9,277.8];
% StrikeK = [11.5,11.5,11.5,11.5,11.5,11.5];
% DATES = {'10/30/2020';'03/05/2021';'05/03/2021';'04/19/2021';'01/11/2021';...
%     '03/31/2021'};
% 
% testind = length(stocks);
% testind = 0;
% for p = 1:testind
%     
%     
%     
%     % Assigning the BSC parameters
%     security = stocks{p};     % Security Ticker
%     warrant = warrants{p};
%     A = AssetA(p);            % Total cash offer in millions
%     K = StrikeK(p);            % Warrant Strike Price
%     
%     % Assigning FLDS for stock, warrant import from BBG
%     start_date = DATES{p};
%     
%     if start_date(2) =='/'
%         start_date = strcat('0',start_date);
%     end
%     
%     testdate = start_date;
%     testdate = replace(testdate,'/','');
%     testdate = replace(testdate,' ','');
%     if floor(log10(str2num('testdate')))<6
%         continue;
%     end
%     
%     end_date = '09/28/2021';
% 
%     if datetime(start_date)==datetime(end_date)
%         continue;
%     end
%     stock_data_req = {'LAST_PRICE';  'CURRENT_SHARES_OUTSTANDING_RT';...
%         'WRT_EXPIRE_DT';'SECOND_UNDERLYING_TICKER'};
%     war_data_req = {'OPT_STRIKE_PX'; 'STRIKE_PRICE_REALTIME'; 'WRT_DAYS_EXPIRE';...
%      'OPT_DAYS_EXPIRE';'LAST_PRICE';'WRT_IMPLIED_VOLATILITY_LAST';...
%      'DERIVATIVE_STRIKE_PRICE'};
% 
%     %% Importing values from BBG through API call
% 
%     c = blp;
% 
%     % Getting bulk data on the chosen security, warrant
%     [stock_data,~] = getdata(c,security,stock_data_req);
%     %warrant = strcat(stock_data.SECOND_UNDERLYING_TICKER{1},' Equity');
%     
%     [warrant_data,~] = getdata(c,warrant,war_data_req);
%     % Retrieving historical end of day price data for SLCR based on above variables
%     [Price, sec] = history(c,security, 'LAST_PRICE',...
%         start_date, end_date);
%     
%     % closing blp
%     close(c);
% 
%     price = Price(:,2);
% 
%     %% Assigning struct values from the stock/warrant structs to variables
% 
%     S = stock_data.LAST_PRICE;                           % Stock price assignment
%     N = stock_data.CURRENT_SHARES_OUTSTANDING_RT;        % Number of SLCR shares outstanding
% 
%     % Importing Warrant Struct Vals
%     
%     OPT_K  = warrant_data.OPT_STRIKE_PX; 
%     OPT_KRT = warrant_data.STRIKE_PRICE_REALTIME;
%     DER_K = warrant_data.DERIVATIVE_STRIKE_PRICE;
%     WDE = warrant_data.WRT_DAYS_EXPIRE;
%     ODE = warrant_data.OPT_DAYS_EXPIRE;
%     WRT_lpx = warrant_data.LAST_PRICE;
%     WRT_imvol = warrant_data.WRT_IMPLIED_VOLATILITY_LAST;
% 
%     %% Setting up other BSC parameters
% 
%     T = 0;
%     
%     if isnumeric(WDE)&&WDE>=0
%         T = WDE;       % Warrant strike price
%         WDEsucc = 1;
%     elseif isnumeric(ODE)&&ODE>=0
%         warning('Non-valid value for WRT_DAYS_EXPIRE. Using OPT_DAYS_EXPIRE instead.')
%         T = ODE;
%         ODEsucc = 1;
%     end
%     
%     T = T/365.25;     % Years until warrant expiration
%     
%     if ~(WDEsucc||ODEsucc)
%         T = 3;
%         warning('Warrant expiration data missing, using default value (3 years) for maturity.');
%     end
% 
%     r = [1.50];             % Annualized risk free interest rate
%     r = r/100;              % Bl-Scholes function requires decimal form for r 
%     
%     MV = ([100000:10000:300000])'; 
%     %MV = 10^5;
%     
% %     % Computing annualized volatility for SLCR based on 
% %     sq252 = sqrt(252);
% %     
% %     First compute daily returns percentages.
% %     daily_returns = 100*((price(2:end) - price(1:end-1))./price(1:end-1));
% %     log_daily_returns = log(price(2:end)./price(1:end-1));
% %     
% %     dr_mean = mean (daily_returns);
% %     ldr_mean = mean(log_daily_returns);
% %     
% %     dr_sigma = mystdcalc(daily_returns,dr_mean)*sq252;
% %     ldr_sigma = mystdcalc(log_daily_returns,ldr_mean)*sq252*100;
% %     
% %     sigma = std(daily_returns)*sq252;  % Scale by sqrt(252) to annualize sigma
% %     sigma = dr_sigma;
% %     sigma = sigma/100;                  % Bl-Scholes function requires decimal form for sigma 
% % 
% %     [DR, MEAN, SIGMA] = sigmacalc(price);
% %     [daily_returns, dr_mean, dr_sigma] = drsigmacalc(price);
% %     [log_daily_returns, ldr_mean, ldr_sigma] = ldrsigmacalc(price);
% 
%     sigma = SIGMA;
% 
%   %% Printing values from the stock/warrant structs and function parameters used
%     
%     %fprintf('Security Name: %s\n',SPAC{SPACind,2});
%     fprintf('\n\nSecurity Ticker: %s\n', security);
%     fprintf('Warrant Ticker: %s\n', warrant);
%     if isnumeric(WRT_lpx)&&(WRT_lpx>0)
%         fprintf('Warrant Price (in dollars): %.2f\n', WRT_lpx);
%     end
%     fprintf('Stock Price: %.2f\n',S);
%     fprintf('Cash Offer (in millions): %i\n',A);
%     fprintf('Number of shares outstanding (in millions): %i\n', N);
%     fprintf('SPAC Asset Value/share: %.4f\n',A/N);
%     fprintf('Warrant Strike Price: %.2f\n',K);
%     fprintf('Time until maturity expiration: %i days | %.2f years\n',T*365.25,T);
%     fprintf('Annualized Volatility (percentage): %.4f\n', sigma*100);
%     fprintf('Number of trading day prices used to compute volatility: %i\n',length(price));
%     fprintf('Annualized risk free rate (percentage): %.4f\n',r*100);
%     fprintf('SPAC IPO Date: %s\n',start_date);
% 
%     % Printing Stock Struct Values
%     format bank
%     fprintf('\nStock Struct Parameters:')
%     fprintf('\n');
%     stock_data_req = fieldnames(stock_data);
%     L = length(stock_data_req);
%     for j = 1:L
%         c = eval(['stock_data.',stock_data_req{j}]);
%         if isnumeric(c)
%             fprintf([stock_data_req{j},': %.2f\n'],c);
%         end
%     end
%     
%     fprintf('\n');
%     % Printing Warrant Struct Values
%     fprintf('Warrant Struct Parameters:')
%     fprintf('\n');
%     war_data_req = fieldnames(warrant_data);
%     L = length(war_data_req);
%     for j = 1:L
%         c = eval(['warrant_data.',war_data_req{j}]);
%         if isnumeric(c)
%             fprintf([war_data_req{j},': %.2f\n'],c);
%         end
%     end
%     
%     fprintf('\n');    
%     fprintf('M = # of Warrants Exercised\tComputed (Call) Warrant Price\n');
%     
%     for k = 1:length(MV)
%         M = MV(k);
%         [call,put] = blsprice(A/N,K,r,T,sigma);  % Running Black-Scholes       
%         % Normalizing the prices to reflect dilution after warrants are
%         % exercised.
%         q = M/(N*10^6);
%         call = call/(1+q);                     
%         put = put/(1+q);
%         
%         % printing more results
%         fprintf('%i\t\t\t\t',M);
%         fprintf('$%.4f\n',call);
%     end
% 
% end

function NP = normcdfprime(x)
    NP = exp(-x^2/2)/sqrt(2*pi);
end

function d1 = D1(S,K,r,sigma,T)
    d1 = log(S/K) + (r+sigma^2/2)*T;
    d1 = d1/(sigma*sqrt(T));
end

function d2 = D2(d1,sigma,T)
    d2 = d1 - sigma*sqrt(T);
end

function delta = DELTA(d1) 
    delta = normcdf(d1);
end

function theta = THETA(S,K,r,sigma,T,d1,d2)
    theta = -S*normcdfprime(d1)*sigma/(2*sqrt(T));
    theta = theta - r*K*exp(-r*T)*normcdf(d2);
end

function gamma = GAMMA(S,sigma,T,d1)
    gamma = normcdfprime(d1)/(S*sigma*sqrt(T));
end

function vega = VEGA(S,T,d1)
    vega = S*sqrt(T)*normcdfprime(d1);
end

function rho = RHO(K,r,T,d2)
    rho = -K*T*exp(-r*T)*normcdf(-d2);
end

function avg = IMVOL(mybls,callval,tol)
    a = rand(1); b = rand(1); 
    callhat = 0; mcounter = 10^3; 
    while mybls(a)>callval
        a = a/2;
    end
    while mybls(b)<callval
        b = b*2;
    end
    M = max(a,b); m = min(a,b); avg = mean([m,M]); counter = 1;   
    while true
        callhat = mybls(avg); diff = callhat - callval;
        if diff>tol
            M = avg;
        elseif diff<-tol
            m = avg;
        else
            break;
        end
        avg = mean([m,M]); counter = counter + 1;
        if counter>mcounter
            warning('Implied volatility calculation did not converge.');
            break;
        end
    end
end

function stdval = mystdcalc(u,ubar)
    stdval = sqrt((sum((u-ubar).^2)/(length(u) - 1)));
end

function [DR, MEAN, SIGMA] = sigmacalc(PRICE)
    DR = ((PRICE(2:end) - PRICE(1:end-1))./PRICE(1:end-1));
    MEAN = mean(DR);
    SIGMA = std(DR)*sqrt(252);
end

function [daily_returns, dr_mean, dr_sigma] = drsigmacalc(PRICE)
    daily_returns = ((PRICE(2:end) - PRICE(1:end-1))./PRICE(1:end-1));
    dr_mean = mean(daily_returns);
    dr_sigma = mystdcalc(daily_returns,dr_mean)*sqrt(252);
end

function [log_daily_returns, ldr_mean, ldr_sigma] = ldrsigmacalc(PRICE)
    log_daily_returns = log(PRICE(2:end)./PRICE(1:end-1));
    ldr_mean = mean(log_daily_returns);
    ldr_sigma = mystdcalc(log_daily_returns,ldr_mean)*sqrt(252);
end

% Recording the column names from the worksheet
% ColHeads = cell(n,1);
% for k = 1:n
%     ColHeads{k} = SPAC{1,k};
% end

% Finding indices of the SPAC equities, eliminating noise, and 
% applying BBG naming conventions.

% Working with a more manageable SPAC cell, but may not need this.;

% SPACindices = [1, 2, 5, 7, 12:14, 18:21, 23];
% coln = length(SPACindices);
% SPACcell = cell(m,coln);
% k = 0;

% for j = 1:coln
%     k = SPACindices(j);
%     SPACcell(:,j) = SPAC(:,k);
% end


    % trouble_date = '#N/A N/A';
    %     if all(size(start_date)==size(trouble_date))
    %         if start_date == '#N/A N/A'
    %             if API_SDsucc
    %                 start_date = API_start_date;
    %             else
    %                 warning('Invalid start date entry %s for SPAC %s or IPO date unavailable.\n',...
    %                     trouble_date, security);
    %                 continue;
    %             end
    %         end
    %     end
    
    %  if isempty(start_date)
%         if ~isempty(API_start_date)
%             start_date = API_start_date;
%         else
%             warning('SPAC %s IPO date unavailable.\n',security);
%         end
%     end
% 
%     testdate = start_date;
%     testdate = replace(testdate,'/','');
%     testdate = replace(testdate,' ','');
%    
%     if floor(log10(str2num('testdate')))<6
%         warning('Invalid start date entry form for SPAC %s\n',security);
%         continue;
%     end

%     if ~(WDEsucc||ODEsucc)
%         T = 3;
%         warning('Warrant expiration data missing, using default value (3 years) for maturity.');
%     end

%     %% Computing annualized volatility for SLCR based on 
%     mystdcalc = @(u,ubar) sqrt((sum((u-ubar).^2)/(length(u) - 1)));
%     sq252 = sqrt(252);
%     
%     % First compute daily returns percentages.
%     daily_returns = 100*((price(2:end) - price(1:end-1))./price(1:end-1));
%     log_daily_returns = log(price(2:end)./price(1:end-1));
%     
%     dr_mean = mean (daily_returns);
%     ldr_mean = mean(log_daily_returns);
%     
%     dr_sigma = mystdcalc(daily_returns,dr_mean)*sq252;
%     ldr_sigma = mystdcalc(log_daily_returns,ldr_mean)*sq252*100;
%     
%     %sigma = std(daily_returns)*sq252;  % Scale by sqrt(252) to annualize sigma
%     sigma = dr_sigma;
%     sigma = sigma/100;                  % Bl-Scholes function requires decimal form for sigma 
    