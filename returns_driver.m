   clear;
close all;
clc;

javaaddpath('C:\blp\BloombergWindowsSDK\JavaAPI\v3.15.1.1\lib\blpapi3.jar')

% SECMAT  = {'AAPL US Equity', 'MSFT US Equity', 'AMZN US Equity',...
%             'GOOG US Equity', 'FB US Equity', 'NFLX US Equity',};
%             'ADBE US Equity', 'CRM US Equity', 'NOW US Equity',...
%             'TSLA US Equity','NVDA US Equity', 'AMD US Equity',...
%             'SHOP US Equity', 'PYPL US Equity', 'SQ US Equity',...
%             'MA US Equity', 'V US Equity', 'JPM US Equity'};
 SECMAT ={'IBM US Equity'};
% SECMAT = {'COUP US Equity', 'CRWD US Equity', 'DDOG US Equity',...
%     'DOCU US Equity', 'ESTC US Equity', 'FSLY US Equity',...
%     'MDB US Equity', 'OKTA US Equity', 'NET US Equity',...
%     'NOW US Equity', 'TEAM US Equity', 'TWLO US Equity',...
%     'ZS US Equity', 'ZM US Equity'};
% START_DATES = {'03/07/2014', '06/11/2019', '09/18/2019',...
%     '04/26/2018','10/04/2018','05/16/2019','10/18/2017',...
%     '04/07/2017','09/12/2019','06/29/2012','12/11/2015',...
%     '06/24/2016','03/16/2018','04/19/2019'};
BENCH = 'QQQ US Equity';
start_date = '07/20/2011'; 
end_date = '07/20/2021';
tdw = 10;
trigger = 5;
fdwMAT = [1:15];
elim_crit = max(fdwMAT);
incBENCH = 1;
%%
for k = 1:length(SECMAT)
    security = SECMAT{k};
%     start_date = START_DATES{k};
        Peak_Trough_analysis(security,BENCH,start_date,end_date,tdw,trigger,fdwMAT,elim_crit,incBENCH);
        snapnow;
end

%%



% plot(datetime(new_dates(1111:1195)), new_prices(1111:1195))
% hold on
% scatter(datetime(new_dates(locs(find(locs==1111):find(locs==1132)))),new_prices(locs(find(locs==1111):find(locs==1132))),'ro')








% new_dateS = datetime(char(dates{tdw:end}));

% [price, dates, raw] = xlsread('FB_Prices.xlsx');
% Dates = datetime(char(dates{:}));


% fprintf('Date\t\t\tPrice\tCurrent Recent Max\t%% Current Change\t')
% 
% for l = 1:fdwL
%    fprintf('+%i Session %% Change\t\t',fdwMAT(l)) 
% end
% 
% fprintf('\n');
% 
% for m = 1:length(locs)
%     k = locs(m);
% %     fprintf('%d\t',ind(k))
%     fprintf('%s\t\t',string((new_dates(k))))
%     fprintf('$%.2f\t$%.2f\t\t\t\t%.2f%%\t\t\t\t',price(k+tdw), tdw_max_prices(k), tdw_pctchng(k)); 
%     for l = 1:fdwL
%        fprintf('%.2f%%\t\t\t\t\t',reactMAT(k,l)); 
%     end
%     fprintf('\n');
% end
% 
% 
% % 
% % [pks, pklocs, pw, pp] = findpeaks(price,dates);                               % Original Peaks
% % [troughs,trlocs, tw, tp] = findpeaks(-price,dates);
% % troughs = -troughs;
% % 
% % 
% % plot(dates,price)
% % 
% % hold on
% % 
% % pause;
% % 
% % for k = 1:length(pklocs)
% %     pause(1/5);
% %     scatter(pklocs(k), pks(k),'r')
% % end
% % 
% % pause; 
% % 
% % for k = 1:length(trlocs)
% %     pause(1/5);
% %     scatter(trlocs(k),troughs(k),'b')
% % end
% % 
% % if trlocs(1)<pklocs(1)
% %     trlocs = trlocs(2:end);
% %     troughs = troughs(2:end);
% % end
% % 
% % for k = 1:min(length(trlocs),length(pklocs))
% %     
% %     
% %     
% % end
% %     
% 

