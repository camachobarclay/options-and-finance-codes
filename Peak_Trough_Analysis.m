%%
function Peak_Trough_analysis(security,BENCH,start_date,end_date,tdw,trigger,fdwMAT,elim_crit,incBENCH)

javaaddpath('C:\blp\BloombergWindowsSDK\JavaAPI\v3.15.1.1\lib\blpapi3.jar')

% Starting the BLP session to call for data

c = blp;

% Getting historical price data for securities of interest

[price, sec] = history(c,security, 'LAST_PRICE',...
    start_date, end_date);
[BENCHprice, BENCHsec] = history(c,BENCH, 'LAST_PRICE',...
    start_date, end_date);

% Close the BLP session

close(c);

%%

% Renaming the security character variables so the non-ticker component is
% dropped.

c = zeros(1,10); 

for k = 1:length(security)
    if security(k)~=' '
        c(k) = security(k);
    else
        break;
    end
end


security = char(nonzeros(c)');


c = zeros(1,10);   

for k = 1:length(BENCH)
    if BENCH(k)~=' '
        c(k) = BENCH(k);
    else
        break;
    end
end

BENCH = char(nonzeros(c)');

%%

% Get the dates and prices for the security in acceptable form

dateS = datestr(price(:,1));            % The original char array of dates for the security
dates = string(dateS);                  % Now in string form.
price = price(:,2);                     % Keeping the stock prices for the security

BENCHdateS = datestr(price(:,1));       % Do the same thing for the beenchmark security.
BENCHdates = string(BENCHdateS);
BENCHprice = BENCHprice(:,2);

% Get the lengths of the securities and make sure it is the same number of
% trading days

n = length(price);
BENCHn = length(BENCHprice);

if n~=BENCHn
%     warning(['Difference in trading days for ',BENCH,' and ',security,'.']);
%     warning([BENCH,'Ignored, Peak_Trough_analysis will produce results for ',security,' only.']);
%     incBENCH = 0;
end

% BENCH = 'QQQ US Equity';
% start_date = '09/24/2015';
% end_date = '09/24/2020';
% tdw = 15;
% trigger = 10;
% fdwMAT = 5:5:15;

%%

j = tdw-1;                              % Use this for indexing into.
tdw_max_prices = zeros(n-j,1);          % This will record the rolling max prices.
tdw_pctchng = tdw_max_prices;           % This will record the drawdown decline from the high.

new_dates = dates(tdw:end);             % Working with the newly indexed dates that skip through to 
new_prices = price(tdw:end);            % Newly indexed prices

BENCH_new_dates = BENCHdates(tdw:end);
BENCH_new_prices = BENCHprice(tdw:end);

pctchngcalc = @(x,y)(y/x - 1)*100;      % Use to compute % changes

fdwM = max(fdwMAT);                     % Max day for post session change record
%elim_crit = fdwM;                      % Elimination criteria for clusters
M = n-(j+fdwM);                         % M is the number of days that we will compute statistics for.                                        
                                        % Note M < n - fdw = number of days with recorded decline < n = number of recorded days.
fdwL = length(fdwMAT);                  % Number of sessions we are considering.
reactMAT = zeros(M,fdwL*2);             % Where we record the % session changes
ind = (1:M)';                           % Index vector

mylocs = zeros(M,1);                    % Where relevant locations are stored
a = (1:fdwL);                           % Use a and b to interlace indices for security/BENCH comparisons.
b = (fdwL+1:2*fdwL);

reactIND = reshape([a;b],size(a,1),[])';

%% Compute rolling twd max price, percent drawdown, store results in vectors.

for k = tdw:n
    tdw_max_prices(k-j) = max(price(k-j:k));                    % Record the max price
    cur_max_price = tdw_max_prices(k-j);                        % Helper variables
    cur_price = price(k);
    tdw_pctchng(k-j) = pctchngcalc(cur_max_price,cur_price);    % Use helper variables to compute the percent change decline
end

%%

loccomp = find(abs(tdw_pctchng(1:end-fdwM))>trigger);
L = length(loccomp);

% Have to output only specific days
% decide how much to outpace by

cur_ind =  min(loccomp);


for i = 1:M
    if (abs(tdw_pctchng(i))>trigger)
        if i==cur_ind ||(i>cur_ind && (i - cur_ind)>elim_crit)
            cur_ind = i;
            mylocs(i) = i;
            reactMAT(i,1:fdwL) = pctchngcalc(price(i+j),price(i+j+fdwMAT));
            reactMAT(i,fdwL+1:2*fdwL) = pctchngcalc(BENCHprice(i+j),BENCHprice(i+j+fdwMAT));
        end
    end
end

reactMAT = reactMAT(:,reactIND);
mylocs = nonzeros(mylocs);

BigArray = [ind price(tdw:n - fdwM) tdw_max_prices(1:end-fdwM) tdw_pctchng(1:end-fdwM) reactMAT];
locs = find(reactMAT(:,1));

% normlocs = norm(mylocs - loccomp,1)
%%

format

% fprintf('Norm test is %f\n\n',norm(locs - loccomp,1));

T = table(ind(locs), datetime(new_dates(locs)), new_prices(locs), tdw_max_prices(locs),...
    tdw_pctchng(locs),'VariableNames',{'Index','Date','Close PX','Prev roll. high','% High Dec.'});
    beat = zeros(size(locs));
    secwin = 0;
    secloss = 0;
    lbeat = length(beat);
AT = table([1:length(locs)]','VariableNames',{'Index'});

fprintf('\n\n-----------------------------\n\n');
fprintf('Starting algorithm Peak_Trough_analysis.m\n\n');
fprintf(['Security = ',security,'\n']);
fprintf(['Benchmark Comparison = ',BENCH,'\n']);
fprintf(['Start Date = ',start_date,'\n']);
fprintf(['End Date = ',end_date,'\n']);
fprintf('Threshold trigger = %f\n',trigger);
fprintf('Trailing day window of %d trading sessions.\n', tdw);
fprintf('Post reaction days');

for h = 1:fdwL
    fprintf(' +%d',fdwMAT(h));
    if h<fdwL
        fprintf(',');
    else
        fprintf('\n');
    end
end

color1 = rand(1,3);
color2 = rand(1,3);
fprintf('Trigger events within %d day period of initial trigger are ignored.\n\n',elim_crit);

fprintf('-----------------------------\n\n');

%%

for l=1:fdwL

    col1 = reactMAT(locs,l*2-1);
    col2 = reactMAT(locs,l*2);
    
    beat = col1>col2 + 0.5*(col1==col2);
    col1sorted = sort(col1,'descend');
    col2sorted = sort(col2,'descend');
    
    minval1 = min(col1);
    maxval1 = max(col1); 
    
    minval2 = min(col2);
    maxval2 = max(col2); 

    range1 = range(col1sorted);
    range2 = range(col2sorted);
    secwin = sum(beat);
    secloss = lbeat - secwin;
    
    dval = char(num2str(fdwMAT(l)));
    
    T = addvars(T,col1,'NewVariableNames',char([security,' + ',dval,' d % chng ']));
    T = addvars(T,col1>col2 + 0.5*(col1==col2),'NewVariableNames',char([security, ' beat ', BENCH,' + ',dval, ' d']));
    T = addvars(T,col2,'NewVariableNames',char([BENCH,' + ',dval,' d % chng']));
    
    AT = addvars(AT,col1sorted,'NewVariableNames',char([security,' + ',dval,' d % chng sorted']));
    AT = addvars(AT,col2sorted,'NewVariableNames',char([BENCH,' + ',dval,' d % chng sorted']));

    fprintf([security,' + ',dval,' d %% change data has min value %f, max value %f, and range %f.\n'], minval1, maxval1, range1);
    fprintf([BENCH,' + ',dval,' d %% change data has min value  %f, max value %f, and range %f.\n\n'], minval2, maxval2, range2);
    
    fprintf(['For + ',dval,' sessions ']);
    
    if secwin>secloss
        fprintf([security,' beat ',BENCH,' %.1f to %.1f'], secwin,secloss);
    elseif secwin<secloss
        fprintf([BENCH,' beat ',security,' %.1f to %.1f'], secloss,secwin);
    elseif secwin==secloss
        fprintf([security,' and ',BENCH,' tied at %.1f each.'], secwin,secloss);
    end

    fprintf('\n\n-----------------------------\n\n');
    
 %%   
    figure;
    x = categorical(string(datestr(new_dates(locs),'yyyy/mm/dd'))); 
    b = bar(x,[col1,col2],1);
    set(b(1),'FaceColor',color1,'EdgeColor',color1);
    set(b(2),'FaceColor',color1*.5,'EdgeColor',color1*.5);

    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels1 = string(round(b(1).YData,2));
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
    'VerticalAlignment','top')
    
    xtips2 = b(2).XEndPoints;
    ytips2 = b(2).YEndPoints;
    labels2 = string(round(b(2).YData,2));
    text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')
    set(gcf, 'Position',  [0, 0, 1550, 650])
%     set(gcf, 'Position',  [0, 0, 2000, 1000])

    
    grid on;
    hold on;
   
    legend(security,BENCH,'Location','southeast','NumColumns',2);
    xlabel('Time');
    ylabel('% Reaction');
    title([security,', Benchmark ', BENCH,' % Change Reaction after + ',char(num2str(fdwMAT(l))),' Day Session']);

    hold off
    
end

boxcat = cell(1,fdwL);

for k = 1:fdwL
    boxcat{k} = ['+',num2str(fdwMAT(k)),' d'];
end

%%

figure;
% mybox = boxplot(reactMAT(locs,1:2:end-1),'Labels',boxcat, 'Whisker', inf);
mybox = boxplot(reactMAT(locs,1:2:end-1),'Labels',boxcat);

xlabel('column');
ylabel('% range');
title([security,' Boxplot']);
 set(gcf, 'Position',  [0, 0, 1550, 650])
% set(gcf, 'Position',  [0, 0, 2000, 1000])

grid on

%%

disp(T)
disp(AT)

fprintf(['\n',security,' individual column median values:\n'])
medianMAT = median(reactMAT(locs,1:2:end-1));
disp(medianMAT)

fprintf(['\n',security,' individual column mean values:\n'])
meanMAT = mean(reactMAT(locs,1:2:end-1));
disp(meanMAT);

fprintf('\n-----------------------------\n\n');

diary off;

end

