javaaddpath('C:\blp\BloombergWindowsSDK\JavaAPI\v3.15.1.1\lib\blpapi3.jar');
security = 'SPLK US EQUITY';

start_date = '12/03/2020';
end_date = '11/09/2021';

c = blp;

[Price, sec] = history(c,security, 'LAST_PRICE',...
        start_date, end_date);

close(c);
PT = [mean(Price(:,2)), median(Price(:,2))];