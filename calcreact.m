function calcreact(security, priorclose, postmarket, premarket, open, close, commentary,valuation)

percentcalc = @(x,y) 100*(y/x - 1);
fprintf('%s US Equity\n',security);
fprintf('Current Valuation = %s\n\n',valuation);
initpostreact = percentcalc(priorclose,postmarket);
overnightmove = percentcalc(postmarket,premarket);
nearopenmove = percentcalc(premarket,open);
endofdaymove = percentcalc(open,close);

fprintf('Initial reaction post market (4pm to 8pm): %.2f%% move from $%.2f to $%.2f\n',initpostreact, priorclose, postmarket);
fprintf('Postmarket to premarket overnight move (8pm to 6am): %.2f%% from $%.2f to $%.2f\n',overnightmove,postmarket, premarket);
fprintf('Premarket to open move (6am to 9:30am): %.2f%% from $%.2f to $%.2f\n',nearopenmove, premarket, open);
fprintf('Next day open to close move (9:30am to 4pm): %.2f%% from $%.2f to $%.2f\n', endofdaymove, open, close);
fprintf('%s\n\n', commentary)

end