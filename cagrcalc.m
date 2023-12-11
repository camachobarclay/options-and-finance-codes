function CAGR = cagrcalc(baseval, futureval, t0, tn)

elapseT = tn - t0;
CAGR =  (futureval/baseval)^(1/elapseT) - 1;

end