y = 1.9;

[call, ~] = blsprice(100, 95, 0.1,0.25,y);
clear sigma
mybl = @(x) blsprice(100,95,0.1,0.25,x);

a = rand(1);
b = rand(1);
tol = 10^(-6);
    fprintf('entering a loop\n')
while mybl(a)>call
    a = a/2;
%     fprintf('a = %f\n',a);
%     fprintf('call of a = %f\n',mybl(a));
%     fprintf('orig call val = %f\n',call);
end
    fprintf('entering b loop\n');
while mybl(b)<call
    b = 2*b;
%     fprintf('b = %f\n',b);
%     fprintf('call of b = %f\n',mybl(b));
%     fprintf('orig call val = %f\n',call);
end

M = max(a,b);
m = min(a,b);
avg = mean([m,M]);
counter = 1;

while true
    callhat = mybl(avg);
    diff = callhat - call;
    if diff>tol
        M = avg;
    elseif diff<-tol
        m = avg;
    else
        break;
    end
    avg = mean([m,M]);
    counter = counter + 1;
    
end