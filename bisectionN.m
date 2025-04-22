function c = bisectionN(a,b,f,tol)
e =tol;
itermax=100;
for(i=1:itermax)
    c=(a+b)/2;
    if((abs(f(c)))<e)
        break
    else
        if(f(a)*f(c)<0)
            b=c;
        else
            a=c;
        end
    end
end
