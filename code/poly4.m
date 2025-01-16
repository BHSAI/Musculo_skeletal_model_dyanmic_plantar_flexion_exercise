function [params,rsq,curve]=poly3(x,y)
    fo=fitoptions('Method','NonLinearLeastSquares','lower',[-Inf,-Inf,-Inf,-Inf,-Inf],'upper',[Inf,Inf,Inf,Inf],'StartPoint',[1,1,1,1,1]);
    ft=fittype('(a*(x^4))+(b*(x^3))+(c*(x^2))+(d*(x^1))+e','options',fo);
    [curve,gof] = fit(x,y,ft);
    a=curve.a;
    b=curve.b;
    c=curve.c;
    d=curve.d;
    e=curve.e;
    params=[a;b;c;d;e];
    rsq=gof.rsquare;
end