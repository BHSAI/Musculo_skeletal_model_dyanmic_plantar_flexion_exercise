function [a,b,c,rsq,curve]=classic_hill(x,y)
fo=fitoptions('Method','NonLinearLeastSquares','lower',[0,0,0],'upper',[Inf,Inf,Inf],'StartPoint',[1,1,1]);
ft=fittype('(a*(x^b))/(x+c)','options',fo);
[curve,gof] = fit(x,y,ft);
a=curve.a;
b=curve.b;
c=curve.c;
rsq=gof.rsquare;
end