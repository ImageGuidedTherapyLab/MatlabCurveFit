
clear all
close all


load test.mat

options = optimset('display','off');
ydata = double(image(50:100,50:100,:));
xdata = 40.*[1:12]';
% allocate space
%x0 = zeros(size(x,1),size(x,2),3);
x0 = rand(size(ydata,1),size(ydata,2),3);
%[x,resnorm,residual,exitflag, output]=lsqcurvefit(@t2decay,x0,xdata,ydata,LB,UB,options);
[x,resnorm,residual,exitflag, output]=lsqcurvefit(@t2decay,x0,xdata,ydata);
