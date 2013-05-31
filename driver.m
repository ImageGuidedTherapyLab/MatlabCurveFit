clear all
close all
load test.mat

options = optimset('display','iter','Jacobian','on');
ydata = double(image(50:100,50:100,:));
xdata = 40.*[1:size(ydata,3)]';
% allocate space
%   array shape wrt jacobian of matrix function and variables, see:
% http://www.mathworks.com/help/optim/ug/writing-objective-functions.html#brkjub4
x0 = rand(3,size(ydata,1),size(ydata,2));
UB = inf(3,size(ydata,1),size(ydata,2));
LB = -inf(3,size(ydata,1),size(ydata,2));
[x,resnorm,residual,exitflag, output]=lsqcurvefit(@t2decay,x0,xdata,ydata,LB,UB,options);
%[x,resnorm,residual,exitflag, output]=lsqcurvefit(@t2decay,x0,xdata,ydata);
