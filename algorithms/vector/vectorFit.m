function solution = vectorFit(xdata, ydata)

    options = optimset('display','iter','jacobian','on','MaxIter',40,'MaxFunEvals',40,'TolFun',1e-6);
    %EchoTime = [ 1.468, 3.764, 6.06, 8.356, 10.652, 12.948, 15.244, 17.54, 19.836, 22.132, 24.428, 26.724];
    % allocate space
    %   array shape wrt jacobian of matrix function and variables, see:
    % http://www.mathworks.com/help/optim/ug/writing-objective-functions.html#brkjub4
    x0 = rand(3,size(ydata,1),size(ydata,2));
    UB =  ones(3,size(ydata,1),size(ydata,2));
    % TODO - Better bounding values
    UB(1,:,:) = 200 * UB(1,:,:) ;
    UB(2,:,:) =  60 * UB(2,:,:) ;
    UB(3,:,:) = 100 * UB(3,:,:) ;
    LB =  zeros(3,size(ydata,1),size(ydata,2));
    [solution,resnorm,residual,exitflag, output]=lsqcurvefit(@vectorT2Decay,x0,xdata,ydata,LB,UB,options);
    %[x,resnorm,residual,exitflag, output]=lsqcurvefit(@t2decay,x0,xdata,ydata);

end
