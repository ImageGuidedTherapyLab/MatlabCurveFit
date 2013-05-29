
function [ModelVector,ModelJacobian] = t2decay(x,xdata)
% objective function values at x
ModelVector = zeros(size(x,1),size(x,2),size(xdata,1));
% build array of model predicted values
for iii = 1:size(xdata,1)
 echotime = xdata(iii,1)
 ModelVector(:,:,iii) =  x(:,:,1).*exp(-echotime./x(:,:,2)) + x(:,:,3);
end
if nargout > 1   % two output arguments
   ModelJacobian = 0; % Jacobian of the function evaluated at x
end
