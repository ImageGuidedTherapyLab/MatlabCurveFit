function error = rmse( x, y )
    
    error = sqrt( sum( (x(:) - y(:)).^2 ) / numel(y) );
    
end
