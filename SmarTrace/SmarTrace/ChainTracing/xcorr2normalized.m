function c = xcorr2normalized(x, y)
%     % method 1
%     c=normxcorr2(y,x);
% 
%     % only keep the center
%     d1=floor((size(c)-size(x))/2);
%     d2=ceil((size(c)-size(x))/2);
%     c = c(d1(1)+1:end-d2(1), d1(2)+1:end-d2(2));
    
    % see 'coeff' flag in xcorr matlab command
%     % Compute autocorrelations at zero lag
%     %cxx0 = sum(abs(x).^2);
%     cxx0 = sum(sum(abs(x).^2));
%     cyy0 = sum(sum(abs(y).^2));
%     scale = sqrt(cxx0*cyy0);
%     c = c./scale;

    % method 2
    c=xcorr2(y,x);
    c=flipud(fliplr(c));

    % only keep the center
    d1=floor((size(c)-size(x))/2);
    d2=ceil((size(c)-size(x))/2);
    c = c(d1(1)+1:end-d2(1), d1(2)+1:end-d2(2));
end

