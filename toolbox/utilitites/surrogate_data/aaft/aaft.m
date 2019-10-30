function [zV] = aaft(xV, nsur, xV2)
% [zV] = aaft(xV, nsur, xV2)
% AAFT (Amplitude Adjusted Fourier Transform)
% For more details see: D. Kugiumtzis, “Surrogate data test for nonlinearity including monotonic transformations,” Phys. Rev. E, vol. 62, no. 1, 2000.
% generates surrogate data for a given time series 'xV' using the AAFT  
% 
% INPUT:
% - xV        : The original time series (column vector)
% - nsur      : The number of surrogate data to be generated
% - xV2       : Additional timeseries (binary) that limits 
%               arragements to use (useful to avoid similar stimuli arragements)
% OUTPUT:
% - wM        : The CAAFT surrogate data (matrix of 'nsur' columns)
%
% modified 03/10/2018, added a extra timeseries for further controlling
% random phase arragements (useful to avoid same stimuli structure)

n = length(xV);

% The following gives the rank order, ixV
[oxV, T1] = sort(xV);
[T, ixV] = sort(T1);

% ===== AAFT algorithm 
zV = zeros(n, nsur); 
k = 1;

% plot magnitude and angles
%axH(1) = subplot(2, 1, 1);
%axH(2) = subplot(2, 1, 2);

while k <= nsur
    
    % Rank ordering white noise with respect to xV
    rV = randn(n, 1);
    [orV, T] = sort(rV);
    yV = orV(ixV); % Y
    
    % >>>>> Phase randomisation (Fourier Transform): yV -> yftV 
    if rem(n, 2) == 0
      n2 = n/2;
    else
      n2 = (n-1)/2;
    end
    tmpV = fft(yV, 2*n2);
    magnV = abs(tmpV);
    fiV = angle(tmpV);
    
    % plot magnitude and angles
    %plot(magnV, 'Parent', axH(1)); hold on
    %plot(fiV, 'Parent', axH(2)); hold on
    
    rfiV = rand(n2-1, 1) * 2 * pi;
    nfiV = [0; rfiV; fiV(n2+1); -flipud(rfiV)];
    
    % New Fourier transformed data with only the phase changed
    tmpV = [magnV(1:n2+1)' flipud(magnV(2:n2))']';
    tmpV = tmpV .* exp(nfiV .* i); 
    
    % Transform back to time domain
    yftV = real(ifft(tmpV, n)); % 3-step AAFT;
    
    % <<<<<
    [T, T2] = sort(yftV); % Rank ordering xV with respect to yftV 
    [T, iyftV] = sort(T2);
    
    if exist('xV2', 'var') && ~isempty(xV2)
        if ~ismember(xV2, xV2(iyftV), 'rows') 
            zV(:, k) = oxV(iyftV);  % zV is the AAFT surrogate of xV
            k = k + 1;
        end
    else
        zV(:, k) = oxV(iyftV);  % zV is the AAFT surrogate of xV
        k = k + 1;
    end
    
end

end
