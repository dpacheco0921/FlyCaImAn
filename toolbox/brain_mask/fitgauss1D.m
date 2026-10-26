function hDat = fitgauss1D(Data1D, ...
    fittype, fSigma, Idx, bin, xif, yif)
% fitgauss1D: fit gauss distributions to 1D data
%
% Usage:
%   hDat = fitgauss1D(Data1D, ...
%       fittype, fSigma, Idx, bin, xif, yif)
%
% Args:
%   Data1D: vector to fit
%   fittype: number of gaussian to fit (gauss1, gauss2, etc)
%       (default, 'gauss2')
%   fSigma: number of standard-deviation-like widths (coef c of the fit)
%       around the center of gaussian Idx used to set the thresholds.
%       One value is used on both sides, two values are [left right]
%       (default, 5)
%   Idx: gaussian number to retrieve
%       (default, 1)
%   bin: bining interval, in the units of Data1D
%       (default, 1)
%   xif: input x-variable of histogram
%   yif: input y-variable of histogram
%       (xif and yif are computed from Data1D if either is empty)
%
% Returns:
%   hDat: structure with fields
%       xif, yif: histogram x values (units of Data1D) and counts
%       coef: fitted coefficients [a1 b1 c1 a2 b2 c2 ...] (a: height,
%           b: center, c: width)
%       yfit1, yfit2, ...: fitted curve of each gaussian over xif
%       yths: y values (1:max count) used to draw the threshold lines
%       Idx: gaussian number used for the thresholds
%       xths: [left right] threshold (units of Data1D) repeated for each
%           row of yths

if ~exist('fittype','var')
    fittype = 'gauss2';
end

if ~exist('fSigma','var')
    fSigma = 5;
end

if ~exist('Idx','var')
    Idx = 1;
end

if ~exist('bin','var')
    bin = 1;
end

if ~exist('xif','var')
    xif = [];
end

if ~exist('yif','var')
    yif = [];
end

Data1D = double(Data1D);
if ~isempty(xif) && ~isempty(yif)
    hDat.xif = xif;
    hDat.yif = yif;
else
    hDat.xif = 0:bin:max(Data1D)*1.2;
    [hDat.yif, ~] = hist(Data1D, hDat.xif);
    hDat.yif(1) = 0;
end

% fit n-gaussians
Fitpar = fit(double(hDat.xif'), double(hDat.yif'), fittype);
hDat.coef = coeffvalues(Fitpar);
ngauss = (numel(hDat.coef)/3);
hDat.yths = 1:max(hDat.yif);
hDat.Idx = Idx;

for ng = 1:ngauss
    
    eval(['hDat.yfit',num2str(ng), ...
        ' = hDat.coef(', num2str(1+3*(ng-1)), ...
        ')*exp(-((hDat.xif-hDat.coef(', ...
        num2str(2+3*(ng-1)),'))/hDat.coef(', ...
        num2str(3+3*(ng-1)),')).^2);']);
    
    if ng == Idx
        
        if numel(fSigma) == 1
            
        	eval(['hDat.xths = [hDat.coef(', ...
                num2str(2+3*(ng-1)),')-fSigma*hDat.coef(', ...
                num2str(3+3*(ng-1)),') hDat.coef(', ...
                num2str(2+3*(ng-1)),')+fSigma*hDat.coef(', ...
                num2str(3+3*(ng-1)),')];'])
            
        elseif numel(fSigma) == 2
            
            eval(['hDat.xths = [hDat.coef(', ....
                num2str(2+3*(ng-1)),')-fSigma(1)*hDat.coef(', ...
                num2str(3+3*(ng-1)),') hDat.coef(', ...
                num2str(2+3*(ng-1)),')+fSigma(2)*hDat.coef(', ...
                num2str(3+3*(ng-1)),')];'])
            
        else
            
            % stop here, otherwise hDat.xths would be missing for the caller
            error('fitgauss1D: fSigma has more than 2 values')
            
        end
        
        hDat.xths = repmat(hDat.xths,[numel(hDat.yths) 1 1]);
        
    end
    
end

end
