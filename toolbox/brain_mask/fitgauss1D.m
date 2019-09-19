function hDat = fitgauss1D(Data1D, ...
    fittype, fSigma, Idx, bin, xif, yif)
% fitgauss1D: fit gauss distributions to 1D data
%
% Usage:
%   fitgauss_to_fluohist(wDat, ...
%       lnbis, gaussnumber, sigmaRed, ...
%       gauss2use, image2replace, figH, ...
%       axH, plot_flag)
%
% Args:
%   Data1D: vector to fit
%   fittype: number of gaussian to fit (gauss1, gauss2, etc) 
%   fSigma: sigma
%   Idx: gaussian number to retrieve
%   bin: bining interval
%   xif: input x-variable of histogram
%   yif: input y-variable of histogram

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
            
            fprintf('Error fSigma has more than 2 values\n');
            return
            
        end
        
        hDat.xths = repmat(hDat.xths,[numel(hDat.yths) 1 1]);
        
    end
    
end

end
