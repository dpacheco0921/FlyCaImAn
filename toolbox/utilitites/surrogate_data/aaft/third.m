function [XXX,waxis]=third(y,nlag,centre,outmax);
% Third order moment;
%       function [xxx,waxis]=third(y,nlag,centre);
%       y       - time-series
%       nlag    - number of lags, max=n
%       centre  - centre series by subtracting mean (1=Yes, 2=No)
%       outmax  - display lags of maxima and minima (1=Yes, default=No)

    [nsamp width]=size(y);
if nsamp<10;
   disp('warning n<10')
end;

if exist('outmax')~=1;
   outmax=0;
end;

% ---------------- cumulants in non-redundant region -----------------
     XXX = zeros(nlag+nlag+1,nlag+nlag+1);
     if centre==1;
        mn=mean(y);
        centred=y-mn; 
     else
        centred=y; 
     end;
     for d=0:nlag;
         for k=d:nlag;
             large=max([d;k;0]);
             XXX(d+nlag+1,k+nlag+1)=sum(centred(1:nsamp-large).*centred(1+d:nsamp-large+d)...
                                     .*centred(1+k:nsamp-large+k))/nsamp;
             XXX(nlag+1+k,nlag+1+d)=XXX(d+nlag+1,k+nlag+1); % Symmetry
             XXX(nlag+1-d,nlag+1+k-d)=XXX(d+nlag+1,k+nlag+1); % Symmetry
             XXX(nlag+1+k-d,nlag+1-d)=XXX(d+nlag+1,k+nlag+1); % Symmetry
             XXX(nlag+1+d-k,nlag+1-k)=XXX(d+nlag+1,k+nlag+1); % Symmetry
             XXX(nlag+1-k,nlag+1+d-k)=XXX(d+nlag+1,k+nlag+1); % Symmetry
         end;
     end;
  
     waxis=(-nlag:nlag);

% Lags of minima and maxima
if outmax==1;
   [maxmaxr,maxmaxc]=find(XXX==max(max(XXX)));
   [minminr,minminc]=find(XXX==min(min(XXX)));
   disp('Maximum at');
   output=cat(2,maxmaxr(1)-nlag-1,maxmaxc(1)-nlag-1);
   disp(output);
   disp('Minimum at');
   output=cat(2,minminr(1)-nlag-1,minminc(1)-nlag-1);
   disp(output);
end;

return;

