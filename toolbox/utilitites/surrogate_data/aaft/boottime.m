function[region,jackstats]=boottime(X,nlag,R,alphax,figs)
% function[region,jackstats]=boottime(X,nlag,R,alpha,figs)
% Time domain statistic for non-linearity
% Bootstraps using AAFT methods then compares observed third order moment with expected
% For more details see: Barnett & Wolff, “A Time-Domain Test for Some Types of Nonlinearity,” IEEE transactions on signal processing, Vol 53, No 1, January 2005
%% Inputs
% X - series as a column vector
% nlag - number of third order moment lags to compute
% R - number of bootstrap replications
% alpha - level of test
% figs - plot figure? 1=Yes, 0=No
%% Outputs
% jackstats - test of null hypothesis using double bootstrap
%             [outside,standardised,upperlimit,pvalue,H0(accept/reject)]
% (outside - total area outside of null hypothesis limits)
%% Other macros called
% third - 3rd order moment
% aaft - AAFT algorithm

% Initial variables;
[n width]=size(X);

% Difference the series;
Xmean=mean(X);
Xdiff=(X-Xmean);

% Get the series 3rd order moment for later use;
[Xtemp]=third(Xdiff,nlag,0);
reglags=nlag+1:nlag+nlag+1;
Xthird=Xtemp(reglags,reglags);
Xthird(1,1)=0; % Remove skewness;

% Get R*3 surrogates using the AAFT method
% First R for initial limits 2nd & 3rd R for bootstrap limits
[aaftsers]=aaft(Xdiff,R*3);

% Run each series through the third order moment;
aaftthird=zeros(nlag+1,nlag+1,R*3);
for k=1:R*3;
   [sersthird]=third(aaftsers(:,k),nlag,0);
   aaftthird(:,:,k)=sersthird(reglags,reglags);
end;
aaftthird(1,1,:)=0; % Remove skewness;

% Get the (1-alpha)th centile at each coordinate and difference from the series;
clevel_l=100*(alphax/2);
clevel_u=100*(1-(alphax/2));
mcent_l1=zeros(nlag+1);
mcent_u1=zeros(nlag+1);
mcent_l2=zeros(nlag+1);
mcent_u2=zeros(nlag+1);
for r=0:(nlag);
   for s=r:(nlag);
      if ((r+s)>0); % Exclude mu(0,0);
         pts1=reshape(aaftthird(r+1,s+1,1:R),R,1); % First limits
         pts2=reshape(aaftthird(r+1,s+1,R+1:(2*R)),R,1); % Second limits
         mcent_l1(r+1,s+1)=prctile(pts1,clevel_l);
         mcent_u1(r+1,s+1)=prctile(pts1,clevel_u);
         mcent_l2(r+1,s+1)=prctile(pts2,clevel_l);
         mcent_u2(r+1,s+1)=prctile(pts2,clevel_u);
      end;
   end;
end;
diff_l=triu(Xthird-mcent_l1); % Just get for s<r;
diff_u=triu(Xthird-mcent_u1); % Just get for s<r;

% Show points significantly higher or lower than limits;
region_u=max(cat(3,diff_u,zeros(nlag+1)),[],3);
region_l=min(cat(3,diff_l,zeros(nlag+1)),[],3);
region=region_u+region_l;
% Total area exceeding limits;
outside=sum(sum(abs(region)));
total=((nlag+1)*(nlag+2)/2)-1; % Total number of points

% Double bootstrap statistic using 2nd set of limits on first set of data;
% 3rd series - limits from 2nd;
jackstat=zeros(R,1);
for jack=(2*R)+1:R*3;
%% Percentile statistic;
   diffjack_u=triu(aaftthird(:,:,jack)-mcent_u2);
   diffjack_l=triu(aaftthird(:,:,jack)-mcent_l2);
   jregion_u=max(cat(3,diffjack_u,zeros(nlag+1)),[],3);
   jregion_l=min(cat(3,diffjack_l,zeros(nlag+1)),[],3);
   jregion=jregion_u+jregion_l;
% Total area exceeding limits;
   jackstat(jack-(2*R),1)=sum(sum(abs(jregion)));
end;

jackstd=std(jackstat);
stan=outside/jackstd;
lowerjack=prctile(jackstat,5);
upperjack=prctile(jackstat,95);
medianjack=median(jackstat);
pjack=length(find(outside<jackstat))/R;
output=cat(2,outside,stan,medianjack,upperjack,pjack);
disp('Double bootstrap statistic');
disp('   observed     obs/std     median-null    95%-null    p-value');
disp(output),
format short g;
testjack=0;
if outside>upperjack
   testjack=1;
end;
jackstats=[outside,stan,upperjack,pjack,testjack];

%% Output maximum difference;
% Upper limit
noplot=0;
if max(max(diff_u))>0;
   [r,c]=find(diff_u==max(max(diff_u)));
   disp('Maximum upper difference at lag');
   disp([r-1,c-1]);
   disp('Size of maximum difference');
   maxdiff=diff_u(r(1),c(1));
   disp(maxdiff);
else;
   disp('Maximum upper difference is zero');
   noplot=noplot+1;
end;

% Lower limit
if min(min(diff_l))<0;
   [r,c]=find(diff_l==min(min(diff_l)));
   disp('Maximum lower difference at lag');
   disp([r-1,c-1]);
   disp('Size of maximum difference');
   maxdiff=diff_l(r(1),c(1));
   disp(maxdiff);
else;
   disp('Maximum lower difference is zero');
   noplot=noplot+1;
end;

% Plot (only points that exceed limits)
if ((figs==1)&(noplot<2));
   blot=repmat(NaN,nlag,nlag);
   blotter=zeros(nlag+1);
   blotter(2:nlag+1,1:nlag)=tril(blot);
   regtoplot=region+blotter;
   surf([0:nlag],[0:nlag],regtoplot);
   shading interp;
   colormap('default');
   set(gcf,'Name','Points of 3rd order moment that exceed limits');
   set(gca,'FontSize',24);
   set(gca,'FontName','Times');
   if ((r-1>0)&(r-1<nlag))
      set(gca,'ytick',[0 r-1 nlag]);
   else;
      set(gca,'ytick',[0 nlag]);
   end;
   if ((c-1>0)&(c-1<nlag))
      set(gca,'xtick',[0 c-1 nlag]);
   else;
      set(gca,'xtick',[0 nlag]);
   end;
   xlabel('r');
   ylabel('s');
   zlabel('UL(r,s)')
elseif ((figs==1)&(noplot==2));
   disp('Plot asked for but no points exceed limits');
end;

return;
