function res = daroo(t,d,derr0,k,teq,period,timeException,ismore)
% Detection And Removal of Offsets and Outliers (DAROO)
% [tu,du,iu,io,istep,dresi] = daroo(t,d,derr0,k,teq,period,timeException)
%Inputs
% required
% t,time£¨yyyy.ddd)
% d,observations£¨mm£©
% optional
% derr0,formal error of observations
% k, factor for	MAD£¬default=3
% teq, given event occurrence times
% period, seasonal periods, 1, 0.5 for annual and semi-annual, respectively£¬ default=[]
% timeException, time windows to be excluded, in format [t11 t12; t21 t22;...]. such as postseismic stage 
% ismore, whether refined jump detection should be performed on residuals, default=0
%Outputs
%res, a struct array with fields as follows
% res.istep, index of Jump epoches
% res.io, index of outlier epoches
% res.iu, index of usable epoches
% res.tu, time of usable epoches
% res.du, data of usable epoches
% res.dstep, time series of all steps
% res.trend, trend fitted by usable epoches
% res.dperiod, period term for usable epoches
% res.dresi, residuals after removing offsets and outliers
%originally on 2025-08-12
%modified on 2025-09-28
%Keliang Zhang
%Institute of Geology, CEA
if nargin<3 || isempty(derr0),    derr = ones(size(t));end
if nargin<4 || isempty(k),           k = 3;           end     
if nargin<5 || isempty(teq),       teq = [];          end  
if nargin<6 || isempty(period), period = 1;           end
if nargin<8 || isempty(ismore), ismore = 1;           end
if nargin<7 || isempty(timeException), timeException = []; end

warning('off');
%==================================================================================
%---------- 0. preprocessing -----------------------------------------------------%
%==================================================================================
t1      = t(:);%              
t00     = t1(:);
d00     = d(:);
nt00    = numel(t00);
its00   = (1:nt00)';
isa     = find(isnan(d00));              
if nargin<3 || isempty(derr0)
    derr00 = ones(nt00,1);
else
    derr00 = derr0(:);
end
%==================================================================================
%Step1.1:Jump detection by differencing displacement and instantaneous velocities % 
%==================================================================================                         
istep=[];istep1=[];
istep1= dualStepIndex(t00, d00);
%NMS  remove nearest
istep = nmsSteps(istep1, t00, d00, 10);                                                                                                                                             
    
teq0=teq;
teqi=[];
if ~isempty(istep)
    teqi=t(istep);
end
teq=unique([teq0;teqi]);
isn   = find(~isnan(d00));
xsn   = tsdetrend(t00(isn),d00(isn),derr00(isn),period,teq);  
%==================================================================================
%Step1.2. <optional> Jump detection from residuals after large offsets            % 
%==================================================================================   
ismore = 0;
if ismore
   [idxAll, chi2Val] = detectAfterLargeJump(t00(isn), xsn);
   idx=isn(idxAll);%global index
   istep2 = nmsSteps(idx, t00, d00, 10);
   istep=unique([istep;istep2]);
   istep = nmsSteps(istep, t00, d00, 10);
   teq2=[];
   if ~isempty(istep)
       teq2=t(istep);% the final result of Jumps
   end
   teq=unique([teq0;teq2]);
   xsn=tsdetrend(t00(isn),d00(isn),derr00(isn),period,teq);  
end
%==================================================================================
%Step2. outlier detection and removal % 
%==================================================================================   
%-------first round outliers removal  ----------
i0      = gross(xsn,k);                            % 3¦Ò+leave-one-out
io0     = isn(i0);                                 % global indices
iu00    = setdiff(its00,[io0;isa]);                % sustained indices
t0      = t00(iu00);
d0      = d00(iu00);
derr0   = derr00(iu00);

%-2. offsets fit after outliers removed  ----------
dres0 = tsdetrend(t0,d0,derr0,period,teq); 

%---------- 3. exception time windows ----------
if nargin<5 || isempty(timeException)
    if isempty(teq0)
        timeException = [];
    else
        % time window for postseismic stage default 3 months 
        timeException = [];
        for i = 1:length(teq0)
            teq3 = teq0(i)+90/365;
            timeException(i,:) = [teq0(i) teq3];
        end
    end
end
ixs = []; dixs = []; tixs = []; derrixs = [];
if ~isempty(timeException)
    tEx = timeException;
    if size(tEx,1)==2 && size(tEx,2)==1, tEx = tEx.'; end
    for i = 1:size(tEx,1)
        i1 = find(t0>=tEx(i,1) & t0<=tEx(i,2));
        ixs = [ixs; i1];
    end
end
its = setdiff((1:length(t0))',ixs);

%---------- 4. second round outlier removal ----------
t  = t0(its);
y  = d0(its);
x  = dres0(its);
xerr = derr0(its);
i1 = gross(x,k); 
io1 = iu00(its(i1));

%---------- 5. iterative until no outliers----------
i2 = [];io2=[];
while ~isempty(i1)
    its   = setdiff(its, i1);
    tu    = t0(its);
    du    = d0(its);
    duerr = derr0(its);
    dres1= tsdetrend(tu,du,duerr,period,teq);
    x1    = dres1;
    i10   = gross(x1,k);
    i1    = its(i10);
    i2    = unique([i2(:); i1(:)]);
end
if ~isempty(i2)
io2 = iu00(i2);
end

%Results ----------
io  = unique([io0; io1; io2]);
iu  = setdiff(its00, [io; isa]);
tu  = t00(iu);
du  = d00(iu);
duerr=derr00(iu);
[dresi,G,m]=tsdetrend(tu,du,duerr,period,teq);
res.raw=[t00 d00 derr00];
res.iu=iu;
res.io=io;
res.tu=tu;
res.du=du;
res.istep=istep;
res.istepNewInd=[];
if ~isempty(istep)
res.istepNewInd=find(ismember(iu,istep));%istep is index in original data,  
end
res.dresi=dresi;
res.trend=G(:,1:2)*m(1:2,:);

res.dperiod=[];
if ~isempty(period)
    np=length(period);
    idxp=2+(1:np*2);
    yp=G(:,idxp)*m(idxp,:);
    res.dperiod=yp;
end
res.predStep=[];
res.dstep=[];
res.eqdstep=[];
res.eqStep=[];
res.noneqdstep=[];
res.noneqStep=[];
nG=size(G,2);
nstep=nG-2-np*2;
if nstep>0
    idxc=(2+2*np+1):nG; %==(2+2*np)+(1:nstep);
    ys=G(:,idxc)*m(idxc,:); 
    res.dstep=ys;%total
    res.predStep=m(idxc,:);
if ~isempty(teq0)
    ieq=find(ismember(teq,teq0));
    if ~isempty(ieq)
        ys=G(:,idxc(ieq))*m(idxc(ieq),:);    
        res.eqdstep=ys;%due to Eq
        res.eqStep=m(idxc(ieq),:);
    if nstep>length(ieq)%due to other type of Jumps
        is=setdiff(1:nstep,ieq);
          ys=G(:,idxc(is))*m(idxc(is),:);    
        res.noneqdstep=ys;
        res.noneqStep=m(idxc(is),:);
    end
    end
end
end

end

function keepIdx = nmsSteps(stepIdx, t, y, win)
if isempty(stepIdx),
    keepIdx=[];return;
else
n = numel(y);
tStep = t(stepIdx);
[~, ord] = sort(tStep);
stepIdx  = stepIdx(ord);
tStep    = tStep(ord);
keep = true(size(stepIdx));

for k = 1:numel(stepIdx)
    if ~keep(k), continue; end
    tooClose = abs(tStep - tStep(k)) <= win/365;
    idxLocal = stepIdx(tooClose);      
  
    varVec = zeros(numel(idxLocal), 1);
    for j = 1:numel(idxLocal)
        lj = idxLocal(j);
        
        y0 = mean(y(1:lj));
        y1 = mean(y(lj+1:end));
        res = [y(1:lj) - y0; y(lj+1:end) - y1];
        varVec(j) = var(res);
    end
    [~, best] = min(varVec);           
    keepIdxLocal = find(tooClose);
    keep(keepIdxLocal([1:(best-1) (best+1):end])) = false;
end
keepIdx = stepIdx(keep);
end
end
function [idxStep,res] = dualStepIndex(t, y,k,winDays)
if nargin<3|| isempty(k),       k = 3;        end
if nargin<4 || isempty(winDays),winDays = 10; end

t=t(:);y=y(:);nt=length(t);
idxStep = [];

[f, P] = lombscargle(y, t, 'normalized',  0.5, 4);
idx = f > 0.05 & f < 0.5;               % 0.05¨C0.5 cpy
slope = polyfit(log10(f(idx)), log10(P(idx)), 1);
flickerFlag = abs(slope(1) + 1) < 0.3; 

dy  = diff(y);%dy=[dy;dy(nt-1)];
DmadBase = mad(dy)/0.6745;
thres = DmadBase * (2.0 + 0.5 * ~flickerFlag); 

dt=diff(t);igap=[];
if min(dt)<1/360
    dt1=dt*365;
    igap=find(dt1>10);
   vRate=diff(y)./dt1;
   dtm=min(dt1);
else
   vRate=diff(y)./dt;
   igap=find(dt>10);
   dtm=min(dt);
end
if ~isempty(igap)
   for i=1:length(igap)
       t1=t(1:igap(i));y1=y(1:igap(i));
       A=[ones(igap(i),1) t1 period2x(t1,1)];
       x=A\y1;
       t2=[t(1:igap(i));t(igap(i)+1)-min(dt)];%expolation
       A2=[ones(igap(i)+1,1) t2 period2x(t2,1)];
       y1pred=A2*x;
       D=y(igap(i)+1)-y1pred(igap(i)+1);
       vRate(igap(i))=D/dtm;
       dy(igap(i))=D;
   end
end

loc  = [find(dy > thres);find(dy < -thres)];
locD = loc + 1;               
locD=setdiff(locD,[1:3 nt-2:nt]);
[idxD,dD, Dratio, DpVal] = jumpModelFTest(t, y, locD, winDays);
res.dD=dD;res.Dratio=Dratio;res.DpVal=DpVal;

VmadBase = mad(vRate)/0.6745;
thresV  = VmadBase * (2.0 + 0.5 * ~flickerFlag);
locv  = [find(vRate > thresV); find(vRate < -thresV)];
locV=locv+1;
locV=setdiff(locV,[1:3 nt-2:nt]);
[idxV,dV, Vratio, VpVal] = jumpModelFTest(t, y, locV, winDays);
res.dV=dV;res.Vratio=Vratio;res.VpVal=VpVal;
idxStep=union(idxD,idxV);
end

function [iStep,d, varRatio, pVal] = jumpModelFTest(t, y, idxStep, winDays, dThres,alpha)
% good = (abs(d) > dThres) & (varRatio < 0.8) & (pVal < alpha);
if nargin<4 || isempty(winDays),winDays = 30; end
if nargin<5 || isempty(dThres),dThres = 2; end
if nargin<6 || isempty(alpha),alpha = 0.05;end
if isempty(idxStep)
    iStep=[];d=[];varRatio=[]; pVal=[];
else
id1=diff([idxStep(1);idxStep]);
in1=find(id1~=1);
idxStep=idxStep(in1);
id3=diff([idxStep(1)-3;idxStep]);
in3=find(id3>3);
idxStep=idxStep(in3);
t=t(:);y=y(:);
dt0=diff(t);
if median(dt0)<1/300,dt0=dt0*365;end
dt  = min(dt0);
w   = winDays/round(dt);
n   = numel(t);
X=[];
if winDays>1000,X=period2x(t,1);end
isStep=zeros(numel(idxStep),1);
d=isStep;
varRatio=[]; pVal=[];
nidx=numel(idxStep);
%automatically change the time window length
% to avoid the influence of jumps on nearby jumps
for k = 1:nidx
%     disp(k)
    idxC=idxStep(k);
    i1  = max([1, idxC-w]);
    i2  = min([n idxC+w]);
    if nidx>1
    if k==1
        if idxStep(2)>idxStep(1)+3
           i2  = min([n idxC+w idxStep(2)-1]);
        end
    elseif k==nidx
        if idxStep(k)>idxStep(k-1)+3
           i1  = max([1, idxC-w idxStep(k-1)+1]);
        end
    else
        if idxStep(k)>idxStep(k-1)+3
          i1  = max([idxC-w idxStep(k-1)+1]);
        end
        if idxStep(k+1)>idxStep(k)+3
           i2  = min([idxC+w idxStep(k+1)-1]);
        end
    end
    end
    iks=i1:i2;
tWin = t(iks);
yWin = y(iks);
m    = numel(tWin);
% 1. H0: there is no jump in the time window: linear fit
X0 = [ones(m,1), tWin];
if m>1000
    X0 = [X0 X(iks,:)];
end
b0 = X0\yWin;
res0 = yWin - X0*b0;
RSS0 = sum(res0.^2);    
df0  = m - 2;            % degree of freedom
% 2. H1: there is a jump in the time window: segment fit + jump
ni1=length(i1);ni2=length(i2);
y2=y(i2);y1=y(i1);t1=t(i1);t2=t(i2);
meany2=nanmean(y2);meany1=nanmean(y1);
if t1(ni1)-t1(1)<5
    if t2(ni2)-t2(1)<5
        d(k,1)=meany2-meany1;
        res1=[y1-meany1;y(i2)-meany2];
    else
        X2=[ones(ni2,1) t2];
        p=polyfit(t2,y2);
        y2pred=polyval(p,t2);
        d(k,1)=y2pred(1)-meany1;
        res1=[y1-meany1;y(i2)-y2pred];
    end
else
    if ni2<5
        X1=[ones(ni1,1) t1];
        p=polyfit(t1,y1);
        y1pred=polyval(p,t1);
        y1predR=polyval(p,t2(1)-min(dt0));%one day before the break
        d(k,1)=meany2-y1predR;
        res1=[y1-y1pred;y(i2)-meany2];
    else
      dum  = double(tWin >= tWin(idxC-i1+1));
      X1   = [ones(m,1), tWin, dum];
       if m>1000
          X1 = [X1 X(iks,:)];
       end
       b1 = X1\yWin;
       res1 = yWin - X1*b1;
       d(k,1) = b1(3);                 % jump 
    end
end
df1  = m - 3;    
RSS1 = sum(res1.^2);

% 3. jump & variation ratio
varRatio(k,1) = RSS1/RSS0;

F(k,1) = (RSS0-RSS1)/(RSS0/df0);
pVal(k,1) = 1-fcdf(F(k,1),1,df0); 

isStep(k,1) = (abs(d(k,1))>dThres) & (varRatio(k,1)<0.6) & (pVal(k,1)<alpha);
end
iStep1=idxStep(find(isStep));
 [~, ord] = sort(iStep1, 'ascend');
iStep1 = iStep1(ord);
keep = true(size(iStep1));
for k = 2:numel(iStep1)
    if abs(iStep1(k) - iStep1(k-1)) < winDays, keep(k) = false; end
end
iStep = iStep1(keep);
end
end

function [dres,G,m]=tsdetrend(t,d,derr,period,teq)
% [dres,v,vr]=tsdetrend(t,d,derr,period,teq)
%inputs
%      t: observation time,nt in format of 2020.0014
%      d: observations of dimension nd, [nt,nd]=size(d);
%   derr: observation formal error of same dimension as d
% period:
%    teq:
%outputs
%   dres: residuals after removing the linear trend from observations d
%      v: linear rates
%     vr: confidence interval of linear rates
warning('off');
[nt,nd]=size(d);
if nargin<3 || isempty(derr),derr=ones(nt,nd);end
if nargin<4 || isempty(period),period=[];np=0;end
if nargin<5 || isempty(teq),teq=[];end
neq=length(teq(:));
dt=t(:)-t(1);
G0=[dt ones(nt,1)];

for kk=1:nd
    i1=find(isnan(d(:,kk)));
    if length(i1)>=1
        nti=nt-length(i1);
        iu=setdiff(1:nt,i1);
    else
        iu=1:nt;
        nti=nt;
    end
    
dt=t(iu)-t(iu(1));
G0=[dt ones(nti,1)];
Gp=[];
if ~isempty(period)
    np=length(period);
    Gp=period2x(t(iu),period);
end
eps=10^-6;
Gc=[];
if ~isempty(teq)
    for i=1:neq 
        teqi=teq(i);
        i0=find(t<teqi);
        i1=find(t>=teqi);
        ni0=length(i0);bzero=zeros(ni0,1);
        ni1=length(i1);aones=ones(ni1,1);
        Gc=[Gc [bzero;aones]];
    end
end
G=[G0 Gp Gc];
if any(isnan(derr(iu,kk))),r=ones(nti,1);else r=derr(iu,kk);end
 [m(:,kk),dm(:,kk)]=lscov(G,d(iu,kk),1./r.^2);
end
dpred=G*m;
v=m';
vr=dm';
dres=d-dpred;

end

function X=period2x(t,period)
  dtn=min(diff(t));
period=augperiod(period,dtn);
  X=[];
  for i=1:length(period)
      X=[X sin(2*pi/period(i)*t) cos(2*pi/period(i)*t)];
  end
end
function period=augperiod(period,dtn)
if period
    if dtn==1,
       if period<1
           period=365*period;
       end
    elseif dtn<1/360
        if period>1
           period=period/365;
       end
    end
end
end
function [pxx, f] = lombscargle(x, t, normalized, fmax, ofac)
%   [pxx, f] = lombscargle(x, t);
%   [pxx, f] = lombscargle(x, t, true, 0.5, 4);
%inputs:
%   x         : signal vector
%   t         : time
%   normalized: true 
%   fmax      : max freq (cycles/unit t)£¬[]¡úNyquist
%   ofac      : oversampling factor£¬default 4
%outputs:
%   pxx       : normalized PSD
%   f         : freqency (cycles/unit t)
%
if nargin < 3 || isempty(normalized), normalized = true; end
if nargin < 4 || isempty(fmax),      fmax = [];          end
if nargin < 5 || isempty(ofac),      ofac = 4;           end

ok   = ~isnan(x);
t    = t(ok);  x = x(ok);
n    = numel(t);
if n < 3,  pxx = []; f = []; return; end


Tspan = t(end) - t(1);
fMin  = 1 / (ofac * Tspan);
if isempty(fmax)
    fMax = 1 / (2 * min(diff(sort(t))));
else
    fMax = fmax;
end
f    = fMin : fMin : fMax;
nf   = numel(f);
f=f(:);

x0   = x - mean(x);
varx = mean(x0.^2);


pxx  = zeros(nf,1);
for k = 1:nf
    w = 2*pi*f(k);
    tan2wTau = sum(sin(2*w*t)) / sum(cos(2*w*t));
    tau      = atan(tan2wTau) / (2*w);      
    arg      = w*(t - tau);
    XC       = sum(x0 .* cos(arg));
    XS       = sum(x0 .* sin(arg));
    P        = 0.5 * (XC^2 / sum(cos(arg).^2) + XS^2 / sum(sin(arg).^2));
    pxx(k)   = P / varx;                    
end

end

function [idxAll, chi2Val] = detectAfterLargeJump(t, x, alpha, W)
% Inputs£º
% t(N¡Á1) 
% x(N¡Á1) residuals
% alpha confidence level
% W 
% Outputs£º
%  idxAll(K¡Á1) 
%  chi2Val(K¡Á1) ¦Ö2 value
if nargin < 3, alpha = 0.01; end
if nargin < 4, W = 1000; end
N = numel(x);  
t = t(:);  x = x(:);  
dt = median(diff(t));

sigma = robust_sigma(x);                          
chi0=chi2inv(1 - alpha, 1);

chi2Val = NaN(N,1);  idxAll=[];tauSig = [];
for k = 2+1:N-2
    idx = k;                                        
    winLoc = max(1,idx-W):min(N,idx+W);            
    tLoc   = t(winLoc);  xLoc = x(winLoc);
     Nk=length(winLoc);
     k1=winLoc(1):k;
     nk1=length(k1);

    locMAD = median(abs(xLoc - median(xLoc)));
    if locMAD == 0, continue; end                   


    mu0 = mean(xLoc(1:find(tLoc==t(idx),1)));       
    mu1 = mean(xLoc(find(tLoc==t(idx),1)+1:end));   
    A     = abs(mu0 - mu1);                         
    chi2  = sqrt(nk1*(Nk-nk1)/Nk) *A/locMAD;  %Williams 2003        


    if chi2 > chi0
        chi2Val(idx) = chi2;
        idxAll(end+1,1)=idx;
    end
end
idxAll = unique(sort(idxAll));                      
end

function [io,iu,dnew]=gross(dresi,k,indException)
%dresi, one column = one variable
%k, factor
%indException: exception(s) not applied in the gross error check, for example the early postseismic period
%              ind11 ind12
%              ind21 ind22
%               ...   ...
%     i1=find(t>=20210522 & t<=20210630);
%     i2=find(t>=20220108 & t<=20220131);
%     ind=[i1(1) i1(end);i2(1) i2(end)];
%     [io1,iu,dnew]=gross(dresi(:,2),K,ind);
if nargin<2 || isempty(k),k=3;end
nx0=length(dresi);
iu0=(1:nx0)';iu1=iu0;
if nargin<3,indException=[];end
if ~isempty(indException)%except for time span during early postseismic stage
    ind=indException;
    inds=[];
    for i=1:size(ind,1)
         ind1=find(iu0>=ind(1) & iu0<=ind(2));
         inds=[inds;ind1];
    end
    iu1=setdiff(iu0,inds);
end
x1=dresi(iu1);dnew=dresi;
isn=find(~isnan(x1));
x=x1(isn);
nx=length(x);ix=1:nx;
dx=abs(x-nanmean(x));
sd=std(x);
i01=find(dx>=k*sd);
i02=io_leaveoneout(x,k);%Index Outliers with leave-one-out cross-validation
i0=unique([i01;i02]);
io1=i0;io=io1;iu=[];
iu1=setdiff(ix,io1);
%if there are any outliers, then find outliers in the remaining data,until there is no outliers
while ~isempty(i0)  
    iu1=setdiff(ix,io1);%remove io from whole
    x1=dx(iu1);
    dx1=abs(x1-mean(x1));
    sd1=std(x1);
    i1=find(dx1>k*sd1);
    i0=iu1(i1); 
    io1=unique([io1;i0(:)]);
end
iu=1:nx0;
if ~isempty(io1)
    io=iu0(isn(io1));
    iu=setdiff(1:nx0,io);
    dnew(io)=mean(dresi(iu));
end

end

function [io,iu] = io_leaveoneout(residuals,k)
if nargin<2 || isempty(k), k = 3; end
x=residuals(:);
nx  = length(x);
ix  = 1:nx;

med = median(x);
dx  = iqr(x);
io1 = find(x > med+k*dx | x < med-k*dx);
% ---------- 2.  MAD ----------
[~, ord]   = sort(x);
xSort = x(ord);
n     = numel(xSort);

io2 = [];
for i = 1:n
    if i == 1
        xRem = xSort(2:end);
    elseif i == n
        xRem = xSort(1:n-1);
    else
        xRem = xSort([1:i-1 i+1:end]);
    end
    mr  = median(xRem);                
    mad = median(abs(xRem - mr)) / 0.6745; 
    z   = abs(xSort(i) - mr) / (mad + eps);
    if z >= k
        io2 = [io2; ord(i)]; %
    end
end

io = unique([io1; io2]);
if ~isempty(io)
    iu = setdiff(ix, io);
else
    iu = ix;
end
end

function s = robust_sigma(x)
s = mad(diff(x))/sqrt(2);  
end
