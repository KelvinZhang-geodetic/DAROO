% Generating time series with N<=10 offsets
% Adding white noise with  rand walk + flicker + spikes
% run daroo to detect offsets and remove outliers£»
% output and plot
% 1.construct time (2015-01-01 - 2022-06-30)
clear

res=[];
Nsample=1000;
for ik=1:500
% ik=1;
t0 = (2015:1/365.25:2022.5);
t0=t0(:);
nt=length(t0); 
dt=1/365.25;

% 3. adding trend, annual, colored noise and spikes
v       = randiDistinct([-20 20],1,-3:3);  % mm/yr
res(ik).v=v;
trend   = v * (t0 - t0(1)) ;                            
a = 0 + (5-0).*rand(1,1);b = 0 + (2-0).*rand(1,1);
annual  = a * sin(2*pi*t0) + b * cos(2*pi*t0);                   
noise   = 0.8 * randn(nt,1);

% simple colored noise. first order 
for i = 2:nt, noise(i) = 0.25*noise(i-1) + 0.8*randn; end


%adding flicker
rw_level = 1;
rw = rw_level*cumsum(randn(size(t0)))*sqrt(dt);
fl_level = 0.05;%if flicker noise is too large , it may cover the jump
k = -1;                          % 
f = (1:nt).^k;
f = f(:) /std(f)* fl_level;          
rng(1); wn = randn(size(t0));
fl = real(ifft(fft(wn).*fft(f)));

% it=(1:nt)';without gaps
it=[1:1000 1400:1900 1950:nt]';%gaps
nit=length(it);
t0=t0(it);
% 2.  n  offsets
nsteps=10;
    disp(ik)
rng('shuffle');   %
isteps=randi([10 nit-100],nsteps,1);%e.g. [400 585 1656 2003  3029];
minGap=10;
istep = isteps([true; diff(isteps) > minGap]);
istep=unique(istep);
nstep=length(istep);
offsetSize  = randiDistinct([-20,20],nstep,[-4:4]); 
offsetEpoch = t0(istep);%[2017.15, 2019.45, 2021.70];  
              % mm
dstep       = zeros(nit,1);
for j = 1:numel(offsetEpoch)
    ist = find(t0 >= offsetEpoch(j),1,'first');
    dstep(ist:end) = dstep(ist:end) + offsetSize(j);
end
% outlier/spikes
% spikeIdx = [120, 501, 889, 1400, 1900, 2300];
% spikeAmp = [15, -18, 12, 20, -14, 16];
spikeIdx=randiDistinct([10 nt-100],nsteps,isteps);%[400 585 1656 2003  3029];
spikeAmp  = randi([-20,20],nsteps,1);%[8, -12, 6 18 -30];    

spike    = zeros(nt,1);
spike(spikeIdx) = spikeAmp;

fprintf('Input offset epochs: %s\n', mat2str(istep'));
% d=d(it);
trend=trend(it);annual= annual(it);
noise=noise(it);
spike=spike(it);rw=rw(it);fl= fl(it);  
% simulated raw
d        =  trend + annual + dstep + noise + spike +rw + fl;  
derr0    = 3.0 * ones(nit,1);                         

% 4. run daroo detect steps + remove outlier
% t, d, derr, k, teq, period, timeException
k         = 3;            
period    = 1;%[365.25];    
teq       = [];          
timeExcep = [];           
ismore=0;
res0 = daroo(t0, d, derr0, k, teq, period, timeExcep,ismore);
res(ik).data=[t0 d derr0];
res(ik).istep=istep;
res(ik).istepDetected=res0.istep;
res(ik).Step=offsetSize;
res(ik).StepPred=res0.predStep;
res(ik).iu=res0.iu;
res(ik).io=res0.io;
res(ik).tu=res0.tu;
res(ik).du=res0.du;
res(ik).ab=[a b];
res(ik).dresi=res0.dresi;
res(ik).trend=res0.trend;
res(ik).dperiod=res0.dperiod;
res(ik).dstep=res0.dstep;
res(ik).offsetSize=offsetSize;
[frac,TP,FP,FN]=iftrue(res(ik).istepDetected,res(ik).istep);
fraction(ik,1:3)=frac;

% 5. result printing
fprintf('Detected offset epochs: %s\n', mat2str(res0.istep'));

figure
set(gcf,'units','centimeters',...
        'PaperType','a4',...%papersize [21 29.7]
        'PaperPositionMode','manual',...
        'PaperOrientation','portrait',...
        'PaperUnits','centimeters',...
        'PaperPosition',[1,1,12,15]);
subplot('position',[0.12 0.76 0.85 0.19])
plot(t0,trend + annual,'o', 'Color', [0.5 0.5 0.51],'markersize',2);
ylabel('Cleaned (mm)','fontsize',10,'fontname','Arial')
axis tight
subplot('position',[0.12 0.53 0.85 0.19])
plot(t0,  trend + annual + noise, 'o','markersize',2, 'Color', [0.5 0.5 0.51]);
hold on
box on;axis tight
ylabel('Whitened (mm)','fontsize',10,'fontname','Arial')
subplot('position',[0.12 0.30 0.85 0.19])
plot(t0,trend + annual + dstep + noise + spike, 'Color', [0.5 0.5 0.51])
box on;axis tight
ylabel('Steps+Spikes (mm)','fontsize',10,'fontname','Arial')
subplot('position',[0.12 0.07 0.85 0.19])
plot(t0,d, 'Color', [0.5 0.5 0.51])
xlabel('Time (year)','fontsize',10,'fontname','Arial')
ylabel('Color-Noiseed (mm)','fontsize',10,'fontname','Arial')
box on;axis tight
if ik==1 || mod(ik,20)==0
file1=['daroo_example', num2str(ik),'_signal'];
print('-depsc','-r2000',file1);
print('-djpeg','-r600',file1);

end
close
% 6. ploting comparison
figure
set(gcf,'units','centimeters',...
        'PaperType','a4',...%papersize [21 29.7]
        'PaperPositionMode','manual',...
        'PaperOrientation','portrait',...
        'PaperUnits','centimeters',...
        'PaperPosition',[1,1,12,15]);
subplot('position',[0.12 0.76 0.85 0.19])
plot(t0,d,'o', 'Color', [0.5 0.5 0.51],'markersize',2);
ylabel('Raw (mm)','fontsize',10,'fontname','Arial')
axis tight
subplot('position',[0.12 0.53 0.85 0.19])
plot(t0, d, 'o','markersize',2, 'Color', [0.8 0.8 0.81]);
hold on
plot(res0.tu, res0.trend+res0.dstep+res0.dperiod, '-','LineWidth', 1, 'Color', [0.2 0.2 0.8]);
plot(t0(res0.io), d(res0.io),'o','markersize',3,'markeredgeColor',[0.4 0.4 0.41],'markerfaceColor',[0.2 0.5 0.81]);
plot(t0(res0.istep), d(res0.istep),'v','markersize',6,'markeredgeColor',[0.8 0.2 0.21],'markerfaceColor',[0.8 0.2 0.21]);
p=res0.trend(10)-res0.trend(1);
hl=legend('Raw','Trend','Detected Steps','Detected Outliers');
if p>0
    set(hl,'location','southeast','fontsize',4)
else
    set(hl,'location','southwest','fontsize',4)
end
box on;axis tight
ylabel('Detection (mm)','fontsize',10,'fontname','Arial')
subplot('position',[0.12 0.30 0.85 0.19])
plot(t0,trend + annual + dstep + noise + spike, 'Color', [0.5 0.5 0.51])
hd
plot(t0(res0.iu),res0.dstep+res0.trend, 'Color', [0.8 0.2 0.21]);
box on;axis tight
ylabel('Trend+step (mm)','fontsize',10,'fontname','Arial')
subplot('position',[0.12 0.07 0.85 0.19])
plot(t0(res0.iu),res0.dresi, 'Color', [0.5 0.5 0.51]);
axis tight
yy=max(abs(res0.dresi));set(gca,'ylim',yy*[-2 2]);
xlabel('Time (year)','fontsize',10,'fontname','Arial')
ylabel('Residual (mm)','fontsize',10,'fontname','Arial')
box on;
if ik==1 || mod(ik,20)==0
file1=['daroo_example' num2str(ik)];
print('-depsc','-r2000',file1);
print('-djpeg','-r600',file1);
end
close
save(['daroo_test',nowdate,'.mat'])
end

% for i=1:nexp
%     [frac,TP,FP,FN]=iftrue(res(ik).istepDetected,res(ik).istep);
% fraction(ik,1:3)=frac;
