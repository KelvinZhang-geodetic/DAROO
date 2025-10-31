function [fraction,TP,FP,FN]=iftrue(istepDetected,istep)
%true positive: in istep, detected
%true false: out of istep, detected
%false positive: in istep, not detected 
idec=istepDetected(:);istep=istep(:);
nds=length(idec);
if nargin<2
    TP=nds; TF=0;FP=0;
else
ntp=0;ntf=0;nfp=0;
TP=intersect(idec,istep);
nstep=length(istep);
TP1=[];
for i=1:nstep
    i1=find((idec-istep(i))<=3);
    if ~isempty(i1)
        TP1=[TP1;istep(i)];
    end
end
FP=setdiff(istep,TP);
FN=setdiff(idec,[TP(:); FP(:)]);
fraction=zeros(1,3);
% fraction(1)=length(TP)/length(istep);
fraction(1)=length(TP1)/length(istep);
fraction(2)=length(FP)/length(istep);
fraction(3)=length(FN)/length(istepDetected);

end