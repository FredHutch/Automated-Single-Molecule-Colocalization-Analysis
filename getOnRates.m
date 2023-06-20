function [POccupied, PNaked]=getOnRates(dataFromAnalysis,framerate)

% AndrewonRates computes pulses time post- and ante- binding. It takes as
% input the 'data' output from ASCMA function. SRC is the interrogated
% channel, and SGC is the 'chaperone' channel.
%
% INPUT:        datafromAnalysis- structure obtained from thr ascma
% function ('data' output)
%               framerate- scalar for the frame per second resolution of
%               the movie
%
% OUTPUT:       POccupied- n-by-2 matrix listing the duration in seconds of
% ON pulses of the red channel occupied by green pulses (first column) as
% well as censorship status (second column)
%               PNaked- m-by-2 matric listing the duration in seconds of ON
%               pulses of the red channels not occupied by green pulses as
%               well as censorship status

R=dataFromAnalysis.SRC.binary;
G=dataFromAnalysis.SGC.binary;
[nF, nS]=size(R);

RG=R+G;
tf=max(RG,[],1)==2;
k=find(tf);
nO=numel(k);
rn=randi(nO,1);
r=k(rn);

B=true(nF,3);
WRo=[]; %overlapping Rchannel lapse
WRnaked=[]; %non overlapping lapse

for i=1:nS
    RS=R(:,i);
    GS=G(:,i);
    dR=[0;diff(RS)];
    dG=[0;diff(GS)];

    Rend=dR==-1;
    Gend=dG==-1;
    Gend=Gend & RS;
    k=find(Gend)+1;
    Gendfill=Gend | Rend;
    B(:,2)=Gendfill;
    Bfill=imfill(B,k+nF);
    MB=Bfill(:,2);
    MB(Rend)=false;
    GG=GS | MB;
    %G(:,i)=GG;

    Rrec=imreconstruct(double(GG),RS,4);
    Grec=imreconstruct(RS,double(GG),4);
    dRr=diff([0;Rrec]);
    dGr=diff([0;Grec]);

    Delay=getDelay(dRr,dGr,GG);
    Delay=Delay*framerate;
    OccCens=NaN(size(Delay));
    [wo,co]=getPulse(Rrec,framerate);
    tf=Delay>0;
    OccCens(tf)=Delay(tf);
    wo(tf)=wo(tf)-OccCens(tf);


    Ro=RS & GG;
    Rnaked=RS & ~Rrec;
    %[wo,co]=getPulse(Ro,framerate);
    [wn,cn]=getPulse(Rnaked,framerate);

    WRo=[WRo;[wo co Delay OccCens]];

    WRnaked=[WRnaked;[wn cn]];
    if i==r
        testRS=RS;
        testGS=GS;
        testGG=GG;
        testRo=Ro;
        testRnaked=Rnaked;
        disp([wo co Delay OccCens]);
        disp([wn cn]);
    end

    

end
POccupied=WRo;
POccupied=array2table(POccupied,'VariableNames',{'RedCh01OccupiedPulseLength','Censorship','GreenCh02ArrivalDelay','RedCh02NakedCensoredPulseLength'});
PNaked=WRnaked;
PNaked=array2table(PNaked,'VariableNames',{'RedCh01NakedPulseLength','Censorship'});

if nO>0

Rbb=diff([0;testRS;0]);
kend=find(Rbb==-1);
kfirst=find(Rbb==1);
kk=[kfirst kend];
kkk=[kk fliplr(kk)];
kkkr=kkk'*framerate;
yyy=zeros(size(kkkr));
yyy(size(yyy,1)/2+1:end,:)=0.75;
figure;
patch(kkkr,yyy,[0.8510    0.3255    0.0980]);

Gbb=diff([0;testGS;0]);
kend=find(Gbb==-1);
kfirst=find(Gbb==1);
kk=[kfirst kend];
kkk=[kk fliplr(kk)];
kkkr=kkk'*framerate;
yyy=zeros(size(kkkr));
yyy(1:size(yyy,1),:)=1;
yyy(size(yyy,1)/2+1:end,:)=1.75;
patch(kkkr,yyy,[0.4667    0.6745    0.1882]);

GGbb=diff([0;testGG;0]);
kend=find(GGbb==-1);
kfirst=find(GGbb==1);
kk=[kfirst kend];
kkk=[kk fliplr(kk)];
kkkr=kkk'*framerate;
yyy=zeros(size(kkkr));
yyy(1:size(yyy,1),:)=2;
yyy(size(yyy,1)/2+1:end,:)=2.75;
patch(kkkr,yyy,[0.7020    0.9412    0.5451]);

Robb=diff([0;testRo;0]);
kend=find(Robb==-1);
kfirst=find(Robb==1);
kk=[kfirst kend];
kkk=[kk fliplr(kk)];
kkkr=kkk'*framerate;
yyy=zeros(size(kkkr));
yyy(1:size(yyy,1),:)=3;
yyy(size(yyy,1)/2+1:end,:)=3.75;
patch(kkkr,yyy,[0.5686    0.0471    0.4314]);

Rnbb=diff([0;testRnaked;0]);
kend=find(Rnbb==-1);
kfirst=find(Rnbb==1);
kk=[kfirst kend];
kkk=[kk fliplr(kk)];
kkkr=kkk'*framerate;
yyy=zeros(size(kkkr));
yyy(1:size(yyy,1),:)=4;
yyy(size(yyy,1)/2+1:end,:)=4.75;
patch(kkkr,yyy,[1.0000    0.8627    0.4902]);

ylim([0 5]);xlabel('Time');ylabel('Pulses');
yticks([0.25 1.25 2.25 3.25 4.25]);
yticklabels({'Rbinary';'Gbinary';'Gadjusted';'Roccupied';'Rnaked'});
set(gcf,'Position',[584   756   690   245]);
set(gca,'Color',[0.1 0.1 0.1])
end




end



function [pwidth, pcensor] = getPulse(trace,framerate)

% INPUT:    -trace numframe-by-1 vector containing the ON/OFF binarized
%           timeseries -framerate scalar second per frame
%
% OUTPUT:   -pwidth vector containing durations of pulses
%           -pcensor vector (same length as pwidth) describing whether
%           pulse is censored (1) or not (0)

etrace = [0; trace; 0];
dtrace = diff(etrace);
kup = find(dtrace == 1);
kdown = find(dtrace == -1);

pwidth = (kdown - kup) * framerate;
if ~isempty(pwidth)
    pcensor = zeros(size(pwidth));
    pcensor(1) = -trace(1); pcensor(end) = trace(end);
else
    pcensor = [];
end

end

function Delay=getDelay(dRr,dGr,GG)

kr=find(dRr==1);
Delay=NaN(numel(kr),1);
kg=find(dGr==1);
for j=1:numel(kr)
    if GG(kr(j))==0
        tf=kg>kr(j);
        kgg=kg(tf);
        kggs=sort(kgg-kr(j),'ascend');
        Delay(j)=kggs(1);
    elseif GG(kr(j))==1
        tf=kg<=kr(j);
        kgg=kg(tf);
        kggs=sort(kgg-kr(j),'ascend');
        Delay(j)=kggs(end);
    end
end
end

