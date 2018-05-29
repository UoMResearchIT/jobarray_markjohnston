                                                % function [nRuns,runToggle,nHalts,haltToggle]=runhalt(position)
% segment track into halts and runs
% requires position along contour as a function of time
%(1) scan through track
%(2) look for stationary sections 'halts', where particle moves less than Lpix in more than Tscale frames
%(3) cut out halts 
%(4) examine remaining sections 'runs', discard if too short in distance or time
function [nRuns,runToggle,nHalts,haltToggle]=runhalt(position)

Lpix=1;
LtooLittle=1;
LtooMuch=inf;
haltTscale=14; %Changed from 5 to 14 90304
%runTscale=5;
tSmootherScale=1;
T=1;
runHaltToggle=false(size(position));
discardToggle=[false; diff(position)>LtooMuch];
positionConv=conv(position,ones(tSmootherScale,1)/tSmootherScale);
positionConv=positionConv(1:end-tSmootherScale+1);

while T<=length(positionConv)
    if ~isnan(positionConv(T))
        indices=cumsum(abs(positionConv-positionConv(T))>Lpix);
        toggle=indices==indices(T);
        indices=find(toggle);
        indices=indices(2:end);
        if length(indices)>=haltTscale
            runHaltToggle(indices)=1;
        end
    end
    %T=find(diff(toggle)==-1)+1;
    T=T+1;
end
%[haltToggle nHalts]=bwlabel(runHaltToggle&~discardToggle);
[runToggle nRuns]=bwlabel(~runHaltToggle&~discardToggle);
% check runs again... reject any that are too short in time or distance
for N=1:nRuns
    indices=find(runToggle==N);
    %if (length(indices)<runTscale)&&(position(indices(end))-position(indices(1))<LtooLittle)
    if (position(indices(end))-position(indices(1))<LtooLittle)
        runHaltToggle(indices)=1;
    end
end
[haltToggle nHalts]=bwlabel(runHaltToggle&~discardToggle);
[runToggle nRuns]=bwlabel(~runHaltToggle&~discardToggle);
