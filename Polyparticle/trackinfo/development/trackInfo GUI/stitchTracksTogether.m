%% stitch tracks together... output new tr
function tr=stitchTracksTogether(trInput, tThresh, pixThresh)

% tThresh=5;
% pixThresh=5;

[lifetime,birth,death]=calclifetime(trInput);
stitchPrev=zeros(1,length(trInput));
stitchNext=zeros(1,length(trInput));
for n=1:length(trInput)
    indices=find((birth-death(n)<tThresh)&(birth-death(n)>0));
    tally=[];
    best=inf;
    for N=indices
        dt=birth(N)-death(n);
        dr=sqrt((trInput{N}(1,1)-trInput{n}(end,1))^2+(trInput{N}(1,2)-trInput{n}(end,2))^2);
        test=dt/tThresh+dr/pixThresh;
        if test < 1
            if test<best
                best=test;
                tally=N;
            end
        end
    end
    if ~isempty(tally)
        stitchNext(n)=tally;
        stitchPrev(tally)=n;
    end
end
%
%disp('calculating sequences')
tally=[];
sequences=cell(0);
for n=1:length(trInput)
    if ~ismember(n,tally)
        N=stitchPrev(n);
        sequencePrev=[];
        while N
            sequencePrev=[N sequencePrev];
            N=stitchPrev(N);
        end
        N=stitchNext(n);
        sequenceNext=[];
        while N
            sequenceNext=[sequenceNext N];
            N=stitchNext(N);
        end
        sequences=[sequences {[sequencePrev n sequenceNext]}];
        tally=[tally sequences{end}];
    end
end
%
tr=cell(size(sequences));
%disp('stitching tracks')
for n=1:length(sequences);
    tr{n}=[];
    for N=1:length(sequences{n});
        if N==1
            gap=0;
        else
            gap=birth(sequences{n}(N))-death(sequences{n}(N-1))-1;
        end
        tr{n}=[tr{n}; NaN(gap,10); trInput{sequences{n}(N)}];
    end
end
