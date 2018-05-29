% CALCLIFETIME
% ---
% usage lifetime=calclifetime(tr)
% or [lifetime,birth,death]=calclifetime(tr)
function varargout=calclifetime(tr)

birth=[];
death=[];
lifetime=[];

for i=1:length(tr)
    if ~isempty(tr{i})
        birth(i)=tr{i}(1,3);
        death(i)=tr{i}(end,3);
        lifetime(i)=death(i)-birth(i);
    else
        birth(i)=0;
        death(i)=0;
        lifetime(i)=0;
    end
end

if nargout>=1
    varargout{1}=lifetime;
end
if nargout>=3
    varargout{2}=birth;
    varargout{3}=death;
end
