function S = setTimeProp(S,StartDate, DateFormat)
%SETTIMEPROP Set TimeProperies for all timeseries objects contained in S.
if isstruct(S)
    Prop = fieldnames(S);
else
    Prop = properties(S);
end
for np = 1:length(Prop)
    if isa(S.(Prop{np}),'timeseries')
        S.(Prop{np}).TimeInfo.StartDate = StartDate;
        S.(Prop{np}).TimeInfo.Format = DateFormat;
    end
end
end

