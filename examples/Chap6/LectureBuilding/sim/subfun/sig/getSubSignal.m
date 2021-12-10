function ts = getSubSignal(ts,ix)
%Returns timeseries with only one dimensional Data (Data(:,ix))
if nargin <= 1
    ix = 1;
end
ts.Data = ts.Data(:,ix);
end

