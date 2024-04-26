function stim = deltaStimExtend(bt, nT, v)
% Returns a sparse vector with events at binned timings

if bt<0
    bt=-bt;
    invert=1;
else
    invert=0;
end

bidx = bt <= nT;
bt = bt(bidx); % ignore the events after nT bins

o = ones(numel(bt), 1);

if nargin < 3
    v = o;
else
    v = v(bidx);
end

assert(numel(o) == numel(v));
if invert
    stim = -sparse(bt, o, v, nT, 1);
else
    stim = sparse(bt, o, v, nT, 1);
end