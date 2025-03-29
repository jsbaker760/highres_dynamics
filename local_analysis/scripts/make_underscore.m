function out=make_underscore(in)

out=arrayfun(@(x) ['_' in(x)],1:numel(in),'UniformOutput',false);
out=horzcat(out{:});