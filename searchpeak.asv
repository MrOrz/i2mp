% searchpeak()
% alternative to the slow matlab findpeaks()

function [pks, loc] = searchpeak(vec, dist)
  v1 = [vec(2:end); vec(end)];
  v2 = [vec(1); vec(1:end-1)];
  logical = vec>v1&vec>v2;
  pks = vec(logical);
  loc = find(logical);
  
  % sorting
  [pks, idx] = sort(pks, 'descend');
  loc = loc(idx);

  if (nargin > 1)
    % dist
    idx = ones(length(loc),1) > 0;
    for i = 1:length(loc)
      v = abs(loc-loc(i)) > dist;
      v(1:i) = 1;
      idx = idx & v;
    end

    pks = pks(idx);
    loc = loc(idx);
  end
end