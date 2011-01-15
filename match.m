function [R,L] = match(lmRef,lmQ)
  IX = 1;

  hashRef = landmark2hash(lmRef);
  recordhash(hashRef);
  hashQ = landmark2hash(lmQ);  
  hashQ = unique(hashQ,'rows');  
  %disp(['landmarks ',num2str(size(lmQ,1)),' -> ', num2str(size(hashQ,1)),' hashes']);
  
  Rt = get_hash_hits(hashQ);
  nr = size(Rt,1);
  if nr > 0
    R = zeros(3);   
    % Find the most popular time offset
    [dts,xx] = unique(sort(Rt(:,1)),'first');
    dtcounts = 1+diff([xx',size(Rt,1)]);
    [vv,xx] = max(dtcounts);
    %    [vv,xx] = sort(dtcounts, 'descend');
    R = [vv(1),dts(xx(1)),size(Rt,1)];
    R = [sum(abs(Rt(:,1)-dts(xx(1)))<=1),dts(xx(1)),size(Rt,1)];
  
    % Sort by descending match count
    [vv,xx] = sort(R(:,1),'descend');
    R = R(xx,:);

    % Extract the actual landmarks
    H = Rt( (Rt(:,1)==R(IX,2)),:);
    % Restore the original times
  %  for i = 1:size(H,1)
  %    hix = find(hashQ(:,3)==H(i,2));
  %    hix = hix(1);  % if more than one...
  %    H(i,1) = H(i,1)+hashQ(hix,2);
  %    L(i,:) = hash2landmark(H(i,:));
  %  end

    % Return no more than 100 hits, and only down to 10% the #hits in
    % most popular
    maxrtns = 100;
    if size(R,1) > maxrtns
      R = R(1:maxrtns,:);
    end
    maxhits = R(1,1);
    nuffhits = R(:,1)>(0.1*maxhits);
    %R = R(nuffhits,:);    
    R
  else
    R = [];
    disp('*** NO HITS FOUND ***');
end
function R = get_hash_hits(H)

global HashTable HashTableCounts
nhtcols = size(HashTable,1);

if min(size(H))==1
  H = [zeros(length(H),1),H(:)];
end

TIMESIZE=16384;

Rsize = 100;  % preallocate
R = zeros(Rsize,2);
Rmax = 0;

for i = 1:length(H)
  hash = H(i,2);
  htime = double(H(i,1));
  nentries = min(nhtcols,HashTableCounts(hash+1));
  htcol = double(HashTable(1:nentries,hash+1));  
  times = round(htcol-TIMESIZE);
  while (Rmax+nentries > Rsize)
    R = [R;zeros(Rsize,2)];
    Rsize = size(R,1);
  end
  dtimes = times-htime;
  R(Rmax+[1:nentries],:) = [ dtimes, repmat(double(hash),nentries,1)];
  Rmax = Rmax + nentries;
end
R = R(1:Rmax,:);
