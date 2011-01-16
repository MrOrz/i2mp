function T = match(lmRef,lmQ)
  % Find the time offset of two landmark;
  
  % Step 1: Construct hash table by the reference audio
  hashRef = unique(landmark2hash(lmRef),'rows');
  recordhash(hashRef);
  
  % Step 2: Find the hash of query audio
  hashQ = unique(landmark2hash(lmQ),'rows');  
  
  % Step 3: Check in the hash table
  Rt = checkHash(hashQ);  
  % Step 4: If something match then print out
  if size(Rt,1) > 0
    R = zeros(2);  
    % Unique the offser value
    [dt,dt_x] = unique(sort(Rt(:,1)),'first');    
    dtcounts = 1 + diff([dt_x',size(Rt,1)]);    
    % Sort the counts    
    [dt_count,dt_count_x] = sort(dtcounts, 'descend'); 
    % The top 10 
    for i = 1:10
      % R = [Time difference count ,time difference]
      R = [dt_count(i),dt(dt_count_x(i))];        
      % Sort by descending match count      
      [sort_dt_count,sort_dt_count_x] = sort(R(:,1),'descend');    
      R = R(sort_dt_count_x,:);
      display(['The rank ',int2str(i),' time offset is ',int2str(R(2)),' .']);
      display(['Which has ',int2str(R(1)),' matched landmark.']);
      display(' ')
      if(i==1)
          T=R(2);
      end
    end    
  else    
    T=0;
    disp('No Match');
  end
  
function R = checkHash(hash)
  % Find the set of the offset
  
  global HashTable HashTableCounts
  nhtcols = size(HashTable,1);
  Rsize = 100;  % preallocate
  R = zeros(Rsize,1);
  Rmax = 0;
  for i = 1:length(hash)    
    key2 = 1 + hash(i,2);
    cur_t = double(hash(i,1));
    numberOfEntry = min(nhtcols,HashTableCounts(key2));
    htcol = double(HashTable(1:numberOfEntry,key2));  
    ref_t = round(htcol);
    while (Rmax + numberOfEntry > Rsize)
      R = [R;zeros(Rsize,1)];
      Rsize = size(R,1);
    end    
    R(Rmax+[1:numberOfEntry],:) = ref_t - cur_t;
    Rmax = Rmax + numberOfEntry;
  end
  R = R(1:Rmax,:);
