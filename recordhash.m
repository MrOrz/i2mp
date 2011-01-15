function recordhash(hash)
  global HashTable HashTableCounts   
  for i=1:size(hash,1)
    key2 = 1 + hash(i,2);  % Avoid problems with key == 0    
    numberOfEntry =  HashTableCounts(key2) + 1;
    if numberOfEntry <= size(HashTable,1)	  
	  r = numberOfEntry; % put entry in next available slot
    else      
      r = ceil(numberOfEntry * rand(1)); % If full, choose a random size
    end
    % Put it in the hash table
    if r <= size(HashTable,1)   
      hashValue = int32(round(hash(i,1)));      
      HashTable(r,key2) = hashValue;      
    end
    HashTableCounts(key2) = numberOfEntry; % Update the HashTableCounts value
  end