function hashInitial()
  global HashTable HashTableCounts  
  
  maxnentries = 100;  
  nhashes = 2^20;  % 20bits hash
    
  HashTable = zeros(maxnentries, nhashes, 'uint32');
  HashTableCounts = zeros(1, nhashes);