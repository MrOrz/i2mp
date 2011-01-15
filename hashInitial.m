function hashInitial()

global HashTable HashTableCounts

%if exist('HashTable','var') == 0
%   HashTable = [];
%end

nhashes = 2^20;

% 1M hashes x 32 bit entries x 100 entries = 400MB in core
%maxnentries = 100;
maxnentries = 20;

disp(['Max entries per hash = ',num2str(maxnentries)]);

%if length(HashTable) == 0
  HashTable = zeros(maxnentries, nhashes, 'uint32');
  HashTableCounts = zeros(1, nhashes);
%end