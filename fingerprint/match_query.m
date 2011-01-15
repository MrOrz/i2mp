function [R,L] = match_query(D,SR,IX)
% [R,L] = match_query(D,SR,IX)
%     Match landmarks from an audio query against the database.
%     Rows of R are potential maxes, in format
%      songID  modalDTcount modalDT
%     i.e. there were <modalDTcount> occurrences of hashes 
%     that occurred in the query and reference with a difference of 
%     <modalDT> frames.
%     L returns the actual landmarks that this implies for IX'th return.
% 2008-12-29 Dan Ellis dpwe@ee.columbia.edu

if nargin < 3
  IX = 1;
end


% Target query landmark density
% (reference is 7 lm/s, query can be maybe 4x denser?)
dens = 20;

% collapse stereo
if size(D,2) == 2
  D = mean(D,2);
end

%Rt = get_hash_hits(landmark2hash(find_landmarks(D,SR)));
Lq = find_landmarks(D,SR, dens);

%%Lq = fuzzify_landmarks(Lq);
% Augment with landmarks calculated half-a-window advanced too
landmarks_hopt = 0.032;
offset = 32;
SR
for i = 1:offset-1
  Lq = [Lq;find_landmarks(D(round(landmarks_hopt*i/offset*SR):end),SR, dens)];
end
% add in quarter-hop offsets too for even better recall

Hq = unique(landmark2hash(Lq), 'rows');

disp(['landmarks ',num2str(size(Lq,1)),' -> ', num2str(size(Hq,1)),' hashes']);
Rt = get_hash_hits(Hq);

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
  for i = 1:size(H,1)
    hix = find(Hq(:,3)==H(i,2));
    hix = hix(1);  % if more than one...
    H(i,1) = H(i,1)+Hq(hix,2);
    L(i,:) = hash2landmark(H(i,:));
  end

  % Return no more than 100 hits, and only down to 10% the #hits in
  % most popular
  maxrtns = 100;
  if size(R,1) > maxrtns
    R = R(1:maxrtns,:);
  end
  maxhits = R(1,1);
  nuffhits = R(:,1)>(0.1*maxhits);
  %R = R(nuffhits,:);

else
  R = [];
  disp('*** NO HITS FOUND ***');
end
