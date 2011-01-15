function H = landmark2hash(L,S)
% input : [t f1 f2 dt] = [起始時間 頻率一 頻率二 間隔時間]
% output : 
%  Convert a set of 4-entry landmarks <t1 f1 f2 dt> 
%  into a set of <songid time hash> triples ready to store.
%  S is a scalar songid, or one per landmark (defaults to 0)

% Hash value is 20 bits: 8 bits of F1, 6 bits of delta-F, 6 bits of delta-T

if nargin < 2
  S = 0;
end
if length(S) == 1
  S = repmat(S, size(L,1), 1);
end

H = uint32(L(:,1));
% Make sure F1 is 0..255, not 1..256
F1 = rem(round(L(:,2)-1),2^8);
DF = round(L(:,3)-L(:,2));
if DF < 0
  DF = DF + 2^8;
end
DF = rem(DF,2^6);
DT = rem(abs(round(L(:,4))), 2^6);
H = [S,H,uint32(F1*(2^12)+DF*(2^6)+DT)];
