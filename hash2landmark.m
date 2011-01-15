function L = hash2landmark(H)
% L = hash2landmark(H)
%  Convert a set of <time hash> pairs ready from store 
%  into a set of 4-entry landmarks < f1 f2 dt t1>.

% Hash value is 20 bits: 8 bits of F1, 6 bits of F2-F1, 6 bits of delta-T
H1 = H(:,1);
H2 = double(H(:,2));
F1 = floor(H2/(2^12));
H2 = H2 - (2^12)*F1;
F1 = F1 + 1;
DF = floor(H2/(2^6));
H2 = H2 - (2^6)*DF;
if DF > 2^5
  DF = DF-2^6;
end
F2 = F1+DF;
DT = H2;

L = [F1,F2,DT,H1];


