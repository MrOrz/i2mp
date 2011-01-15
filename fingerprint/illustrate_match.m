function [DM,SRO,TK,T] = illustrate_match(DQ,SR,FL,IX)
% [DM,SRO,TK,T] = illustrate_match(DQ,SR,FL,IX)
%     Show graphically which landmarks led to a match.
%     DQ @ SR is the waveform of the query.
%     FL is cell array giving the filenames for the elements in the
%     database.
%     Runs the query, then shows specgrams of query (with
%     landmarks) and the IX'th hit (with landmarks), highlighting the
%     matches.  IX defaults to 1.
%     DM returns the waveform of the matching audio, at sampling
%     rate SRO.
% 2008-12-30 Dan Ellis dpwe@ee.columbia.edu

if nargin < 4;  IX = 1; end

% Run the query
[R,Lm] = match_query(DQ,SR,IX);
% Lm returns the matching landmarks in the original track's time frame

% Recalculate the landmarks for the query
dens = 7;
Lq = find_landmarks(DQ,SR,dens);
% Plot them
subplot(211)
show_landmarks(DQ,SR,Lq);

% Recalculate landmarks for the match piece
tbase = 0.032;  % time base of analysis
%matchtrk = R(IX,1);
matchdt = R(IX,2);
matchtrk =1;
[d,SRO] = mp3read(FL{matchtrk});
Ld = find_landmarks(d,SRO,dens);
% Plot them, aligning time to the query
subplot(212)
show_landmarks(d,SRO,Ld,matchdt*tbase + [0 length(DQ)/SR]);
tmp = length(DQ)*SRO/SR;
if tmp < 1
    tmp = 1;
end
tmp2 = round(matchdt*tbase*SRO);
if (tmp2 < 0)
    tmp2 = 0;
end
DM = d(tmp2+[1:tmp]);
[p,name,e] = fileparts(FL{matchtrk});
name(find(name == '_')) = ' ';
title(['Match: ',name,' at ',num2str(matchdt*tbase),' sec']);

% Highlight the matching landmarks
show_landmarks([],SRO,Lm,[],'o-g');
subplot(211)
show_landmarks([],SRO,Lm-repmat([matchdt 0 0 0],size(Lm,1),1),[],'o-g');
title('Query audio')

TK = matchtrk;
T = matchdt;

function show_landmarks(D,SR,L,T,C)
% show_landmarks(D,SR,L[,T,C])
%    Display the landmarks superimposed on a spectrogram.
%    Rows of L are landmark pairs <t1 f1 f2 dt>
%    T is optional 2-element time range selector (empty for all)
%    C is optional graphic mode specifier (defaults to 'o-r')
%    If D is 2-D, it is taken as the spectrogram up to 4 kHz from
%    find_landmarks (and SR is ignored).
% 2008-12-30 Dan Ellis dpwe@ee.columbia.edu

if nargin < 4
  T = [];
end
if nargin < 5
  C = 'o-r';
end

targetSR = 8000;
% We use a 64 ms window (512 point FFT) for good spectral resolution
fft_ms = 64;
nfft = round(targetSR/1000*fft_ms);

% win overlap in fingerprint calculation
win_overlap = 2;
% win overlap for display here
my_win_overlap = 4;

fbase = targetSR/nfft;
tbase = fft_ms/win_overlap/1000;
my_tbase = fft_ms/my_win_overlap/1000;

if (size(D,1)<3) || (size(D,2)<3)
  % we have an actual soundfile

   if length(D) > 0

   %%%%%%%% vvvvvvvvv Copied from find_landmarks

     if size(D,1) > size(D,2)
       D = D';
     end
     if size(D,1) == 2;
       D = mean(D);
     end
     
     % Resample to target sampling rate
     if (SR ~= targetSR)
       srgcd = gcd(SR, targetSR);
       D = resample(D,targetSR/srgcd,SR/srgcd);
     end

     % Take spectral features
     S = abs(specgram(D,nfft,targetSR,nfft,nfft-nfft/my_win_overlap));

   %%%%%%%% ^^^^^^^^ Copied from find_landmarks
   
   end

else
  S = D;
end

if length(D) > 0
   
  [nr,nc] = size(S);
  tt = [1:nc]*my_tbase;
  ff = [0:nr-1]*fbase;

  imagesc(tt,ff,20*log10(S));
  axis xy
  ca = caxis;
  caxis([-80 0]+ca(2));

end
  
hold on

for i = 1:size(L,1);
  lrow = L(i,:);
  t1q = lrow(1);
  f1q = lrow(2);
  f2q = lrow(3);
  dtq = lrow(4);
  t2q = t1q+dtq;
  t1 = t1q*tbase;
  t2 = t2q*tbase;
  f1 = f1q*fbase;
  f2 = f2q*fbase;
  plot([t1 t2],[f1 f2],C);
end

hold off

if length(T) == 2
  a = axis;
  axis([T(1) T(2) a(3) a(4)]);
end
