%% Open file and read human genome
clc; clear;
f1='humhbb.txt'; 
f=fopen(f1);        
seq=fread(f,'*char')';

if (ispc | ismac)
    eol = 13;   % pc's and mac's use 13 (carriage return) for end of line
else
    eol = 10;   % others use 10 (line feed)
end

% Find beginning of DNA sequence and read to end-of-line (cr='13')
x=findstr(seq,'ORIGIN');   seq=seq(x:end);
x=findstr(seq,eol);         seq=seq(x:end);

seq=seq(isletter(seq));   % just take letters (drop numbers and spaces).
if strcmp(f1,'humhbb.txt')
    seq = seq(62205:63628);
end
s1a=double((seq=='a'));   % find all of the letter 'a' and replace with 1.
s1g=double((seq=='g'));   % find all of the letter 'g' and replace with 1.
s1t=double((seq=='t'));   % find all of the letter 't' and replace with 1.
s1c=double((seq=='c'));   % find all of the letter 'c' and replace with 1.

fclose(f);
%% Open file and read sequence A
f1='seqa.txt'; 
f=fopen(f1);        
seq=fread(f,'*char')';

if (ispc | ismac)
    eol = 13;   % pc's and mac's use 13 (carriage return) for end of line
else
    eol = 10;   % others use 10 (line feed)
end

% Find beginning of DNA sequence and read to end-of-line (cr='13')
x=findstr(seq,'ORIGIN');   seq=seq(x:end);
x=findstr(seq,eol);        seq=seq(x:end);

seq=seq(isletter(seq));   % just take letters (drop numbers and spaces).
s2a=double((seq=='a'));   % find all of the letter 'a' and replace with 1.
s2g=double((seq=='g'));   % find all of the letter 'g' and replace with 1.
s2t=double((seq=='t'));   % find all of the letter 't' and replace with 1.
s2c=double((seq=='c'));   % find all of the letter 'c' and replace with 1.

fclose(f);

%% Open file and read sequence B

% Open file and read sequence.
f1='seqb.txt'; 
f=fopen(f1);        
seq=fread(f,'*char')';

if (ispc | ismac)
    eol = 13;   % pc's and mac's use 13 (carriage return) for end of line
else
    eol = 10;   % others use 10 (line feed)
end

% Find beginning of DNA sequence and read to end-of-line (cr='13')
x=findstr(seq,'ORIGIN');   seq=seq(x:end);
x=findstr(seq,eol);        seq=seq(x:end);

seq=seq(isletter(seq));   % just take letters (drop numbers and spaces).
s3a=double((seq=='a'));   % find all of the letter 'a' and replace with 1.
s3g=double((seq=='g'));   % find all of the letter 'g' and replace with 1.
s3t=double((seq=='t'));   % find all of the letter 't' and replace with 1.
s3c=double((seq=='c'));   % find all of the letter 'c' and replace with 1.

fclose(f);

%% Perform Correlations
r1a=xcorr(s1a,s2a);
r1a=r1a*(1/min(length(s1a), length(s2a)));

r2a=xcorr(s1a,s3a);
r2a=r2a*(1/min(length(s1a), length(s3a)));

%subplot(4,1,1)
%stem(r1a)

r1g=xcorr(s1g,s2g);
r1g=r1g*(1/min(length(s1g), length(s2g)));
r2g=xcorr(s1g,s3g);
r2g=r2g*(1/min(length(s1g), length(s3g)));

%subplot(4,1,2)
%stem(r1g)

r1c=xcorr(s1c,s2c);
r1c=r1c*(1/min(length(s1c), length(s2c)));
r2c=xcorr(s1c,s3c);
r2c=r2c*(1/min(length(s1c), length(s3c)));

%subplot(4,1,3)
%stem(r1c)

r1t=xcorr(s1t,s2t);
r1t=r1t*(1/min(length(s1t), length(s2t)));
r2t=xcorr(s1t,s3t);
r2t=r2t*(1/min(length(s1t), length(s3t)));

%subplot(4,1,4)
%stem(r1t)

%% Sum Correlations and Plot
seq1corr=r1a+r1t+r1g+r1c;
seq2corr=r2a+r2t+r2g+r2c;

subplot(2,1,1);
stem(seq1corr)
title('Correlation of Human Hemoglobin vs. SeqA.txt');
subplot(2,1,2);
stem(seq2corr)
title('Correlation of Human Hemoglobin vs. SeqB.txt');
hold on