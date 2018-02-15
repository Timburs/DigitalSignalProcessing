%% Open file and read human genome
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
x=findstr(seq,eol);        seq=seq(x:end);

seq=seq(isletter(seq));   % just take letters (drop numbers and spaces).
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
r2a=xcorr(s1a,s3a);

r1g=xcorr(s1g,s2g);
r2g=xcorr(s1g,s3g);

r1c=xcorr(s1c,s2c);
r2c=xcorr(s1c,s3c);

r1t=xcorr(s1t,s2t);
r2t=xcorr(s1t,s3t);

%% Sum Correlations and Plot
seq1corr=r1a+r1t+r1g+r1c;
seq2corr=r2a+r2t+r2g+r2c;

x=0:length(r1a)-1;
subplot(2,1,1); 
plot(x,seq1corr)
title('Genome Sequence 1');
axis([70000 146651 100 400])
subplot(2,1,2);
plot(x,seq2corr)
title('Genome Sequence 2');
axis([70000 146651 100 400])