% GenerateFilters.m
%
%	Methodology:
%	This program makes use of the Oct3dsgn function written by Christophe Couvreur
%	In order to get the poles to lie within the unit circle for lower frequencies,
%	the fitlers were designed at smapling frequencies lower than 44100 Hz. 
%	Therefore, lower frequency 1/3 Octave Bands have their spectrums truncated
%	according to the fs/2 of the adjusted sampling frequencies.
%
%	Syntax:
%	[H, err]=GenerateFilters(f, Fs)
%	
%	Variables:
%	INPUT
%	f		=	Frequency (Hz) of the sectrum to be filtered, used to determine df
%	
%	WORKING
%	df		=	Frequency resolution of spectrum (Hz) 
%	Fs		=	Assumed Sampling Frequency (Hz) of Signal
%	FiltOrd	=	Order of filter to be designed
%	Fc		=	Center frequencies (Hz) for filters to be designed
%	q		=	Integer resampling coefficient
%	FsNew	=	New Sampling Frequency (Hz)
%	B#		=	Filter Polynmomial Coefficient
%	A#		=	Filter Polynomial Coefficient
%	fo		=	Starting Frequency (Hz)
%	f2		=	Filter Frequency Vector (Hz)
%	H1		=	Temporary storage for Frequency Response of Filter
%	H		=	Fitler Set
%
%	OUTPUT
%	H		=	Frequency Responses of the Filter Set
%	err   = 	Value for an error return
%               0 = no error
%               1 = unkown error


%	Author:	Aaron Hastings, Herrick Labs, Purdue University
%	Date Started:	31 Oct 00
%	Last Revision:	20 March 01 (Added Fs input)
%	Status: No Known Bugs (Just "Features"  :-)

function[H, err]=GenerateFilters(f, Fs)

%%	Begin function

err=1;

%Fs=44100;

Fc=[25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, ...
      800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, ...
      10000, 12500, 16000];

%%	Filters with fc < 220 will be resampled and then zero padded

%%	Prototype
%%-----------------------------------------------------------------------------
%	q=4; 	
%	FsNew=Fs/q;
%	[B2,A2] = Oct3dsgn(Fc(1),FsNew);
%	df=f(2)-f(1);
%	fo=0;
%	f2=fo:df:FsNew/2;
%	H2=freqz(B2,A2,2*pi*f2/FsNew)';
%	H2sq=abs(H2).^2;
%	FiltLevel(1)=10*log10(sum((10.^(Yxx(1,1:4411)/10)).*(abs(H2.^2))'))
%	loglog(f2,abs(H2.^2))
%%-----------------------------------------------------------------------------

df=f(2)-f(1);
fo=0;
FiltOrd=3;

%%	Filter resampled for fc = 25 Hz.
ink=1;
q=16; 	
FsNew=Fs/q;
[B,A] = Oct3dsgn(Fc(ink),FsNew,FiltOrd);
f2=fo:df:FsNew/2;
H1=freqz(B,A,2*pi*f2/FsNew)';
H(ink,:)=[H1' (10e-37)*ones(1,length(f)-length(H1))];

%%	Filter resampled for fc = 31.5 Hz.
ink=2;
q=8; 	
FsNew=Fs/q;
[B,A] = Oct3dsgn(Fc(ink),FsNew,FiltOrd);
f2=fo:df:FsNew/2;
H1=freqz(B,A,2*pi*f2/FsNew)';
stop=length(H1);
H(ink,:)=[H1' (10e-37)*ones(1,length(f)-length(H1))];

%%	Filter resampled for fc = 40 Hz.
ink=3;
q=8; 	
FsNew=Fs/q;
[B,A] = Oct3dsgn(Fc(ink),FsNew,FiltOrd);
f2=fo:df:FsNew/2;
H1=freqz(B,A,2*pi*f2/FsNew)';
H(ink,:)=[H1' (10e-37)*ones(1,length(f)-length(H1))];

%%	Filter resampled for fc = 50 Hz.
ink=4;
q=8; 	
FsNew=Fs/q;
[B,A] = Oct3dsgn(Fc(ink),FsNew,FiltOrd);
f2=fo:df:FsNew/2;
H1=freqz(B,A,2*pi*f2/FsNew)';
H(ink,:)=[H1' (10e-37)*ones(1,length(f)-length(H1))];

%%	Filter resampled for fc = 63 Hz.
ink=5;
q=4; 	
FsNew=Fs/q;
[B,A] = Oct3dsgn(Fc(ink),FsNew,FiltOrd);
f2=fo:df:FsNew/2;
H1=freqz(B,A,2*pi*f2/FsNew)';
H(ink,:)=[H1' (10e-37)*ones(1,length(f)-length(H1))];

%%	Filter resampled for fc = 80 Hz.
ink=6;
q=4; 	
FsNew=Fs/q;
[B,A] = Oct3dsgn(Fc(ink),FsNew,FiltOrd);
f2=fo:df:FsNew/2;
H1=freqz(B,A,2*pi*f2/FsNew)';
H(ink,:)=[H1' (10e-37)*ones(1,length(f)-length(H1))];

%%	Filter resampled for fc = 100 Hz.
ink=7;
q=4; 	
FsNew=Fs/q;
[B,A] = Oct3dsgn(Fc(ink),FsNew,FiltOrd);
f2=fo:df:FsNew/2;
H1=freqz(B,A,2*pi*f2/FsNew)';
H(ink,:)=[H1' (10e-37)*ones(1,length(f)-length(H1))];

%%	Filter resampled for fc = 125 Hz.
ink=8;
q=2; 	
FsNew=Fs/q;
[B,A] = Oct3dsgn(Fc(ink),FsNew,FiltOrd);
f2=fo:df:FsNew/2;
H1=freqz(B,A,2*pi*f2/FsNew)';
H(ink,:)=[H1' (10e-37)*ones(1,length(f)-length(H1))];

%%	Filter resampled for fc = 160 Hz.
ink=9;
q=2; 	
FsNew=Fs/q;
[B,A] = Oct3dsgn(Fc(ink),FsNew,FiltOrd);
f2=fo:df:FsNew/2;
H1=freqz(B,A,2*pi*f2/FsNew)';
H(ink,:)=[H1' (10e-37)*ones(1,length(f)-length(H1))];

%%	Filter resampled for fc = 200 Hz.
ink=10;
q=2; 	
FsNew=Fs/q;
[B,A] = Oct3dsgn(Fc(ink),FsNew,FiltOrd);
f2=fo:df:FsNew/2;
H1=freqz(B,A,2*pi*f2/FsNew)';
H(ink,:)=[H1' (10e-37)*ones(1,length(f)-length(H1))];


if Fs==44100;
   for ink=11:29
      q=1;
      FsNew=Fs/q;
      [B,A] = Oct3dsgn(Fc(ink),FsNew,FiltOrd);
      f2=fo:df:FsNew/2;
      H1=freqz(B,A,2*pi*f2/FsNew)';
      H(ink,:)=[H1' (10e-37)*ones(1,length(f)-length(H1))];
   end   
elseif Fs==22050;
   %disp('Warning:  You are using a low sample frequency.  Sounds should not have content above 10 kHz!!!')
   %%  At some time I should get around to making some rectangular filters for these...
   for ink=11:26
      q=1;
      FsNew=Fs/q;
      [B,A] = Oct3dsgn(Fc(ink),FsNew,FiltOrd);
      f2=fo:df:FsNew/2;
      H1=freqz(B,A,2*pi*f2/FsNew)';
      H(ink,:)=[H1' (10e-37)*ones(1,length(f)-length(H1))];
   end   
   for k=27:29 H(ink,:)=[(10e-37)*ones(1,length(f))]; end
elseif Fs==16000
   for ink=11:25
      q=1;
      FsNew=Fs/q;
      [B,A] = Oct3dsgn(Fc(ink),FsNew,FiltOrd);
      f2=fo:df:FsNew/2;
      H1=freqz(B,A,2*pi*f2/FsNew)';
      H(ink,:)=[H1' (10e-37)*ones(1,length(f)-length(H1))];
   end
   for ink=26:29 H(ink,:)=[(10e-37)*ones(1,length(f))]; end
elseif Fs==8000
  for ink=11:22
      q=1;
      FsNew=Fs/q;
      [B,A] = Oct3dsgn(Fc(ink),FsNew,FiltOrd);
      f2=fo:df:FsNew/2;
      H1=freqz(B,A,2*pi*f2/FsNew)';
      H(ink,:)=[H1' (10e-37)*ones(1,length(f)-length(H1))];
  end
  for ink=23:29 H(ink,:)=[(10e-37)*ones(1,length(f))]; end
end   

err=0;

return
