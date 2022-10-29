function c=melcepst(s,fs,w,nc,p,n,inc,fl,fh)
%function c=melcepst(s,fs)
%usage by MARINI 12122012:   c=melcepst(s,250,'M',12,14,1024,204);  //ok!!
%usage by MARINI 09012013:   c=melcepst(s,250,'M',12,14,1024,256);

% s = spectrogram(x,window,noverlap,F,fs);
% w_c3=spectrogram(f10,1024,768,256);
% wt_c3=spectrogram(emoNo,64,50,1024);  assumption 
%usage by MARINI 09012013 ver2 3.23pm:   c=melcepst(s,250,'M',12,14,64,14);
% theta   fct_f10=melcepst([f10(:,1),f10(:,3),f10(:,5),f10(:,7),f10(:,9),f10(:,11),f10(:,13),f10(:,15)],250,'M',12,14,64,14);
% alpha   fca_f10=melcepst([f10(:,2),f10(:,4),f10(:,6),f10(:,8),f10(:,10),f10(:,12),f10(:,14),f10(:,16)],250,'M',12,14,64,14);
% 
%MELCEPST Calculate the mel cepstrum of a signal C=(S,FS,W,NC,P,N,INC,FL,FH)
%
%
% Simple use: c=melcepst(s,fs)	% calculate mel cepstrum with 12 coefs, 256 sample frames
%				  c=melcepst(s,fs,'e0dD') % include log energy, 0th cepstral coef, delta and delta-delta coefs
%
% Inputs:
%     s	 speech signal
%     fs  sample rate in Hz (default 11025)
%     nc  number of cepstral coefficients excluding 0'th coefficient (default 12)
%     n   length of frame in samples (default power of 2 < (0.03*fs))
%     p   number of filters in filterbank (default: floor(3*log(fs)) = approx 2.1 per ocatave)
%     inc frame increment (default n/2)
%     fl  low end of the lowest filter as a fraction of fs (default = 0)
%     fh  high end of highest filter as a fraction of fs (default = 0.5)
%
%		w   any sensible combination of the following:
%
%				'R'  rectangular window in time domain
%				'N'	Hanning window in time domain
%				'M'	Hamming window in time domain (default)
%
%		      't'  triangular shaped filters in mel domain (default)
%		      'n'  hanning shaped filters in mel domain
%		      'm'  hamming shaped filters in mel domain
%
%				'p'	filters act in the power domain
%				'a'	filters act in the absolute magnitude domain (default)
%
%			   '0'  include 0'th order cepstral coefficient
%				'E'  include log energy
%				'd'	include delta coefficients (dc/dt)
%				'D'	include delta-delta coefficients (d^2c/dt^2)
%
%		      'z'  highest and lowest filters taper down to zero (default)
%		      'y'  lowest filter remains at 1 down to 0 frequency and
%			   	  highest filter remains at 1 up to nyquist freqency
%
%		       If 'ty' or 'ny' is specified, the total power in the fft is preserved.
%
% Outputs:	c     mel cepstrum output: one frame per row. Log energy, if requested, is the
%                 first element of each row followed by the delta and then the delta-delta
%                 coefficients.
%

% BUGS: (1) should have power limit as 1e-16 rather than 1e-6 (or possibly a better way of choosing this)
%           and put into VOICEBOX
%       (2) get rdct to change the data length (properly) instead of doing it explicitly (wrongly)

%      Copyright (C) Mike Brookes 1997
%      Version: $Id: melcepst.m 1683 2012-03-22 12:55:45Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<2 fs=11025; end
if nargin<3 w='M'; end
if nargin<4 nc=12; end
if nargin<5 p=floor(3*log(fs)); end
if nargin<6 n=pow2(floor(log2(0.03*fs))); end
if nargin<9
   fh=0.5;   
   if nargin<8
     fl=0;
     if nargin<7
        inc=floor(n/2);
     end
  end
end

if isempty(w)
   w='M';
end
if any(w=='R')
   z=enframe(s,n,inc);
elseif any (w=='N')
   z=enframe(s,hanning(n),inc);
else
   z=enframe(s,hamming(n),inc);
end
f=rfft(z.');
[m,a,b]=melbankm(p,n,fs,fl,fh,w);
pw=f(a:b,:).*conj(f(a:b,:));
pth=max(pw(:))*1E-20;
if any(w=='p')
   y=log(max(m*pw,pth));
else
   ath=sqrt(pth);
   y=log(max(m*abs(f(a:b,:)),ath));
end
c=rdct(y).';
nf=size(c,1);
nc=nc+1;
if p>nc
   c(:,nc+1:end)=[];
elseif p<nc
   c=[c zeros(nf,nc-p)];
end
if ~any(w=='0')
   c(:,1)=[];
   nc=nc-1;
end
if any(w=='E')
   c=[log(sum(pw)).' c];
   nc=nc+1;
end

% calculate derivative

if any(w=='D')
  vf=(4:-1:-4)/60;
  af=(1:-1:-1)/2;
  ww=ones(5,1);
  cx=[c(ww,:); c; c(nf*ww,:)];
  vx=reshape(filter(vf,1,cx(:)),nf+10,nc);
  vx(1:8,:)=[];
  ax=reshape(filter(af,1,vx(:)),nf+2,nc);
  ax(1:2,:)=[];
  vx([1 nf+2],:)=[];
  if any(w=='d')
     c=[c vx ax];
  else
     c=[c ax];
  end
elseif any(w=='d')
  vf=(4:-1:-4)/60;
  ww=ones(4,1);
  cx=[c(ww,:); c; c(nf*ww,:)];
  vx=reshape(filter(vf,1,cx(:)),nf+8,nc);
  vx(1:8,:)=[];
  c=[c vx];
end
 
if nargout<1
   [nf,nc]=size(c);
   t=((0:nf-1)*inc+(n-1)/2)/fs;
   ci=(1:nc)-any(w=='0')-any(w=='E');
   imh = imagesc(t,ci,c.');
   axis('xy');
   xlabel('Time (s)');
   ylabel('Mel-cepstrum coefficient');
	map = (0:63)'/63;
	colormap([map map map]);
	colorbar;
end


function [f,t,w]=enframe(x,win,inc,m)
%ENFRAME split signal up into (overlapping) frames: one per row. [F,T]=(X,WIN,INC)
%
% Usage:  (1) f=enframe(x,n)     % split into frames of length n
%
%         (2) f=enframe(x,hamming(n,'periodic'),n/4)     % use a 75% overlapped Hamming window of length n
%
%         (3) frequency domain frame-based processing:
% 
%             S=...;                              % input signal
%             OV=2;                               % overlap factor of 2 (4 is also often used)
%             INC=20;                             % set frame increment in samples
%             NW=INC*OV;                          % DFT window length
%             W=sqrt(hamming(NW,'periodic'));     % omit sqrt if OV=4
%             W=W/sqrt(sum(W(1:INC:NW).^2));      % normalize window
%             F=rfft(enframe(S,W,INC),NW,2);      % do STFT: one row per time frame, +ve frequencies only
%             ... process frames ...
%             X=overlapadd(irfft(F,NW,2),W,INC);  % reconstitute the time waveform (omit "X=" to plot waveform)
%
%  Inputs:   x    input signal
%          win    window or window length in samples
%          inc    frame increment in samples
%            m    mode input:
%                  'z'  zero pad to fill up final frame
%                  'r'  reflect last few samples for final frame
%                  'A'  calculate window times as the centre of mass
%                  'E'  calculate window times as the centre of energy
%
% Outputs:   f    enframed data - one frame per row
%            t    fractional time in samples at the centre of each frame
%            w    window function used
%
% By default, the number of frames will be rounded down to the nearest
% integer and the last few samples of x() will be ignored unless its length
% is lw more than a multiple of inc. If the 'z' or 'r' options are given,
% the number of frame will instead be rounded up and no samples will be ignored.
%
% Example of frame-based processing:
%          INC=20       						% set frame increment in samples
%          NW=INC*2     						% oversample by a factor of 2 (4 is also often used)
%          S=cos((0:NW*7)*6*pi/NW);				% example input signal
%          W=sqrt(hamming(NW),'periodic'));  	% sqrt hamming window of period NW
%          F=enframe(S,W,INC);               	% split into frames
%          ... process frames ...
%          X=overlapadd(F,W,INC);               % reconstitute the time waveform (omit "X=" to plot waveform)

% Bugs/Suggestions:
%  (1) Possible additional mode options:
%        'u'  modify window for first and last few frames to ensure WOLA
%        'a'  normalize window to give a mean of unity after overlaps
%        'e'  normalize window to give an energy of unity after overlaps
%        'wm' use Hamming window
%        'wn' use Hanning window
%        'x'  include all frames that include any of the x samples

%	   Copyright (C) Mike Brookes 1997-2012
%      Version: $Id: enframe.m 2470 2012-11-02 15:27:24Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx=length(x(:));
if nargin<2 || isempty(win)
    win=nx;
end
if nargin<4 || isempty(m)
    m='';
end
nwin=length(win);
if nwin == 1
    lw = win;
    w = ones(1,lw);
else
    lw = nwin;
    w = win(:)';
end
if (nargin < 3) || isempty(inc)
    inc = lw;
end
nli=nx-lw+inc;
nf = fix((nli)/inc);
na=nli-inc*nf;
f=zeros(nf,lw);
indf= inc*(0:(nf-1)).';
inds = (1:lw);
f(:) = x(indf(:,ones(1,lw))+inds(ones(nf,1),:));
if nargin>3 && (any(m=='z') || any(m=='r')) && na>0
    if any(m=='r')
        ix=1+mod(nx-na:nx-na+lw-1,2*nx);
        f(nf+1,:)=x(ix+(ix>nx).*(2*nx+1-2*ix));
    else
        f(nf+1,1:na)=x(1+nx-na:nx);
    end
    nf=size(f,1);
end
if (nwin > 1)   % if we have a non-unity window
    f = f .* w(ones(nf,1),:);
end
if nargout>1
    if any(m=='E')
        t0=sum((1:lw).*w.^2)/sum(w.^2);
    elseif any(m=='E')
        t0=sum((1:lw).*w)/sum(w);
    else
        t0=(1+lw)/2;
    end
    t=t0+inc*(0:(nf-1)).';
end


function y=rfft(x,n,d)
%RFFT     Calculate the DFT of real data Y=(X,N,D)
% Data is truncated/padded to length N if specified.
%   N even:	(N+2)/2 points are returned with
% 			the first and last being real
%   N odd:	(N+1)/2 points are returned with the
% 			first being real
% In all cases fix(1+N/2) points are returned
% D is the dimension along which to do the DFT



%      Copyright (C) Mike Brookes 1998
%      Version: $Id: rfft.m 713 2011-10-16 14:45:43Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=size(x);
if prod(s)==1
    y=x
else
    if nargin <3 || isempty(d)
        d=find(s>1,1);
        if nargin<2
            n=s(d);
        end
    end
    if isempty(n) 
        n=s(d);
    end
    y=fft(x,n,d);
    y=reshape(y,prod(s(1:d-1)),n,prod(s(d+1:end))); 
    s(d)=1+fix(n/2);
    y(:,s(d)+1:end,:)=[];
    y=reshape(y,s);
end

function y=rdct(x,n,a,b)
%RDCT     Discrete cosine transform of real data Y=(X,N,A,B)
% Data is truncated/padded to length N.
%
% This routine is equivalent to multiplying by the matrix
%
%   rdct(eye(n)) = diag([sqrt(2)*B/A repmat(2/A,1,n-1)]) * cos((0:n-1)'*(0.5:n)*pi/n)
%
% Default values of the scaling factors are A=sqrt(2N) and B=1 which
% results in an orthogonal matrix. Other common values are A=1 or N and/or B=1 or sqrt(2). 
% If b~=1 then the columns are no longer orthogonal.
%
% see IRDCT for the inverse transform

% BUG: in line 51 we should do chopping after transform and not before



%      Copyright (C) Mike Brookes 1998
%      Version: $Id: rdct.m 713 2011-10-16 14:45:43Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fl=size(x,1)==1;
if fl x=x(:); end
[m,k]=size(x);
if nargin<2 n=m;
end
if nargin<4 b=1;  
    if nargin<3 a=sqrt(2*n);
    end
    end
if n>m x=[x; zeros(n-m,k)];
elseif n<m x(n+1:m,:)=[];
end

x=[x(1:2:n,:); x(2*fix(n/2):-2:2,:)];
z=[sqrt(2) 2*exp((-0.5i*pi/n)*(1:n-1))].';
y=real(fft(x).*z(:,ones(1,k)))/a;
y(1,:)=y(1,:)*b;
if fl y=y.'; end


function [x,mc,mn,mx]=melbankm(p,n,fs,fl,fh,w)
%MELBANKM determine matrix for a mel/erb/bark-spaced filterbank [X,MN,MX]=(P,N,FS,FL,FH,W)
%
% Inputs:
%       p   number of filters in filterbank or the filter spacing in k-mel/bark/erb [ceil(4.6*log10(fs))]
%		n   length of fft
%		fs  sample rate in Hz
%		fl  low end of the lowest filter as a fraction of fs [default = 0]
%		fh  high end of highest filter as a fraction of fs [default = 0.5]
%		w   any sensible combination of the following:
%             'b' = bark scale instead of mel
%             'e' = erb-rate scale
%             'l' = log10 Hz frequency scale
%             'f' = linear frequency scale
%
%             'c' = fl/fh specify centre of low and high filters
%             'h' = fl/fh are in Hz instead of fractions of fs
%             'H' = fl/fh are in mel/erb/bark/log10
%
%		      't' = triangular shaped filters in mel/erb/bark domain (default)
%		      'n' = hanning shaped filters in mel/erb/bark domain
%		      'm' = hamming shaped filters in mel/erb/bark domain
%
%		      'z' = highest and lowest filters taper down to zero [default]
%		      'y' = lowest filter remains at 1 down to 0 frequency and
%			        highest filter remains at 1 up to nyquist freqency
%
%             'u' = scale filters to sum to unity
%
%             's' = single-sided: do not double filters to account for negative frequencies
%
%             'g' = plot idealized filters [default if no output arguments present]
%
% Note that the filter shape (triangular, hamming etc) is defined in the mel (or erb etc) domain.
% Some people instead define an asymmetric triangular filter in the frequency domain.
%
%		       If 'ty' or 'ny' is specified, the total power in the fft is preserved.
%
% Outputs:	x     a sparse matrix containing the filterbank amplitudes
%		          If the mn and mx outputs are given then size(x)=[p,mx-mn+1]
%                 otherwise size(x)=[p,1+floor(n/2)]
%                 Note that the peak filter values equal 2 to account for the power
%                 in the negative FFT frequencies.
%           mc    the filterbank centre frequencies in mel/erb/bark
%		    mn    the lowest fft bin with a non-zero coefficient
%		    mx    the highest fft bin with a non-zero coefficient
%                 Note: you must specify both or neither of mn and mx.
%
% Examples of use:
%
% (a) Calcuate the Mel-frequency Cepstral Coefficients
%
%       f=rfft(s);			        % rfft() returns only 1+floor(n/2) coefficients
%		x=melbankm(p,n,fs);	        % n is the fft length, p is the number of filters wanted
%		z=log(x*abs(f).^2);         % multiply x by the power spectrum
%		c=dct(z);                   % take the DCT
%
% (b) Calcuate the Mel-frequency Cepstral Coefficients efficiently
%
%       f=fft(s);                        % n is the fft length, p is the number of filters wanted
%       [x,mc,na,nb]=melbankm(p,n,fs);   % na:nb gives the fft bins that are needed
%       z=log(x*(f(na:nb)).*conj(f(na:nb)));
%
% (c) Plot the calculated filterbanks
%
%      plot((0:floor(n/2))*fs/n,melbankm(p,n,fs)')   % fs=sample frequency
%
% (d) Plot the idealized filterbanks (without output sampling)
%
%      melbankm(p,n,fs);
%
% References:
%
% [1] S. S. Stevens, J. Volkman, and E. B. Newman. A scale for the measurement
%     of the psychological magnitude of pitch. J. Acoust Soc Amer, 8: 185–19, 1937.
% [2] S. Davis and P. Mermelstein. Comparison of parametric representations for
%     monosyllabic word recognition in continuously spoken sentences.
%     IEEE Trans Acoustics Speech and Signal Processing, 28 (4): 357–366, Aug. 1980.


%      Copyright (C) Mike Brookes 1997-2009
%      Version: $Id: melbankm.m 713 2011-10-16 14:45:43Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note "FFT bin_0" assumes DC = bin 0 whereas "FFT bin_1" means DC = bin 1

if nargin < 6
    w='tz'; % default options
    if nargin < 5
        fh=0.5; % max freq is the nyquist
        if nargin < 4
            fl=0; % min freq is DC
        end
    end
end
sfact=2-any(w=='s');   % 1 if single sided else 2
wr=' ';   % default warping is mel
for i=1:length(w)
    if any(w(i)=='lebf');
        wr=w(i);
    end
end
if any(w=='h') || any(w=='H')
    mflh=[fl fh];
else
    mflh=[fl fh]*fs;
end
if ~any(w=='H')
    switch wr
                    case 'f'       % no transformation
        case 'l'
            if fl<=0
                error('Low frequency limit must be >0 for l option');
            end
            mflh=log10(mflh);       % convert frequency limits into log10 Hz
        case 'e'
            mflh=frq2erb(mflh);       % convert frequency limits into erb-rate
        case 'b'
            mflh=frq2bark(mflh);       % convert frequency limits into bark
        otherwise
            mflh=frq2mel(mflh);       % convert frequency limits into mel
    end
end
melrng=mflh*(-1:2:1)';          % mel range
fn2=floor(n/2);     % bin index of highest positive frequency (Nyquist if n is even)
if isempty(p)
    p=ceil(4.6*log10(fs));         % default number of filters
end
if any(w=='c')              % c option: specify fiter centres not edges
if p<1
    p=round(melrng/(p*1000))+1;
end
melinc=melrng/(p-1);
mflh=mflh+(-1:2:1)*melinc;
else
    if p<1
    p=round(melrng/(p*1000))-1;
end
melinc=melrng/(p+1);
end

%
% Calculate the FFT bins corresponding to [filter#1-low filter#1-mid filter#p-mid filter#p-high]
%
switch wr
        case 'f'
        blim=(mflh(1)+[0 1 p p+1]*melinc)*n/fs;
    case 'l'
        blim=10.^(mflh(1)+[0 1 p p+1]*melinc)*n/fs;
    case 'e'
        blim=erb2frq(mflh(1)+[0 1 p p+1]*melinc)*n/fs;
    case 'b'
        blim=bark2frq(mflh(1)+[0 1 p p+1]*melinc)*n/fs;
    otherwise
        blim=mel2frq(mflh(1)+[0 1 p p+1]*melinc)*n/fs;
end
mc=mflh(1)+(1:p)*melinc;    % mel centre frequencies
b1=floor(blim(1))+1;            % lowest FFT bin_0 required might be negative)
b4=min(fn2,ceil(blim(4))-1);    % highest FFT bin_0 required
%
% now map all the useful FFT bins_0 to filter1 centres
%
switch wr
        case 'f'
        pf=((b1:b4)*fs/n-mflh(1))/melinc;
    case 'l'
        pf=(log10((b1:b4)*fs/n)-mflh(1))/melinc;
    case 'e'
        pf=(frq2erb((b1:b4)*fs/n)-mflh(1))/melinc;
    case 'b'
        pf=(frq2bark((b1:b4)*fs/n)-mflh(1))/melinc;
    otherwise
        pf=(frq2mel((b1:b4)*fs/n)-mflh(1))/melinc;
end
%
%  remove any incorrect entries in pf due to rounding errors
%
if pf(1)<0
    pf(1)=[];
    b1=b1+1;
end
if pf(end)>=p+1
    pf(end)=[];
    b4=b4-1;
end
fp=floor(pf);                  % FFT bin_0 i contributes to filters_1 fp(1+i-b1)+[0 1]
pm=pf-fp;                       % multiplier for upper filter
k2=find(fp>0,1);   % FFT bin_1 k2+b1 is the first to contribute to both upper and lower filters
k3=find(fp<p,1,'last');  % FFT bin_1 k3+b1 is the last to contribute to both upper and lower filters
k4=numel(fp); % FFT bin_1 k4+b1 is the last to contribute to any filters
if isempty(k2)
    k2=k4+1;
end
if isempty(k3)
    k3=0;
end
if any(w=='y')          % preserve power in FFT
    mn=1; % lowest fft bin required (1 = DC)
    mx=fn2+1; % highest fft bin required (1 = DC)
    r=[ones(1,k2+b1-1) 1+fp(k2:k3) fp(k2:k3) repmat(p,1,fn2-k3-b1+1)]; % filter number_1
    c=[1:k2+b1-1 k2+b1:k3+b1 k2+b1:k3+b1 k3+b1+1:fn2+1]; % FFT bin1
    v=[ones(1,k2+b1-1) pm(k2:k3) 1-pm(k2:k3) ones(1,fn2-k3-b1+1)];
else
    r=[1+fp(1:k3) fp(k2:k4)]; % filter number_1
    c=[1:k3 k2:k4]; % FFT bin_1 - b1
    v=[pm(1:k3) 1-pm(k2:k4)];
    mn=b1+1; % lowest fft bin_1
    mx=b4+1;  % highest fft bin_1
end
if b1<0
    c=abs(c+b1-1)-b1+1;     % convert negative frequencies into positive
end
% end
if any(w=='n')
    v=0.5-0.5*cos(v*pi);      % convert triangles to Hanning
elseif any(w=='m')
    v=0.5-0.46/1.08*cos(v*pi);  % convert triangles to Hamming
end
if sfact==2  % double all except the DC and Nyquist (if any) terms
    msk=(c+mn>2) & (c+mn<n-fn2+2);  % there is no Nyquist term if n is odd
    v(msk)=2*v(msk);
end
%
% sort out the output argument options
%
if nargout > 2
    x=sparse(r,c,v);
    if nargout == 3     % if exactly three output arguments, then
        mc=mn;          % delete mc output for legacy code compatibility
        mn=mx;
    end
else
    x=sparse(r,c+mn-1,v,p,1+fn2);
end
if any(w=='u')
    sx=sum(x,2);
    x=x./repmat(sx+(sx==0),1,size(x,2));
end
%
% plot results if no output arguments or g option given
%
if ~nargout || any(w=='g') % plot idealized filters
    ng=201;     % 201 points
    me=mflh(1)+(0:p+1)'*melinc;
    switch wr
                case 'f'
            fe=me; % defining frequencies
            xg=repmat(linspace(0,1,ng),p,1).*repmat(me(3:end)-me(1:end-2),1,ng)+repmat(me(1:end-2),1,ng);
        case 'l'
            fe=10.^me; % defining frequencies
            xg=10.^(repmat(linspace(0,1,ng),p,1).*repmat(me(3:end)-me(1:end-2),1,ng)+repmat(me(1:end-2),1,ng));
        case 'e'
            fe=erb2frq(me); % defining frequencies
            xg=erb2frq(repmat(linspace(0,1,ng),p,1).*repmat(me(3:end)-me(1:end-2),1,ng)+repmat(me(1:end-2),1,ng));
        case 'b'
            fe=bark2frq(me); % defining frequencies
            xg=bark2frq(repmat(linspace(0,1,ng),p,1).*repmat(me(3:end)-me(1:end-2),1,ng)+repmat(me(1:end-2),1,ng));
        otherwise
            fe=mel2frq(me); % defining frequencies
            xg=mel2frq(repmat(linspace(0,1,ng),p,1).*repmat(me(3:end)-me(1:end-2),1,ng)+repmat(me(1:end-2),1,ng));
    end

    v=1-abs(linspace(-1,1,ng));
    if any(w=='n')
        v=0.5-0.5*cos(v*pi);      % convert triangles to Hanning
    elseif any(w=='m')
        v=0.5-0.46/1.08*cos(v*pi);  % convert triangles to Hamming
    end
    v=v*sfact;  % multiply by 2 if double sided
    v=repmat(v,p,1);
    if any(w=='y')  % extend first and last filters
        v(1,xg(1,:)<fe(2))=sfact;
        v(end,xg(end,:)>fe(p+1))=sfact;
    end
    if any(w=='u') % scale to unity sum
        dx=(xg(:,3:end)-xg(:,1:end-2))/2;
        dx=dx(:,[1 1:ng-2 ng-2]);
        vs=sum(v.*dx,2);
        v=v./repmat(vs+(vs==0),1,ng)*fs/n;
    end
    plot(xg',v','b');
    set(gca,'xlim',[fe(1) fe(end)]);
    xlabel(['Frequency (' xticksi 'Hz)']);
end

function [mel,mr] = frq2mel(frq)
%FRQ2ERB  Convert Hertz to Mel frequency scale MEL=(FRQ)
%	[mel,mr] = frq2mel(frq) converts a vector of frequencies (in Hz)
%	to the corresponding values on the Mel scale which corresponds
%	to the perceived pitch of a tone.
%   mr gives the corresponding gradients in Hz/mel.

%	The relationship between mel and frq is given by [1]:
%
%	m = ln(1 + f/700) * 1000 / ln(1+1000/700)
%
%  	This means that m(1000) = 1000
%
%	References:
%
%     [1] J. Makhoul and L. Cosell. "Lpcw: An lpc vocoder with
%         linear predictive spectral warping", In Proc IEEE Intl
%         Conf Acoustics, Speech and Signal Processing, volume 1,
%         pages 466–469, 1976. doi: 10.1109/ICASSP.1976.1170013.
%	  [2] S. S. Stevens & J. Volkman "The relation of pitch to
%		  frequency", American J of Psychology, V 53, p329 1940
%	  [3] C. G. M. Fant, "Acoustic description & classification
%		  of phonetic units", Ericsson Tchnics, No 1 1959
%		  (reprinted in "Speech Sounds & Features", MIT Press 1973)
%	  [4] S. B. Davis & P. Mermelstein, "Comparison of parametric
%		  representations for monosyllabic word recognition in
%		  continuously spoken sentences", IEEE ASSP, V 28,
%		  pp 357-366 Aug 1980
%	  [5] J. R. Deller Jr, J. G. Proakis, J. H. L. Hansen,
%		  "Discrete-Time Processing of Speech Signals", p380,
%		  Macmillan 1993

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: frq2mel.m 1874 2012-05-25 15:41:53Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent k
if isempty(k)
    k=1000/log(1+1000/700); %  1127.01048
end
af=abs(frq);
mel = sign(frq).*log(1+af/700)*k;
mr=(700+af)/k;
if ~nargout
    plot(frq,mel,'-x');
    xlabel(['Frequency (' xticksi 'Hz)']);
    ylabel(['Frequency (' yticksi 'Mel)']);
end
function [frq,mr] = mel2frq(mel)
%MEL2FRQ  Convert Mel frequency scale to Hertz FRQ=(MEL)
%	frq = mel2frq(mel) converts a vector of Mel frequencies
%	to the corresponding real frequencies.
%   mr gives the corresponding gradients in Hz/mel.
%	The Mel scale corresponds to the perceived pitch of a tone

%	The relationship between mel and frq is given by [1]:
%
%	m = ln(1 + f/700) * 1000 / ln(1+1000/700)
%
%  	This means that m(1000) = 1000
%
%	References:
%
%     [1] J. Makhoul and L. Cosell. "Lpcw: An lpc vocoder with
%         linear predictive spectral warping", In Proc IEEE Intl
%         Conf Acoustics, Speech and Signal Processing, volume 1,
%         pages 466–469, 1976. doi: 10.1109/ICASSP.1976.1170013.
%	  [2] S. S. Stevens & J. Volkman "The relation of pitch to
%		  frequency", American J of Psychology, V 53, p329 1940
%	  [3] C. G. M. Fant, "Acoustic description & classification
%		  of phonetic units", Ericsson Tchnics, No 1 1959
%		  (reprinted in "Speech Sounds & Features", MIT Press 1973)
%	  [4] S. B. Davis & P. Mermelstein, "Comparison of parametric
%		  representations for monosyllabic word recognition in
%		  continuously spoken sentences", IEEE ASSP, V 28,
%		  pp 357-366 Aug 1980
%	  [5] J. R. Deller Jr, J. G. Proakis, J. H. L. Hansen,
%		  "Discrete-Time Processing of Speech Signals", p380,
%		  Macmillan 1993

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: mel2frq.m 1874 2012-05-25 15:41:53Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent k
if isempty(k)
    k=1000/log(1+1000/700); % 1127.01048
end
frq=700*sign(mel).*(exp(abs(mel)/k)-1);
mr=(700+abs(frq))/k;
if ~nargout
    plot(mel,frq,'-x');
    xlabel(['Frequency (' xticksi 'Mel)']);
    ylabel(['Frequency (' yticksi 'Hz)']);
end
