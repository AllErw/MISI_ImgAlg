function out = PCfunctChirp(in,dur,fmin,fmax,fsamp,Tukpar)

Nt = length(in);
Nchirp = ceil(dur*fsamp);
Nscan = size(in,2);

% Define chirp:
chirp = zeros(1,Nchirp);
for tCcnt = 0:Nchirp-1
    tchirp = tCcnt/fsamp;
    chirp(tCcnt+1) = sin(2*pi  *  ((fmax-fmin)/2/dur*tchirp + fmin)*tchirp);
end

% Apply Tukey window for mismatched filtering:
if Tukpar>0 && Tukpar<=1
    for tCcnt = 1:ceil(Nchirp*Tukpar/2)
        Tukwin = 0.5*(1+cos(2*pi * ((tCcnt-1)/(Nchirp-1)/Tukpar-0.5)));
        chirp(tCcnt) = chirp(tCcnt) * Tukwin;
        chirp(Nchirp-tCcnt+1) = chirp(Nchirp-tCcnt+1) * Tukwin;
    end
end

% Cross-correlation for pulse compression:
out = zeros(Nt,Nscan);
for scnt = 1:Nscan
    for tcnt = 1:Nt
        for tCcnt = 1:min(Nchirp,Nt-tcnt)
            out(tcnt,scnt) = out(tcnt,scnt) + in(tcnt-1+tCcnt,scnt)*chirp(tCcnt);
        end
    end
end

% out = real(ifft(    fft(in,Nt+Nchirp-1,1) .* (fft(chirp(end:-1:1)',Nt+Nchirp-1,1) * ones(1,Nscan))    ,[],1));
% out = out(Nchirp:end,:);