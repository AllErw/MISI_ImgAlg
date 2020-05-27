% Parameters
fsamp = 250E6;
fmin = 2E6;
fmax = 35E6;
dur = 4E-6;

SNR = 2.5;%Inf;
Tukpar = 0;%0.25;      % Tukey parameter: 0 is off, 0.25 works well

Nt = 2048;
Nscan = 100;
taxis = (0:Nt-1) / fsamp;


% Set impulse response:
groundtruth = zeros(Nt,Nscan);
for scnt = 1:Nscan
    groundtruth(512 , scnt) = 4+2*(rand-.5);
end

% Define chirp:
t1 = linspace(1/fsamp,dur,dur*fsamp) - 1/fsamp;
chirp1 = sin(2*pi  *  ((fmax-fmin)/2/dur*t1 + fmin).*t1);

% Generate chirp-ed A-scans
Ascan = conv2(chirp1,1,groundtruth,'full');
Ascan = Ascan(1:length(groundtruth),:);
Ascan = Ascan + 2/SNR*(rand(size(Ascan)) - .5);

% Emulate amplifier depletion and "mismatch" the chirp filter:
MMfilt = chirp1.*tukeywin(length(chirp1),Tukpar)';

% Pulse compression:
tic;
PCmatlab = conv2(MMfilt(end:-1:1),1,Ascan,'full');
PCmatlab = PCmatlab(length(MMfilt):end , :);
toc;

tic;
PCfuncti = PCfunctChirp(Ascan,dur,fmin,fmax,fsamp,Tukpar);  
toc;

if ~libisloaded('ImgAlg')
disp('Loading library...');
loadlibrary('D:\Git_code\MISI_ImgAlg\x64\Release\MISI_ImgAlg.dll','D:\Git_code\MISI_ImgAlg\MISI_ImgAlg.h','alias','ImgAlg');
disp('Library loaded.');
end

tic;
PCdll = zeros(size(Ascan));
% [~,PCdll] = calllib('ImgAlg','PulseCompChirp',Ascan,dur,fmin,fmax,fsamp,Tukpar,Nt,Nscan,PCdll);
[~,~,PCdll] = calllib('ImgAlg','PulseCompAny',Ascan,Nt,Nscan,MMfilt,length(MMfilt),PCdll);
toc;


figure;
subplot(2,2,1);
imagesc(Ascan);

subplot(2,2,2);
imagesc(PCdll);

subplot(2,2,[3 4]);
plot( [ PCmatlab(512,:) ; PCfuncti(512,:) ; PCdll(512,:) ]' );


