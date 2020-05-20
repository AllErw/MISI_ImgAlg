load('D:\Git_code\MISI_ImgAlg_MATLAB_demos\MATLAB demo code\LGU_test_data.mat');

orig = rf_data(1:6,1:7);
S = zeros(size(rf_data));

R = fft2(orig);
% R = fftshift(fft2(orig));

[Nt,Nsrc] = size(R);

S = zeros(size(R));
S(1:floor(Nt/2)+1 , :) = R(1:floor(Nt/2)+1 , :);

% Impose symmetry:
% S(floor(Nt/2)+2:Nt , 1) = conj(S(ceil(Nt/2):-1:2 , 1));
S(Nt:-1:floor(Nt/2)+2 , 1) = conj(S(2:ceil(Nt/2) , 1));
S(Nt:-1:floor(Nt/2)+2 , Nsrc:-1:2) = conj(S(2:ceil(Nt/2) , 2:Nsrc));

figure;colormap hot;
subplot(2,2,1);imagesc(real(R))
subplot(2,2,2);imagesc(real(S))
subplot(2,2,3);imagesc(orig);
subplot(2,2,4);imagesc(real(ifft2(S)));

