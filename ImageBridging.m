clc;

% Original Images are 256 pixels X 256 pixels.

fprintf('Reading Image \n');

COMPRESSION_PERCENT = 0.10; % Compressed Signal will be approximately
% n = 256^2 * COMPRESSION_PERCENT dimensional.
snr = .05;
percenterasures = .05;

Original_Image_Double = double(imread('Lena.bmp'));

fprintf('Performing Image Compression \n');

Compressed_Image_Double = fft(reshape(Original_Image_Double,[256*256,1]));
[S,I] = sort(abs(Compressed_Image_Double),'descend');
n = round(COMPRESSION_PERCENT*256*256);
LSC = Compressed_Image_Double(I(n+1:256*256));
Compressed_Image_Double(I(n+1:256*256)) = [];

N = 2*n;
Lsize = round(percenterasures*N);
L = [1:Lsize];
W = [Lsize+1:3*Lsize];
LC = setdiff(1:N,L);

f = Compressed_Image_Double;

fprintf('Creating Frames \n');

A = randn(N,n);
[A,~] = qr(A,0);

F = sqrt(N/n)*A';
G = (n/N)*F;

fprintf('Reconstructing Erasures \n');

FC = G' * f;
FC(L) = zeros(size(L'));
noise = randn(size(LC'));
noise = snr * noise ./ norm(noise) * norm(FC(LC));
FC(LC) = FC(LC) + noise;
f_R = F*FC;



FRCL = G(:,L)' * f_R;
FRCB = G(:,W)' * f_R;
C = pinv(F(:,L)'*G(:,W))*(F(:,L)'*G(:,L));
FC(L) = C' * (FC(W) - FRCB) + FRCL;
g = f_R + F(:,L) * FC(L);

fprintf('Plotting Images \n');

C_f = zeros(256*256,1); % Compressed Image.
C_g = zeros(256*256,1); % Reconstructed Image.

I1 = sort(I(1:n),'ascend');

C_f(I1) = f;
Uncompressed_f = ifft(C_f);
Uncompressed_f = reshape(Uncompressed_f,[256,256]);
J_f = uint8(Uncompressed_f);

C_g(I1(1:n)) = g;
Uncompressed_g = ifft(C_g);
Uncompressed_g = reshape(Uncompressed_g,[256,256]);
J_g = uint8(Uncompressed_g);

figure;

subplot(1,2,1);
imshow(J_f);
title('Compressed Image');

subplot(1,2,2);
imshow(J_g);
title('Reconstructed Image');
