clc;

% Original Images are 256 pixels X 256 pixels.

fprintf('Reading Image \n');

COMPRESSION_PERCENT = 0.20; % Compressed Signal will be approximately
% n = 256^2 * COMPRESSION_PERCENT dimensional.
percenterasures = .10;

Original_Image_Double = double(imread('Pepper.bmp'));

fprintf('Performing Image Compression \n')

Compressed_Image_Double = fft(reshape(Original_Image_Double,[256*256,1]));
[S,I] = sort(abs(Compressed_Image_Double),'descend');
n = round(COMPRESSION_PERCENT*256*256)
Compressed_Image_Double(I(n+1:256*256)) = [];

N = 2*n
snr = .10;

f = Compressed_Image_Double;

fprintf('Creating Frames \n');

F = (1/sqrt(n)) * randn(n,N);
S = F * F';
G = S \ F;

fprintf('Reconstructing Erasures \n');

figure;

C_f = zeros(256*256,1); % Compressed Image.
I1 = sort(I(1:n),'ascend');
C_f(I1) = f;
Uncompressed_f = ifft(C_f);
Uncompressed_f = reshape(Uncompressed_f,[256,256]);
J_f = uint8(Uncompressed_f);

subplot(1,5,1);
imshow(J_f);
title('Compressed Image');

Lsize = round(percenterasures*N);
L = [1:Lsize];
W = [Lsize+1:2*Lsize];
W1 = [Lsize+1:round(3*Lsize)];
LC = setdiff(1:N,L);

FC = G' * f;
noise = randn(size(LC'));
noise = noise / norm(noise) * snr * norm(FC(LC));
FC(LC) = FC(LC) + noise;
FC2 = FC;
FC(L) = zeros(size(L'));
FC1 = FC;
f_R = F*FC;

FRCL = G(:,L)' * f_R;
FRCB = G(:,W)' * f_R;
C = (F(:,L)'*G(:,W))\(F(:,L)'*G(:,L));
FC(L) = C' * (FC(W) - FRCB) + FRCL;
g = f_R + F(:,L) * FC(L);

FRCL1 = G(:,L)' * f_R;
FRCB1 = G(:,W1)' * f_R;
C1 = pinv(F(:,L)'*G(:,W1))*(F(:,L)'*G(:,L));
FC1(L) = C1' * (FC1(W1) - FRCB1) + FRCL1;
g1 = f_R + F(:,L) * FC1(L);

g2 = F * FC2;

C_g2 = zeros(256*256,1); % Just Noise.
C_g2(I1(1:n)) = g2;
Uncompressed_g2 = ifft(C_g2);
Uncompressed_g2 = reshape(Uncompressed_g2,[256,256]);
J_g2 = uint8(Uncompressed_g2);

C_f_R = zeros(256*256,1); % Erased Image with noise.
C_f_R(I1(1:n)) = f_R;
Uncompressed_f_R = ifft(C_f_R);
Uncompressed_f_R = reshape(Uncompressed_f_R,[256,256]);
J_f_R = uint8(Uncompressed_f_R);

C_g = zeros(256*256,1); % Reconstructed Image, Single Bridge.
C_g(I1(1:n)) = g;
Uncompressed_g = ifft(C_g);
Uncompressed_g = reshape(Uncompressed_g,[256,256]);
J_g = uint8(Uncompressed_g);

C_g1 = zeros(256*256,1); % Reconstructed Image Over-Bridged.
C_g1(I1(1:n)) = g1;
Uncompressed_g1 = ifft(C_g1);
Uncompressed_g1 = reshape(Uncompressed_g1,[256,256]);
J_g1 = uint8(Uncompressed_g1);

subplot(1,5,2);
imshow(J_g2);
title('Noisy Image');

subplot(1,5,3);
imshow(J_f_R);
title('Erased Image with Noise');

subplot(1,5,4);
imshow(J_g);
title('Reconstructed with Noise');

subplot(1,5,5);
imshow(J_g1);
title('Reconstructed with Noise Doubly Bridged');