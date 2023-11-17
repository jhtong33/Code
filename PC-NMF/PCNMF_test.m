% Demo of using PC-NMF to separate biological chorus and the other noise sources from long-term spectrograms
% Written by Tzu-Hao Harry Lin (Academia Sinica, Taiwan), 2017/5/11
% Corresponding email: schonkopf@gmail.com

% load data
clear; clc; load('data.mat')

% Assign the number of sound source in "k_source"
% If the number of source is unknow, please specify a range (e.g., k_source=2:5;)
% The clustering algorithm will estimate the number fo cluster according to clustering consistency
k_source=2; 
% assign the number of time frames to cascade input spectrogram, 1 represents no cascade (original setting in this study). 
% In most situations, cascade spectrogram by 3 or 5 frames will result in better separation performance
cascading_length=1;

%%% Running PC-NMF on simulation data
%
% Please note that PC-NMF has its uncertainty, different trials may produce
% slightly different results
%
[output, W, W_cluster, H]=LTSA_PCNMF(matrix_standardization(marine_noisy), k_source, cascading_length);
figure(1);
imagesc(matrix_standardization(data)); axis xy; colormap('jet'); caxis([0 1])
figure(2); 
for n=size(output,3):-1:1
    subplot(size(output,3),1,n); imagesc(output(:,:,n)); axis xy; colormap('jet'); caxis([0 0.5]);
    ylabel('Frequency [Hz]'); set(gca,'YScale','log')

end
title('Separation of periodical and non-periodical signals on simulation data')

% %%% Running PC-NMF on real-world data
% [output, W, W_cluster]=LTSA_PCNMF(matrix_standardization(marine_noisy), k_source, cascading_length);
% figure(2); 
% for n=size(output,3):-1:1
%     subplot(size(output,3),1,n); imagesc(output(:,:,n)); axis xy; colormap('jet'); caxis([0 1])
% end
% title('Separation of biological chorus in noisy marine environment')
% 
% [output, W, W_cluster]=LTSA_PCNMF(matrix_standardization(marine_quiet), k_source, cascading_length);
% figure(3); 
% for n=size(output,3):-1:1
%     subplot(size(output,3),1,n); imagesc(output(:,:,n)); axis xy; colormap('jet'); caxis([0 1])
% end
% title('Separation of biological chorus in quiet marine environment')
% 
% [output, W, W_cluster]=LTSA_PCNMF(matrix_standardization(terrestrial_noisy), k_source, cascading_length);
% figure(4); 
% for n=size(output,3):-1:1
%     subplot(size(output,3),1,n); imagesc(output(:,:,n)); axis xy; colormap('jet'); caxis([0 1])
% end
% title('Separation of biological chorus in noisy terrestrial environment')
% 
% [output, W, W_cluster]=LTSA_PCNMF(matrix_standardization(terrestrial_quiet), k_source, cascading_length);
% figure(5); 
% for n=size(output,3):-1:1
%     subplot(size(output,3),1,n); imagesc(output(:,:,n)); axis xy; colormap('jet'); caxis([0 1])
% end
% title('Separation of biological chorus in quiet terrestrial environment')
