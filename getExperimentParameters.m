clc;
clear all;
%%
sd1 = readbin('single_diff1.bin'); 
sd2 = readbin('single_diff1.bin'); 

lbd1 = readbin('lookuptable_basic_diff1.bin'); 
lbd2 = readbin('lookuptable_basic_diff2.bin'); 

lld1 = readbin('lookuptable_lm_diff1.bin'); 
lld2 = readbin('lookuptable_lm_diff2.bin'); 

ld1 = readbin('lm_diff1.bin'); 
ld2 = readbin('lm_diff2.bin'); 

cld1 = readbin('combined_lookuptable_lm_diff1.bin'); 
cld2 = readbin('combined_lookuptable_lm_diff2.bin'); 

ld_2000_1 = readbin('2000_lm_diff1.bin'); 
ld_2000_2 = readbin('2000_lm_diff2.bin');

cbld1 = readbin('combined_basic_lm_diff1.bin'); 
cbld2 = readbin('combined_basic_lm_diff2.bin');

fprintf('sd1 mean: %f, max: %f, min: %f, std: %f, psnr: %f\n', mean2(sd1), max(max(sd1)), min(min(sd1)), std2(sd1), getPSNR(sd1));
fprintf('sd2 mean: %f, max: %f, min: %f, std: %f, psnr: %f\n', mean2(sd2), max(max(sd2)),min(min(sd2)), std2(sd2), getPSNR(sd2));
fprintf('lbd1 mean: %f, max: %f, min: %f, std: %f, psnr: %f\n', mean2(lbd1), max(max(lbd1)), min(min(lbd1)), std2(lbd1), getPSNR(lbd1));
fprintf('lbd2 mean: %f, max: %f, min: %f, std: %f, psnr: %f\n', mean2(lbd2), max(max(lbd2)),min(min(lbd2)), std2(lbd2), getPSNR(lbd2));
fprintf('lld1 mean: %f, max: %f, min: %f, std: %f, psnr: %f\n', mean2(lld1), max(max(lld1)),min(min(lld1)), std2(lld1), getPSNR(lld1));
fprintf('lld2 mean: %f, max: %f, min: %f, std: %f, psnr: %f\n', mean2(lld2), max(max(lld2)),min(min(lld2)), std2(lld2), getPSNR(lld2));
fprintf('ld1 mean: %f, max: %f, min: %f, std: %f, psnr: %f\n', mean2(ld1), max(max(ld1)),min(min(ld1)), std2(ld1), getPSNR(ld1));
fprintf('ld2 mean: %f, max: %f, min: %f, std: %f, psnr: %f\n', mean2(ld2), max(max(ld2)),min(min(ld2)), std2(ld1), getPSNR(ld2));
fprintf('cld1 mean: %f, max: %f, min: %f, std: %f, psnr: %f\n', mean2(cld1), max(max(cld1)),min(min(cld1)), std2(cld1), getPSNR(cld1));
fprintf('cld2 mean: %f, max: %f, min: %f, std: %f, psnr: %f\n', mean2(cld2), max(max(cld2)),min(min(cld2)), std2(cld1), getPSNR(cld2));
fprintf('ld_2000_1 mean: %f, max: %f, min: %f, std: %f, psnr: %f\n', mean2(ld_2000_1), max(max(ld_2000_1)),min(min(ld_2000_1)), std2(ld_2000_1), getPSNR(ld_2000_1));
fprintf('ld_2000_2 mean: %f, max: %f, min: %f, std: %f, psnr: %f\n', mean2(ld_2000_2), max(max(ld_2000_2)),min(min(ld_2000_2)), std2(ld_2000_2), getPSNR(ld_2000_2));
fprintf('cbld1 mean: %f, max: %f, min: %f, std: %f, psnr: %f\n', mean2(cbld1), max(max(cbld1)),min(min(cbld1)), std2(cbld1), getPSNR(cbld1));
fprintf('cbld2 mean: %f, max: %f, min: %f, std: %f, psnr: %f\n', mean2(cbld2), max(max(cbld2)),min(min(cbld2)), std2(cbld2), getPSNR(cbld2));



%%
sd1 = reshape(sd1, 1, 400 * 256);
figure;
plot(sd1);
sd2 = reshape(sd2, 1, 400 * 256);
figure;
plot(sd2);
lbd1 = reshape(lbd1, 1, 400 * 256);
figure;
plot(lbd1);
lbd2 = reshape(lbd2, 1, 400 * 256);
figure;
plot(lbd2);
lld1 = reshape(lld1, 1, 400 * 256);
figure;
plot(lld1);
lld2 = reshape(lld2, 1, 400 * 256);
figure;
plot(lld2);
ld1 = reshape(ld1, 1, 400 * 256);
figure;
plot(ld1);
ld2 = reshape(ld2, 1, 400 * 256);
figure;
plot(ld2);
cld1 = reshape(cld1, 1, 400 * 256);
figure;
plot(cld1);
cld2 = reshape(cld2, 1, 400 * 256);
figure;
plot(cld2);
ld_2000_1 = reshape(ld_2000_1, 1, 400 * 256);
figure;
plot(ld_2000_1);
ld_2000_2 = reshape(ld_2000_2, 1, 400 * 256);
figure;
plot(ld_2000_2);
cbld1 = reshape(cbld1, 1, 400 * 256);
figure;
plot(cbld1);
cbld2 = reshape(cbld2, 1, 400 * 256);
figure;
plot(cbld2);

%%
%get psnr of reconstruction image
lm_fbp1 = readbin('lm_FBP_1.bin');
lm_fbp2 = readbin('lm_FBP_2.bin');
lm_ml1 = readbin('lm_ml_em_img1.bin');
lm_ml2 = readbin('lm_ml_em_img2.bin');
fprintf('lm_fbp1 psnr: %f\n', getPSNR(lm_fbp1));
fprintf('lm_fbp2 psnr: %f\n', getPSNR(lm_fbp2));
fprintf('lm_ml1 psnr: %f\n', getPSNR(lm_ml1));
fprintf('lm_ml2 psnr: %f\n', getPSNR(lm_ml2));

