clc;
clear all;
%%
plr = load('lm_lookuptable_b1.txt');
phr = load('lm_lookuptable_b2.txt');

B1 = 0 : 0.01 : 9.99;
B2 = 0 : 0.01 : 9.99;


surf(B1, B2, plr);
axis tight;
colormap(hot);
shading interp;
figure;
surf(B1, B2, phr);
axis tight;
colormap(hot);
shading interp;
