%%
function L = watershedTransform(bw)
% Compute the distance transform of the complement of the binary image.
D = bwdist(~bw);
figure, imshow(D,[],'InitialMagnification','fit'), title('Distance transform of ~bw');
%%
% Complement the distance transform, and force pixels that don't belong to the objects to be at Inf .
D = -D;
D(~bw) = Inf;
%%
%Compute the watershed transform and display the resulting label matrix as an RGB image.
L = watershed(D);
L(~bw) = 0;
figure, imshow(L), title('Watershed transform of D');
rgb = label2rgb(L,'jet','w');
figure, imshow(rgb,'InitialMagnification','fit'), title('Watershed transform of D');
end