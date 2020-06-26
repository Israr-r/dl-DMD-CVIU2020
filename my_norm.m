function imax = my_norm(im)


imin = im - min(im(:));
imax = imin/max(imin(:));

