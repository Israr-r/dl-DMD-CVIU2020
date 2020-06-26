
function im2 = im_norm(image)

    min_im = min(image(:));
    im1 = image-min_im;
    im2 = im1/(max(im1(:)));

end