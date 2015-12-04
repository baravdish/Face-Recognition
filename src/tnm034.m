%%%%%%%%%%%%%%%%%%%%%%%%%%
% im: Image of unknown face, RGB-image in uint8 format in the
% range [0,255]
%
% id: The identity number (integer) of the identified person,
% i.e. ‘1’, ‘2’,...,‘16’ for the persons belonging to ‘db1’
% and ‘0’ for all other faces.
%%%%%%%%%%%%%%%%%%%%%%%%%%

function [id, id_false, min_value] = tnm034(im, height, width, threshold, kernel_size, decorr, freqestim)
  balanced_image = colorCorrection(im);
  face_image = detectFace(balanced_image);
  [id, id_false, min_value] = verify(face_image, height, width, threshold, kernel_size, decorr, freqestim);
end
