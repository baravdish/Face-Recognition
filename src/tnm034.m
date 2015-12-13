%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: Image of unknown face, RGB-image in uint8 format in the
% range [0,255]
%
% id: The identity number (integer) of the identified person,
% i.e. ‘1’, ‘2’,...,‘16’ for the persons belonging to ‘db1’
% and ‘0’ for all other faces.
%%%%%%%%%%%%%%%%%%%%%%%%%%

function id = tnm034(im)
  try
    balanced_image = colorCorrection(input);
    face_image = detectFace(balanced_image);
    
    id = verify(face_image);
  catch exception
    id = -1;
  end
end
