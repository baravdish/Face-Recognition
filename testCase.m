function result = testCase(img)

sigma = 0.6;
light = 0.8;
scale = 1.1;
rotation = 5;

imgBlur = applyLP(img, sigma);
result.blurFace = detectFace(imgBlur);
result.blurResult = tnm034(imgBlur);

imgTone = applyTone(img, light);
result.toneFace = detectFace(imgTone);
result.toneResult

imgScale = imresize(img, scale);
result.scaleFace = detectFace(imgTone);
result.scaleResult

imgRot = imrotate(img, rotation, 'crop');
result.rotationFace = detectFace(imgRot);
result.rotationResult


