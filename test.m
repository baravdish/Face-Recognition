function tests = test
  tests = functiontests(localfunctions);
end

function setupOnce(testCase)  % do not change function name
  addpath src

  [access_images, number_of_access_images] = readAllFromDir('access', 'img/access/', '*.jpg');
  [no_access_images, number_of_no_access_images] = readAllFromDir('no_access', 'img/no_access/', '*.jpg');
  [hard_images, number_of_hard_images] = readAllFromDir('hard', 'img/hard/', '*.jpg');

  testCase.TestData.access_images = access_images;
  testCase.TestData.number_of_access_images = number_of_access_images;
  testCase.TestData.no_access_images = no_access_images;
  testCase.TestData.number_of_no_access_images = number_of_no_access_images;
  testCase.TestData.hard_images = hard_images;
  testCase.TestData.number_of_hard_images = number_of_hard_images;

end

function teardownOnce(testCase)  % do not change function name
  rmpath src
end

function testReadAccessImages(testCase)
  verifyEqual(testCase,testCase.TestData.number_of_access_images, 16)
end

function testReadNoAccessImages(testCase)
  verifyEqual(testCase,testCase.TestData.number_of_no_access_images, 4)
end

function testReadHardImages(testCase)
  verifyEqual(testCase,testCase.TestData.number_of_hard_images, 38)
end
