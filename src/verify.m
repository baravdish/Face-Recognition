%%%%%%%%%%%%%%%%%%%%%%%%%%
% input: Image to be verified.
%
% output: The identity number (integer) of the identified person,
% i.e. ‘1’, ‘2’,...,‘16’ for the persons belonging to ‘db1’
% and ‘0’ for all other faces.
%%%%%%%%%%%%%%%%%%%%%%%%%%

function output = verify(input)
  database = load('database.mat');
  sizeOfDatabase = size(database.database);
  [id, id_false, min_value] = testImage(input, database.database, sizeOfDatabase(1));
  
  output = id;
end
