function [] = tsunami(filename)
%TSUNAMI Summary of this function goes here
%   Detailed explanation goes here
[points, elements, wall, farfield] = read_su2(filename);

figure()
trimesh(elements', points(1, :)', points(2, :)', -points(3, :)');

end

