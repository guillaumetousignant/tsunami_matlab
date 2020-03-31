function [] = write_solution(filename, eta, amplitude, omega, points, elements, wall, farfield)
%READ_SU2 Summary of this function goes here
%   Detailed explanation goes here
fid = fopen(filename,'w');
if fid == -1
    error('read_su2:invalidInput', ['Unable to open file "', filename, '".']);
end

fprintf(fid, 'AMPLITUDE= %d\n', amplitude);
fprintf(fid, 'OMEGA= %d\n', omega);
fprintf(fid, 'NDIM= 3\n\n');

fprintf(fid, 'NPOIN= %d\n', size(points, 2));
for i = 1:size(points, 2)
    fprintf(fid, '%- 12.12g %- 12.12g %- 12.12g %- 12.12g %- 12.12g\n', points(1, i), points(2, i), points(3, i), real(eta(i, 1)), imag(eta(i, 1)));
end

fprintf(fid, '\nNELEM= %d\n', size(elements, 2));
for i = 1:size(elements, 2)
    fprintf(fid, '5 %d %d %d\n', elements(1, i), elements(2, i), elements(3, i));
end

fprintf(fid, 'NMARK= 2\n');
fprintf(fid, 'MARKER_TAG= wall\n');
fprintf(fid, 'MARKER_ELEMS= %d\n', size(wall, 2));
for i = 1:size(wall, 2)
    fprintf(fid, '3 %d %d\n', wall(1, i), wall(2, i));
end

fprintf(fid, 'MARKER_TAG= farfield\n');
fprintf(fid, 'MARKER_ELEMS= %d\n', size(farfield, 2));
for i = 1:size(farfield, 2)
    fprintf(fid, '3 %d %d\n', farfield(1, i), farfield(2, i));
end

fclose(fid);
end