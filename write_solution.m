function [] = write_solution(filename, eta, amplitude, omega)
%READ_SU2 Summary of this function goes here
%   Detailed explanation goes here
fid = fopen(filename,'w');
if fid == -1
    error('read_su2:invalidInput', ['Unable to open file "', filename, '".']);
end

fprintf(fid, 'AMPLITUDE= %d\n', amplitude);
fprintf(fid, 'OMEGA= %d\n\n', omega);

fprintf(fid, 'NPOIN= %d\n', size(eta, 1));
for i = 1:size(eta, 1)
    fprintf(fid, '%- 12.12g %- 12.12g\n', real(eta(i, 1)), imag(eta(i, 1)));
end

fclose(fid);
end