function plyWrite(pts, rgb, filename)
% Written by Chenxi cxliu@ucla.edu
% Input: fname: output file name, e.g. 'data.ply'
%        P: 3*m matrix with the rows indicating X, Y, Z
%        C: 3*m matrix with the rows indicating R, G, B

num = size(pts, 2);
header = 'ply\n';
header = [header, 'format ascii 1.0\n'];
header = [header, 'comment written by Chenxi\n'];
header = [header, 'element vertex ', num2str(num), '\n'];
header = [header, 'property float32 x\n'];
header = [header, 'property float32 y\n'];
header = [header, 'property float32 z\n'];
header = [header, 'property uchar red\n'];
header = [header, 'property uchar green\n'];
header = [header, 'property uchar blue\n'];
header = [header, 'end_header\n'];

data = [pts', double(rgb')];

fid = fopen(filename, 'w');
fprintf(fid, header);
dlmwrite(filename, data, '-append', 'delimiter', '\t', 'precision', 3);
fclose(fid);


end