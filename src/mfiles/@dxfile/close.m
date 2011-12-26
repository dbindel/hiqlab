% close(dxf)

function close(dxf)

fprintf(dxf.fid1, 'end\n');
fclose(dxf.fid1);
fclose(dxf.fid2);

