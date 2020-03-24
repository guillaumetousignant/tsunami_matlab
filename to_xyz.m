function xyz = to_xyz(sph)
    xyz = [sph(1)*sin(sph(2))*cos(sph(3)), sph(1)*sin(sph(2))*sin(sph(3)), sph(1)*cos(sph(2))];
end