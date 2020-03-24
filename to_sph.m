function sph = to_sph(xyz)
    r = norm(xyz);
    sph = [r, acos(xyz(3)/r), atan2(xyz(2), xyz(1))]; % Used atan(y/x) before, this is the real thing
    sph(isnan(sph)) = 0; % not needed, atan2 doesn't return NaN like atan % turns out acos still does
end