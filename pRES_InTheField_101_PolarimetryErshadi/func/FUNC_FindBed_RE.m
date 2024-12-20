% Written by Reza Ershadi
% May 2023
% University of Tübingen - Geophysics group
function [zBed,iBed] = FUNC_FindBed_RE(s,z)
    pwr = 20.*log10(abs(s));
    dz = mean(diff(z));
    sdata = smoothdata(pwr, 'movmedian', 100);
    d1 = abs(gradient(sdata, z));

% figure,
% subplot(1,2,1)
% plot(pwr,z)
% set(gca,'YDir','reverse')
% subplot(1,2,2)
% plot(d1,z)
% set(gca,'YDir','reverse')

    try
        [~,maxd1] = max(d1);
        rng = round(50/dz,0);
        i_bedrange = [maxd1-rng maxd1+rng];
        z_bedrange = [z(i_bedrange(1)) z(i_bedrange(2))];
        [~,iBed] = max(pwr(i_bedrange(1):i_bedrange(2)));
        iBed = iBed + i_bedrange(1);
        zBed = z(iBed);
    catch
        zBed = max(z);
        iBed = length(z);
    end
end