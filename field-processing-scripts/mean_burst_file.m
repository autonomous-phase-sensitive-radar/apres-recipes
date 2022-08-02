function vdat = mean_burst_file(filename,AttNum)
% Function that takes the mean of all the bursts found inside a .DAT file
% Need to specify which attenuation setting to take average of

% George Lu, July 2022

max_bursts = 100;
% get initial burst
burst_count = 1;
vdat = Field_load(filename,1);
vdats = Field_burst_split_by_att(vdat);
vdat = Field_burst_mean(vdats(AttNum));

% Now iterate
for k = 2:max_bursts % arbitrary max of 100
    vdat_temp = Field_load(filename,k);
    if vdat_temp.Code == -4 % burst not found in file
        break;
    elseif vdat_temp.Code == -5
        break;
    else
        vdats = Field_burst_split_by_att(vdat_temp);
        burst_count = burst_count + 1;
        vdat_temp = Field_burst_mean(vdats(AttNum));  
        vdat.vif = vdat.vif+vdat_temp.vif;
        vdat.chirpTime = vdat.chirpTime+vdat_temp.chirpTime;
        vdat.chirpAtt = vdat.chirpAtt+vdat_temp.chirpAtt;
        
    end
    
end

vdat.ChirpsInBurst = 1;
vdat.processing = [vdat.processing {'burst mean'}];
vdat.vif = vdat.vif./burst_count;
vdat.chirpTime = vdat.chirpTime./burst_count;
vdat.chirpAtt = vdat.chirpAtt./burst_count;

end