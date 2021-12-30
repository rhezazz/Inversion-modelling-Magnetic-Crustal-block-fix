%Objective function for magnetic data inversion
function [misfit] = misfit_magnetic(dBz_obs, dBz_cal)
    ls = length(dBz_obs);
    for j = 1 : ls
        m(j) = ((dBz_cal(j) - (dBz_obs(j))))^2;
    end
    misfit = sqrt((1/ls)*sum(m));
end
