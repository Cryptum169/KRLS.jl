using Revise

import KRLS

using Plots

filter_param = KRLS.FilterParam(4, 10.0, 6)

spat_data = KRLS.load_sample("Data/Gaussian diff/Gaussian_diff_spat.mat")
(r, c) = size(spat_data)
raw_data = spat_data[:, 1]
vx = KRLS.W * raw_data
vy = spat_data[:, 2]
predicted_data = zeros((r, c))
predicted_data[:, 1] = raw_data

krlsfilter = KRLS.KRLSFilter(vx, vy, '3', filter_param)

println("Something")

vx_arr = KRLS.W * spat_data

@time for k = 2:c - 1
    KRLS.update(vx_arr[:, k], spat_data[:, k + 1], krlsfilter)
    predicted_data[:, k] = KRLS.predict(vx_arr[:, k + 1], krlsfilter)
end

plot(spat_data[4,:])

plot!(predicted_data[4, :])
