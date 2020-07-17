using MAT

W = [    0    0.3314    0.3314    0.1657    0.1482    0.1482;
    0.3314         0    0.1657    0.1482    0.1657    0.1172;
    0.3314    0.1657         0    0.1482    0.1172    0.1657;
    0.1657    0.1482    0.1482         0    0.3314    0.3314;
    0.1482    0.1657    0.1172    0.3314         0    0.1657;
    0.1482    0.1172    0.1657    0.3314    0.1657         0]

function load_sample(mat_file_name::String, var_name::String = "spat_data")
    file = matopen(mat_file_name)
    var = read(file, var_name) # note that this does NOT introduce a variable ``varname`` into scope
    close(file)
    # If No var, then broken
    keywd = split(mat_file_name, "/")[end]
    start_idx = 3
    if keywd == "gaussian_diff"
        start_idx = 4
    end

    SPtrend = 3.6814e-04
    spat_data = var[:, start_idx:end]' .- SPtrend
end


#
#     %% Recursive Prediction
#     spData= zeros(r,N-1);
#
#     sample_ident = 1:6;
#     raw_data = spat_data_detrended(:, 1 : order);
#     raw_label = spat_data_detrended(:, order + 1);
#     vx = W * raw_data;
#     vy = raw_label(sample_ident,:);
#
#     filter_wrapper = KRLS(vx, vy, kernel_mode);
#
# %     some_other_data = W * spat_data_detrended;
#
#     for k=2:N-1
#         vx = W * spat_data_detrended(:, k - order + 1 : k);
#         vy = spat_data_detrended(sample_ident, k + 1);
#         filter_wrapper.update(vx, vy);
#         vx_next_step = W * spat_data_detrended(:, k + 1);
#         spData(sample_ident, k) = filter_wrapper.predict(vx_next_step);
#     end
#
#     %% Plots
#     figure;
#     for ikk=1:6
#        subplot(3,2,ikk);
#        plot(1:N,spat_data(1:N,ikk));
#        hold on;
#        plot(2:N,spData(ikk,:)+SPtrend,'r');
#        hold on;
#        set(gca, 'Fontsize', 24);
#        set(get(gca,'Children'),'linewidth',2.0);
#        title(['S',num2str(ikk)]);
#        grid on;
#     end
#
#     h= suptitle('Prediction of Spatial Component');
#     set(h, 'FontSize',32)
# end
