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
