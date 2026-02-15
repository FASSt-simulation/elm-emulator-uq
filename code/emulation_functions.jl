# Import libraries
using Turing
using StatsPlots
using LinearAlgebra
using Distributions
using MultivariateStats
import MultivariateStats: reconstruct
using GaussianProcesses
using ScikitLearn
using StatsBase
using Suppressor
import Base.rand
import Distributions: logpdf

ProjDir = @__DIR__
cd(ProjDir)

# ========== Utilities ==========
V2M(v) = reshape(v, length(v), 1)

# ========== Principle Component Analysis Methods ==========

project(pca_fit, Y) = projection(pca_fit)' * centralize(Y, pca_fit.mean) # [K x M] matrix; projections of the N observations onto leading K components, where Y is [D x M] array (set Y = X to project original data)

project(pca_fit, Y, J) = (projection(pca_fit)[:,1:J])' * centralize(Y, pca_fit.mean)

#reconstruct_J(pca_fit, W, J) = decentralize(projection(pca_fit)[:,1:J] * W, pca_fit.mean) # as above, but only using leading J <= K components to reconstruct, where W is [J x M]
function reconstruct_J(pca_fit, Z_hat, J)
    T_PCA = projection(pca_fit)[:, 1:J]
    W = Z_hat

    y_hat = decentralize(T_PCA * W, pca_fit.mean)
    return y_hat
end

function calculate_10day_average(data0)
    num_years = 5  
    monthday=[10, 10,11,10,10,8,10,10,11,10,10,10,10,10,11,10,10,10,10,10,11,10,10,11,10,10,10,10,10,11,10,10,10,10,10,11];
    monthly_average = Vector{}()
    end_day=0
    for year in 1:num_years
        start_day=end_day+1
        for month in 1:36
            end_day = start_day + monthday[month]-1
            temp=mean(filter(!isnan,data0[start_day:end_day,:]), dims=1)
            push!(monthly_average,temp[1])
            start_day=end_day+1
        end
    end
    return monthly_average
    
end

# validating accuracy: cumulative proportion of variance explained
cumpropvar(pca_fit) = vec(cumsum(eigvals(pca_fit) ./ sum(eigvals(pca_fit))))

rms_error(pca_fit, Y) = [rmsd(Y, reconstruct_J(pca_fit, project(pca_fit,Y,J), J)) for J = 1:size(projection(pca_fit), 2)]; # compute RMSE for all J <= K

# ========== Emulation Methods ==========
function predict_y_GP(GP, p)
    # Learn a GP model for each PCA component.
    crossdata = GaussianProcesses.KernelData(GP.kernel, GP.x, p)
    Kcross = cov(GP.kernel, GP.x, p, crossdata)
    mx = mean(GP.mean, p)

    # Return means of the regression curve without the covariance matrix for reducing computational complexity.
    return mx + Kcross' * GP.alpha
end

function emulator(p, d, C = 0)
    p_norm = p'         # parameters are already normalized
    N_p = (C == 0) ? size(d.T_PCA, 2) : C           # Determine the number of PCA components for reconstruction.
    Z_std_hat = hcat([predict_y_GP(d.GPs[j], p_norm) for j=1:N_p ]...)';
    Z_hat = Z_std_hat.*d.σ_z[1:N_p] .+ d.μ_z[1:N_p]
    y_hat = reconstruct_J(d.T_PCA, Z_hat, N_p) * 60 * 60 * 24 # convert to units of days
end

function calibration_emulator(p, d, C = 0)
    p_norm = p         # parameters are already normalized
    N_p = (C == 0) ? size(d.T_PCA, 2) : C           # Determine the number of PCA components for reconstruction.
    Z_std_hat_vec = [predict_y(d.GPs[j], p_norm')[1] for j = 1:N_p]
    Z_std_hat = hcat(Z_std_hat_vec...)'
    Z_hat = Z_std_hat .* d.σ_z[1:N_p] .+ d.μ_z[1:N_p]
    y_hat = reconstruct_J(d.T_PCA, Z_hat, N_p) * 60 * 60 * 24 # convert to units of days

end

function calibration_emulator1(p, d, C = 0)
    p_norm = p         # parameters are already normalized
    N_p = (C == 0) ? size(d.T_PCA, 2) : C           # Determine the number of PCA components for reconstruction.
    Z_std_hat_vec = [predict_y(d.GPs[j], p_norm')[1] for j = 1:N_p]
    Z_std_hat = hcat(Z_std_hat_vec...)'
    Z_hat = Z_std_hat .* d.σ_z[1:N_p] .+ d.μ_z[1:N_p]
    y_hat = reconstruct_J(d.T_PCA, Z_hat, N_p) * 60 * 60 * 24 # convert to units of days

    num_years = 5 
    monthday=[10, 10,11,10,10,8,10,10,11,10,10,10,10,10,11,10,10,10,10,10,11,10,10,11,10,10,10,10,10,11,10,10,10,10,10,11];
    monthly_average = Vector{}()
    end_day=0
    for year in 1:num_years
        start_day=end_day+1
        for month in 1:36
            end_day = start_day + monthday[month]-1
            temp=mean(filter(!isnan,y_hat[start_day:end_day,:]), dims=1)
            push!(monthly_average,temp[1])
            start_day=end_day+1
        end
    end
    return monthly_average
end
