function [a_point, b_point, a_sim, b_sim, C_orig_jitCorr, C_surr_jitCorr] = computeBandsV2(C_orig,C_surr,lags)
    
    % Computes CCG jitter acceptance based on Amarasingham 2012
    % Similar to Zach's acceptband.m script but only computes the bands 
    % Written by Atanas Stankov. 
    n_surr = size(C_surr,2);

    % compute surrogate mean
    C_mean_surr = mean(C_surr,2);
    
    % pool original CCG with surrogates
    C_pool = [C_surr C_orig];
    
    % sort in ascending order and retrieve pointwise acceptance bands
    M = size(C_pool,2);
    C_pool_sorted = sort(C_pool,2); 
    
    ind_pt025 = round(0.025*M);
    ind_pt975 = round(0.975*M);
    
    % pointwise bands
    a_point = C_pool_sorted(:,ind_pt025);
    b_point = C_pool_sorted(:,ind_pt975);
    
    % compute pool mean and std
    C_mean_pool = (1/(M-1)).*sum(C_pool_sorted(:,1:M-1),2);
    C_std_pool  = sqrt((1./(M-2)).*sum((C_pool_sorted(:,1:M-1) - C_mean_pool).^2, 2));
    
    % normalize pool 
    C_pool_star = (C_pool - C_mean_pool)./C_std_pool;
    
    % sort in ascending order and retrieve simultaneous acceptance bands
    C_plus  = max(C_pool_star);
    C_minus = min(C_pool_star);
    C_plus  = sort(C_plus);
    C_minus = sort(C_minus);

    a_sim = C_minus(ind_pt025)*C_std_pool + C_mean_pool;
    b_sim = C_plus(ind_pt975)*C_std_pool + C_mean_pool;
    
    % jitter demean
    C_orig_jitCorr = C_orig - C_mean_surr;
    C_surr_jitCorr = C_surr - C_mean_surr;
    
%     a_point = a_point - C_mean_surr;
%     b_point = b_point - C_mean_surr;
%     
%     a_sim = a_sim - C_mean_surr;
%     b_sim = b_sim - C_mean_surr;

    %%

    % plot    
%     patch([lags; lags(end:-1:1); lags(1)],[a_sim',b_sim(end:-1:1)',b_sim(1)'],.6.*ones(1,3),'FaceAlpha',1)
%     hold on
%     patch([lags; lags(end:-1:1); lags(1)],[a_point',b_point(end:-1:1)',b_point(1)'],.9.*ones(1,3),'FaceAlpha',1)
%     plot(lags,C_orig_jitCorr,'b','lineWidth',2);
%     hold off
    
end