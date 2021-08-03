data{
    int<lower=0>L; //rows in dataset
    int<lower=0>S; //num of studies
    int<lower=0>M; //number of patients
    int AGold[L];
    int AGmid[L];
    int NotMild[L];
    //int SevYN[L];
    //int ModYN[L];
    int Male[L];
    //int Estimated[L];
    vector[L] LOD;
    vector[L] logvalue;
    int Day[L];
    int PatID[L];
    int StudyNum[L];
}
parameters{
    vector[S] b1;
    vector[S] a1;
    vector[M] b2;
    vector[M] a2;
    real a0;
    real b0;
    real aGM;
    real bGM;
    //real aA1;
    //real bA1;
    //real aA2;
    //real bA2;
    real aSnM;
    real bSnM;
    //real aSS;
    //real bSS;
    vector<lower=0>[2] sigma_study;
    vector<lower=0>[2] sigma_PatID;
    corr_matrix[2] Rho_study;
    corr_matrix[2] Rho_PatID;
    real<lower=0>sigma;
}
model{
    vector[L] mu;
    vector[L] A;
    vector[L] B;
    sigma ~ cauchy( 0 , 2 );
    Rho_PatID ~ lkj_corr( 2 );
    Rho_study ~ lkj_corr( 2 );
    sigma_PatID ~ cauchy( 0 , 2 );
    sigma_study ~ cauchy( 0 , 2 );
    b0 ~ normal( 0 , 2 );
    a0 ~ normal( 5 , 5 );
    {
    vector[2] YY[M];
    vector[2] MU;
    MU = [ 0 , 0 ]';
    for ( j in 1:M ) YY[j] = [ a2[j] , b2[j] ]';
    YY ~ multi_normal( MU , quad_form_diag(Rho_PatID , sigma_PatID) );
    }
    {
    vector[2] YY[S];
    vector[2] MU;
    MU = [ 0 , 0 ]';
    for ( j in 1:S ) YY[j] = [ a1[j] , b1[j] ]';
    YY ~ multi_normal( MU , quad_form_diag(Rho_study , sigma_study) );
    }
    for ( i in 1:L ) {
        B[i] = b0 + b1[StudyNum[i]] + b2[PatID[i]] + bSnM*NotMild[i] + bGM*Male[i];
    }
    for ( i in 1:L ) {
        A[i] = a0 + a1[StudyNum[i]] + a2[PatID[i]] + aSnM*NotMild[i] + aGM*Male[i];
    }
    for ( i in 1:L ) {
        mu[i] = A[i] + B[i] * Day[i];
    }
    for ( i in 1:L ) 
        if ( logvalue[i] <= 0.001 ) target += normal_lcdf(LOD[i] | mu[i], sigma);
    for ( i in 1:L ) 
        if ( logvalue[i] > 0.001 ) target += normal_lpdf(logvalue[i] | mu[i], sigma);
}
generated quantities{
    vector[L] log_lik;
    vector[L] mu;
    vector[L] A;
    vector[L] B;
    for ( i in 1:L ) {
        B[i] = b0 + b1[StudyNum[i]] + b2[PatID[i]] + bSnM*NotMild[i] + bGM*Male[i];
    }
    for ( i in 1:L ) {
        A[i] = a0 + a1[StudyNum[i]] + a2[PatID[i]] + aSnM*NotMild[i] + aGM*Male[i];
    }
    for ( i in 1:L ) {
        mu[i] = A[i] + B[i] * Day[i];
    }
    for ( i in 1:L ) {
        if ( logvalue[i] <= 0.001 ) log_lik[i] = normal_lcdf(LOD[i] | mu[i], sigma);
    }
    for ( i in 1:L ){ 
        if ( logvalue[i] > 0.001 ) log_lik[i] = normal_lpdf(logvalue[i] | mu[i], sigma);
    }
}
