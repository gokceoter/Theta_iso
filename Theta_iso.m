% PDEW2013  1  5  8 58 19.30  55.3900 -134.6500  10.0 6.4 7.7  
% SOUTHEASTERN ALASKA
% event name:    201301050858A

Mrr=7;
Mtt=2;
Mpp=3;
Mrt=3;
Mrp=5;
Mtp=4;

% build the given full moment tensor
M=[Mrr Mrt Mrp; Mrt Mtt Mtp; Mrp Mtp Mpp];
% give two input parameters M0 of earthquake and C_ratio
M0=2.47e+27;
lambda=1; mu=0.7;
C_eig_ratio=3/2*lambda/mu + 1;

% calculate u.n or cos(theta) from the isotropic part of the moment tensor
costheta=trace(M)/(2*M0*C_eig_ratio);
theta=acosd(costheta); %the angle between u and n in degrees

% calculate the deviatoric part of the moment tensor from the full tensor
Mdev=M - (trace(M)/3)*eye(3,3);
[eigvec_Mdev,eigval_Mdev]=eig(Mdev);
Miso=M-Mdev;
norm(Miso);
D=diag(Mdev);
norm(Miso)/(norm(Miso)+maxabs(D))