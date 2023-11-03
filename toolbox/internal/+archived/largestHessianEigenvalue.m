function lambda = largestHessianEigenvalue(Dxx, Dxy, Dyy)
    %   largestHessianEigenvalue Calculates the eigenvalues of the hessian matrix and returns the largest one
    %
    %   lambda = largestHessianEigenvalue(Dxx, Dxy, Dyy)
    %
    %   Inputs:
    %       Dxx, Dxy, Dyy: Second order derivatives of the image obtained with `hessian2D`
    %
    %   Outputs:
    %       lambda: Largest eigenvalue of the hessian matrix
    %
    %   Note:
    %       Function adapted from https://github.com/timjerman/JermanEnhancementFilter
    %       first by Zhengda Li and then by Franco Tavella
    arguments
        Dxx (:,:) {mustBeNumeric}
        Dxy (:,:) {mustBeNumeric}
        Dyy (:,:) {mustBeNumeric}
    end
    % Compute the eigenvalues of J
    tmp = sqrt((Dxx - Dyy).^2 + 4*Dxy.^2);
    mu1 = 0.5*(Dxx + Dyy + tmp);
    mu2 = 0.5*(Dxx + Dyy - tmp);
    % Get largest eigenvalue
    check = abs(mu1) > abs(mu2);
    lambda = mu2; 
    lambda(check) = mu1(check);
end