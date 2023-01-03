function xt = heaviside_threshold(x, beta, eta)
%THRESHOLd applies the tangent-hyperbolicus approximation of the heaviside
%function to x with parameters beta and eta
xt = (tanh(beta*eta)+tanh(beta*(x-eta)))/...
    (tanh(beta*eta)+tanh(beta*(1-eta)));
end