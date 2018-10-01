function GradL = GradientL(z1)
global z1Star

GradL = 2*0.5*(0.5*z1 - z1Star);

%GradL = 2*0.5*z1 - 2*z1Star;
end