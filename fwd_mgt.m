%Reference : Lowrie, W. 2007. Fundamental of Geophysics. Cambridge University Press, p. 330-332.
function [d_Bz] = fwd_mgt(x,x0,z1,z2,mu0,delta_Mz,m)
conts = (mu0)./(2*pi);
for i = 1 : length(x)
    alpha1(i) = atan(((x(i)-x0)+m)./z1);
    alpha2(i) = atan(((x(i)-x0)-m)./z1);
    alpha3(i) = atan(((x(i)-x0)+m)./z2);
    alpha4(i) = atan(((x(i)-x0)-m)./z2);
    d_Bz(i) = conts.*delta_Mz.*((alpha1(i)-alpha2(i))-(alpha3(i)-alpha4(i)));
end
end