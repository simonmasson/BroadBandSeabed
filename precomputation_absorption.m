function gamma = precomputation_absorption(data, depth_estimation, ff)
    SoundVelocity = data.env.SoundSpeed ;
    pH = data.env.Acidity;
    Sal = data.env.Salinity;
    Temp = data.env.Temperature;

    z=depth_estimation;

    A1=(8.86/SoundVelocity)*10^(0.78*pH-5);
    f1=2.8*((Sal/35)^0.5)*10^(4-1245/(Temp+273));
    P1=1;
    A2=21.44*(Sal/SoundVelocity)*(1+0.025*Temp);
    P2=1-(1.37*10^-4*z)+(6.2*(10^-9)*z^2);
    f2=(8.17*10^(8-(1990)/(Temp+273)))/(1+0.0018*(Sal-35));
    P3=1-(3.83*10^-5*z)+(4.9*10^-10*z^2);
    if Temp<=20
        A3=(4.937e-4)-(2.59*10^-5*Temp)+(9.11*10^-7*Temp^2)-(1.5*10^-8*Temp^3);
    else
        A3=(3.964e-4)-(1.146e-5*Temp)+(1.45e-7*Temp^2)-(6.5e-10*Temp^3);
    end

    ffk = ff/1000;
    alpha = ((A1*P1*f1*ffk.^2)./(ffk.^2+f1^2))+((A2*P2*f2*ffk.^2)./(ffk.^2+f2^2))+(A3*P3*ffk.^2);
    gamma = alpha / 1000 / 20 / log10(exp(1));
end