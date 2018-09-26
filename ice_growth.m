function output = ice_growth(sal,time,time_frame)
%calculates ice growth from salinity
%   sal(depth,time) in psu (2D array)
%   time in days (array)
%   time_frame in days
%   h0 height of watercolumn. usually 100 [m] 
%   ice growth in m/day

rho_0=1027; %density of water [kg m-3]
rho_i=920;  %density of ice [kg m-3]
S_i=10;      %Salinity of ice [psu] [Martin S, Kauffman P (1981) A field and laboratory study of wave damping by greaseice. J Glaciol 27:283–313]
h0=100;

    for t=1:(length(time))-time_frame
            S_0=sum(sal(3:21,t));
            S_f=sum(sal(3:21,t+time_frame));
            output(t) = (((rho_0*h0*S_0-rho_0*h0*S_f)/(rho_i*S_i-rho_0*S_f))/(time(t+time_frame)-time(t)));
    end

end

