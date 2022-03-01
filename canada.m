%% Solar radiation for a one year time for Saguenay, Canada
%_________________________________________________________
%Calculating the irradinace in Klagenfurt city over a year
clear all

L = 56;                                                 %Latitude
LoT = 14;
dc_arr = [];                                            %Declination angle array
hsr_arr = [];                                           %Hourly sunrise angle array
K = [];                                                 %K & A: constants
A = [];
Tt_arr_annual = [];                                     %Annual total irradiance
H_arr_year = [];

for N = 1:1:365                                         %Over the year
    dc = 23.45 * sind((360/365) * (N-81));
    dc_arr = [dc; dc_arr];

    hsr = acosd(-tand(dc) * tand(L)) / 15;
    hsr_arr = [hsr; hsr_arr];
    
    sunrise_hour = 12 - round(hsr);
    
    A = 1160 + (75 * sind((360/365) * (N-275)));
    K = 0.174 + (0.035 * sind((360/365) * (N-100)));
    C = 0.095 + (0.04 * sind((360/365) * (N-100)));

    Tt_arr = [];
    H_arr = [];
    
    day_length = (2/15) * ( acosd(-tand(dc) * tand(L)));
    sunset_hour = sunrise_hour + day_length
 
    for now = sunrise_hour:1:sunset_hour
        H = abs(now-12)*15;
        if now > 12 
            H = H * (-1);
        end
        
        beta = asind((sind(L) * sind(dc)) + (cosd(L) * cosd(dc) * cosd(H)));     
        Ib = A * (exp(-K / sind(beta)));
        C = 0.095 + (0.04 * sind((360/365) * (N-100))); 
        
        Ibc = Ib * sind(beta);                              %Indirect Solar Radiation
        Id = C * Ib;                                        %Direct Solar Radiation
        Tt = Ibc + Id;                                      %Total Irradiance
      
        Tt_arr = [Tt_arr; Tt];
        H_arr = [H_arr,H];
    end
    H_arr;
    Tt_arr_annual = [Tt_arr; Tt_arr_annual];
end

 Tt_arr_annual;
 Tt_arr_annual_i = [];
 
for c = 1:1:length(Tt_arr_annual)
    if Tt_arr_annual(c) < 0
        Tt_arr_annual(c) = 0;
            
    elseif Tt_arr_annual(c) > 1000
        Tt_arr_annual(c) = 1000;
    else
        Tt_arr_annual(c) = Tt_arr_annual(c);
    end
    Tt_arr_annual_i = [Tt_arr_annual(c); Tt_arr_annual_i];
end

Tt_arr_annual_i;
% plot(Tt_arr_annual_i)

%%Stand-alone PV system Design__________________________
%%System Specicifications
array_i = [];

for PV_Wp = 200:50:700                                      %Capacity of the PV arrays
    for MaxCharging_state = 200:50:700                                 %Battery capacity (kWh/day)
        PV_eff=0.16;                                        %PV module efficiency 
        Batt_Voltage=12;                                    %Battery Voltage (kanat V_B)
        Inv_RP=2500;                                        %Inverter ratE_net power
        DOD=0.8;                                            %AllowE_net depth of charge
        Charging_eff=0.8;                                   %Charging eff (kanat Charge_eff)
        Alpha= .05;     
        Wire_eff= 0.98;
        MinCharging_state=MaxCharging_state *(1-DOD);
        T_cell = 25;
       
        %Calculating the load based on an estimate for a small cluster 
        %of 20 houses with similar daily energy consumption
        Load_i = [90 90 90 90 90 150 300 300 450 450 450 400 400 400 200 200 580 580 600 90];
        Load = [];
        for i = 1:1:365
            Load = [Load_i,Load];
        end
        Load;                                                    

        E_array = PV_Wp * (Tt_arr_annual_i/1000);
        E_net = E_array - Load;
        SOCi=MaxCharging_state;
        SOCf=[];
        Deff=[];
        Dampf=[];

        for i=1:length(E_net);
            SOC= E_net(i)+SOCi;
    
            if (SOC > MaxCharging_state)
                Dampi=SOC-MaxCharging_state;
                Defi=0;
                SOCi=MaxCharging_state;
                
            elseif (SOC<MinCharging_state)
                SOCi=MinCharging_state;
                Defi=SOC-MinCharging_state;
                Dampi=0;
            else
                SOCi=SOC;
                Defi=0;
                Dampi=0;
            end
    
            SOCf=[SOCf; SOCi];
            Deff=[Deff; Defi];
            Dampf=[Dampf; Dampi];
        end
        SOCf;
        Deff;
        Dampf;
        SOC_per=SOCf./MaxCharging_state;
        LLP_calculatE_net=abs(sum(Deff))/(sum (L));

        array = [PV_Wp, MaxCharging_state, LLP_calculatE_net];
        array_i = [array; array_i]
    end
end

array_i
plot(array_i)