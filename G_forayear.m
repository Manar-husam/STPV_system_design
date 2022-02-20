%% Solar radiation for a one year time for Klagenfurt city
%____________________________________________________________
L = 47;                                                 %Latitude
dc_arr = [];                                            %Declination angle array
hsr_arr = [];                                           %Hourly sunrise angle array
K = [];                                                 % K & A: constants
A = [];
Tt_arr_annual = [];
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
plot(Tt_arr_annual_i)

%% Stand-alone PV system Design
%%System Specicifications
array_i = [];

for PV_Wp = 200:50:700                                      %Capacity of the PV arrays
    for SOCmax = 200:50:700                                 %Battery capacity (kWh/day)
        PV_eff=0.16;                                        %PV module efficiency 
        Batt_Voltage=12;                                    %Battery Voltage (kanat V_B)
        Inv_RP=2500;                                        %Inverter rated power
        DOD=0.8;                                            %Allowed depth of charge
        Charging_eff=0.8;                                   %Charging eff (kanat Charge_eff)
        Alpha= .05;     
        Wire_eff= 0.98;
        SOCmin=SOCmax *(1-DOD);

        %(3.1) Simualtion of the SAPV system
        % P_Ratio=(PV_Wp*(G/1000))/Inv_RP;
        % Inv_eff=97.644-(P_Ratio*1.995)- (0.445./P_Ratio); %5KW

        G = Tt_arr_annual_i;

        L_i = [90 90 90 90 90 150 300 300 450 450 450 400 400 400 200 200 580 580 600 90];
        L = [];
        for i = 1:1:365
            L = [L_i,L];
        end
        L;                                                    %LLP should be 0.05 (5%)

        E_PV= PV_Wp*(G/1000);
        ED=E_PV-L;
        SOCi=SOCmax;
        SOCf=[];
        Deff=[];
        Dampf=[];

%%(3.2)
        for i=1:length(ED);
            SOC= ED(i)+SOCi;
    
            if (SOC > SOCmax)
                Dampi=SOC-SOCmax;
                Defi=0;
                SOCi=SOCmax;
    %%(3.3)
            elseif (SOC<SOCmin)
                SOCi=SOCmin;
                Defi=SOC-SOCmin;
                Dampi=0;

    %%(3.4)
            else
                SOCi=SOC;
                Defi=0;
                Dampi=0;
            end
    
%%(3.5)
            SOCf=[SOCf; SOCi];
            Deff=[Deff; Defi];
            Dampf=[Dampf; Dampi];
        end
        SOCf;
        Deff;
        Dampf;
        SOC_per=SOCf./SOCmax;
        LLP_calculated=abs(sum(Deff))/(sum (L));

        array = [PV_Wp, SOCmax, LLP_calculated];
        array_i = [array; array_i]
    end
end

array_i
array_i_1 = [];
array_i_2 = [];
array_i_3 = [];

for a = 1:1:length(array_i)
    array_i_1 = [array_i(a,1); array_i_1];
    array_i_2 = [array_i(a,2); array_i_2]; 
    array_i_3 = [array_i(a,3); array_i_3]; 
end

array_i_1;
array_i_2;
array_i_3;
plot3(array_i_1, array_i_2, array_i_3)
array_short = [];


for i = 1:1:length(array_i)
    if array_i(i,3) < .05 && array_i(i,3) > .045
        array_short = [array_i(i,1),array_i(i,2), array_i(i,3); array_short];
        
    end
end
        
array_short

array_ti_1 = [];
array_ti_2 = [];
array_ti_3 = [];

for a = 1:1:length(array_short)
    array_ti_1 = [array_short(a,1); array_ti_1];
    array_ti_2 = [array_short(a,2); array_ti_2]; 
    array_ti_3 = [array_short(a,3); array_ti_3]; 
end

array_ti_1;
array_ti_2;
array_ti_3;
figure
plot3(array_ti_1, array_ti_2, array_ti_3)

CPv = .30;
CBatt = .80;
C_arr = [];

for e = 1:1:length(array_short)
    C_arr = [(array_short(e,1) * CPv) + (array_short(e,2) * CBatt), array_short(e,3); C_arr];
end
C_arr


     
     
     