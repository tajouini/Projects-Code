function [out_DA,out_CT] = load_data(T_DA, T_CT, noise)

%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

load ('../data/8760h-01p/DA_uncertainties.mat');
load ('../data/8760h-01p/CT_uncertainties.mat');
load ('../data/8760h-01p/CT_real.mat');
load ('../data/8760h-01p/DA_real.mat');

if noise == 0  % noise case
       
     
                  PL_CT= PL_CT_real; 
                  PV_CT = PV_CT_real; 
                  WT_CT = WT_CT_real;
                  HL_CT = HL_CT_real;
                  v_CT = v_CT_real; 
                  u_CT = u_CT_real;

                  PL_DA= PL_DA_real; 
                  PV_DA = PV_DA_real; 
                  WT_DA = WT_DA_real;
                  HL_DA = HL_DA_real;
                  v_DA = v_DA_real;  
                  u_DA = u_DA_real;

else 

    if T_DA == 24

                  PL_DA= PL_DA_N24; 
                  PV_DA = PV_DA_N24; 
                  WT_DA = WT_DA_N24;
                  HL_DA = HL_DA_N24;
                  v_DA = v_DA_N24;
                  u_DA = u_DA_N24;
    else 
        
        if T_DA == 48
  
                  PL_DA= PL_DA_N48; 
                  PV_DA = PV_DA_N48; 
                  WT_DA = WT_DA_N48;
                  HL_DA = HL_DA_N48;
                  v_DA = v_DA_N48;  
                  u_DA = u_DA_N48;   

        end   
    end  


     if T_CT == 1
       
                  PL_CT= PL_CT_N4; 
                  PV_CT = PV_CT_N4; 
                  WT_CT = WT_CT_N4;
                  HL_CT = HL_CT_N4;
                  v_CT = v_CT_N4;  
                  u_CT = u_CT_N4;   

     else 
         if T_CT == 12
             % CT measurement data 
                  PL_CT= PL_CT_N48; 
                  PV_CT = PV_CT_N48; 
                  WT_CT = WT_CT_N48;
                  HL_CT = HL_CT_N48;
                  v_CT = v_CT_N48;   
                  u_CT = u_CT_N48;   
         
         end
         if T_CT == 24
              % CT measurement data 
                   PL_CT= PL_CT_N96; 
                  PV_CT = PV_CT_N96; 
                  WT_CT = WT_CT_N96;
                  HL_CT = HL_CT_N96;
                  v_CT = v_CT_N96;   
                  u_CT = u_CT_N96;   

         end
     end
end



 out_DA = [PL_DA; PV_DA; WT_DA;HL_DA;v_DA; u_DA];   

 out_CT = [PL_CT; PV_CT; WT_CT;HL_CT;v_CT; u_CT];   

end
