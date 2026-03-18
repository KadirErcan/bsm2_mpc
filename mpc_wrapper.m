function [KLa4_out, Qintr_out] = mpc_wrapper(S_NH, S_NO, S_O, FtoM)
    % This block runs at the sampling time (e.g., every 15 mins)
    
    % Call the Python function step_mpc from bsm2_mpc.py
    try
        % Pass states to Python
        u_opt = py.bsm2_mpc.step_mpc(S_NH, S_NO, S_O, FtoM);
        
        % Convert Python list back to MATLAB doubles
        u_array = double(u_opt);
        
        KLa4_out = u_array(1);
        Qintr_out = u_array(2);
        
    catch ME
        % Fallback to default PI setpoints if Python fails to solve
        disp('MPC failed to solve, falling back to defaults.');
        KLa4_out = 120.0; % Default reginit_bsm2.m value
        Qintr_out = 61944.0; % Default Qintr
    end
end