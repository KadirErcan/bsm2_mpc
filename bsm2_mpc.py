# bsm2_mpc.py
import sys
from types import ModuleType

# 1. Create a "Chameleon" object that pretends to be anything and everything
class Chameleon(ModuleType):
    def __getattr__(self, name): return Chameleon(name)
    def __call__(self, *args, **kwargs): return Chameleon('call')
    __path__ = [] 

# 2. The Interceptor: This catches ANY import starting with 'matplotlib'
class MatplotlibInterceptor:
    def find_spec(self, fullname, path, target=None):
        if fullname.startswith("matplotlib"):
            return self._gen_spec(fullname)
        return None
    def _gen_spec(self, name):
        from importlib.machinery import ModuleSpec
        return ModuleSpec(name, self)
    def create_module(self, spec):
        return Chameleon(spec.name)
    def exec_module(self, module):
        pass

# 3. Inject the interceptor into the very front of Python's search brain
sys.meta_path.insert(0, MatplotlibInterceptor())

# NOW import the libraries you actually need for the math
import numpy as np
import casadi as ca
import do_mpc

class BSM2_LCA_MPC:
    def __init__(self):
        # 1. Initialize Model
        self.model = do_mpc.model.Model('continuous')
        
        # Define States (_x)
        self.S_NH_4 = self.model.set_variable(var_type='_x', var_name='S_NH_4')
        self.S_NO_4 = self.model.set_variable(var_type='_x', var_name='S_NO_4')
        self.S_O_4  = self.model.set_variable(var_type='_x', var_name='S_O_4')
        self.F_to_M = self.model.set_variable(var_type='_x', var_name='F_to_M')
        
        # Define Manipulated Variables (_u)
        self.KLa4   = self.model.set_variable(var_type='_u', var_name='KLa4')
        self.Q_intr = self.model.set_variable(var_type='_u', var_name='Q_intr')
        
        # --- THE MISSING PART: Define the Dynamics (RHS) ---
        # For now, we tell the MPC that the states are simply controlled by the inputs
        # This satisfies the solver so it can "see" the math connection.

        # --- THE REAL ASM1 BIOLOGY ---
        
        # 1. Biological & Plant Parameters (Standard ASM1 values)
        mu_A = 0.5          # Autotrophic max growth rate (1/day)
        K_NH = 1.0          # Ammonia half-saturation (mg/L)
        K_OA = 0.4          # Oxygen half-saturation (mg/L)
        Y_A  = 0.24         # Autotrophic yield
        i_XB = 0.08         # Nitrogen fraction in biomass
        
        # We assume X_BA (Autotrophic Bacteria) is relatively stable around 150 mg/L in Tank 4
        X_BA = 150.0        
        
        # 2. Plant Physical Constants
        V4 = 1333.0         # Volume of Tank 4 (m3)
        S_O_sat = 8.0       # Oxygen saturation at 15°C (mg/L)
        Q_in = 18446.0      # Average influent flow (m3/day)
        
        # 3. The Biological Engine (Nitrification)
        # This tells the MPC: "Bacteria consume Ammonia and Oxygen simultaneously"
        proc3 = mu_A * (self.S_NH_4 / (K_NH + self.S_NH_4)) * (self.S_O_4 / (K_OA + self.S_O_4)) * X_BA
        
        # 4. Write the Differential Equations (The Rules of the Plant)
        
        # OXYGEN: Added by the blowers (KLa), consumed by the bacteria (proc3)
        rhs_S_O_4 = self.KLa4 * (S_O_sat - self.S_O_4) - ((4.57 - Y_A) / Y_A) * proc3
        
        # AMMONIA: Flows in from Tank 3 (assume ~8.0 mg/L), destroyed by bacteria
        rhs_S_NH_4 = ((Q_in + self.Q_intr) / V4) * (8.0 - self.S_NH_4) - (i_XB + (1.0 / Y_A)) * proc3
        
        # NITRATE: Created by bacteria, pumped away by Internal Recirculation (Q_intr)
        rhs_S_NO_4 = (1.0 / Y_A) * proc3 - (self.Q_intr / V4) * self.S_NO_4
        
        # F_to_M: Keeping your variable stable
        rhs_F_to_M = -0.01 * self.F_to_M

        # 5. Apply to the MPC Model
        self.model.set_rhs('S_NH_4', rhs_S_NH_4)
        self.model.set_rhs('S_O_4',  rhs_S_O_4)
        self.model.set_rhs('S_NO_4', rhs_S_NO_4)
        self.model.set_rhs('F_to_M', rhs_F_to_M)
        
        self.model.setup()
        
        # 2. Initialize Controller
        self.mpc = do_mpc.controller.MPC(self.model)
        setup_mpc = {
            'n_horizon': 24,
            't_step': 0.0104, # 15 minutes in days
            'n_robust': 0,
            'store_full_solution': False,
        }
        self.mpc.set_param(**setup_mpc)
        
        # 3. Objective Function (LCA Minimization)
        # Minimize Aeration (KLa) and Recirculation (Q_intr)
        lterm = (1.0 * self.KLa4) + (0.0001 * self.Q_intr) + (10000.0 * (self.S_NH_4 - 1.5)**2) + (10000.0 * (self.S_NO_4 - 10.0)**2)
        
        mterm = (10000.0 * (self.S_NH_4 - 1.5)**2) + (10000.0 * (self.S_NO_4 - 10.0)**2)
        
        self.mpc.set_objective(mterm=mterm, lterm=lterm)
        self.mpc.set_rterm(KLa4=0.1, Q_intr=0.01)
        
        # 4. Constraints
        # 4. Constraints (Physical Limits)
        # We prevent the MPC from ever picking negative values
        self.mpc.bounds['lower', '_u', 'KLa4'] = 30.0      # Min aeration
        self.mpc.bounds['upper', '_u', 'KLa4'] = 360.0     # Max aeration
        self.mpc.bounds['lower', '_u', 'Q_intr'] = 0.0    # Min recirculation
        self.mpc.bounds['upper', '_u', 'Q_intr'] = 92000.0 # Max recirculation
        
        self.mpc.setup()

    def get_action(self, current_states):
        x0 = np.array(current_states).reshape(-1, 1)
        u0 = self.mpc.make_step(x0)
        # Return KLa4 and Q_intr as a basic Python list
        return [float(u0[0]), float(u0[1])]

# Global instance so MATLAB doesn't rebuild the MPC every 15 minutes
mpc_instance = BSM2_LCA_MPC()

def step_mpc(S_NH_4, S_NO_4, S_O_4, F_to_M):
    states = [S_NH_4, S_NO_4, S_O_4, F_to_M]
    return mpc_instance.get_action(states)