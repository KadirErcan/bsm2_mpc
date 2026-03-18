asm1init_bsm2;
settler1dinit_bsm2;
hyddelayinit_bsm2;
primclarinit_bsm2;
thickenerinit_bsm2;
dewateringinit_bsm2;
storageinit_bsm2;
adm1init_bsm2;  % also includes settings for AS/AD and AD/AS interfaces
reginit_bsm2;
load ./Influent_data/sensornoise_bsm2;
sensorinit_bsm2;

load ./Influent_data/constinfluent_bsm2;
load ./Influent_data/dyninfluent_bsm2;


% General parameter for all subsystems
% TEMPMODEL: 0 - influent wastewater temperature is just passed through process reactors 
%            1 - mass balance for the wastewater temperature is used in
%            process reactors
% Note: thickener and dewatering are ideal models, i.e. no impact since T_out =
% T_in, in flow splitters T_out = T_in and in flow combiners mass balance based heat balance is always used.

TEMPMODEL = [ 1 ]; 

% to activate calculation of  dummy states in settler, AS reactors and storage tank set ACTIVATE = 1
ACTIVATE = [ 0 ];