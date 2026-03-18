disp(' ')
disp(' ')
disp('This script will provide the final values of almost all state variables')
disp('of the last BSM2 simulation.')
disp(' ')
disp(' ')
[m n]=size(reac1);

disp('Influent characteristics')
disp('************************')
disp(['   SI =  ', num2str(in(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(in(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(in(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(in(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(in(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(in(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(in(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(in(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(in(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(in(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(in(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(in(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(in(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(in(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(in(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(in(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(in(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(in(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(in(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(in(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(in(m,21)), ' mg xxx/l']);
end;

disp(' ')
disp('Flow conditions')
disp('***************')
disp(['   Influent flow rate to WWTP = ', num2str(in(m,15)), ' m3/d']);
disp(['   Primary clarifier feed flow rate = ', num2str(primaryin(m,15)), ' m3/d']);
disp(['   Primary clarifier sludge flow rate = ', num2str(primaryout(m,36)), ' m3/d']);
disp(['   Primary clarifier effluent flow rate = ', num2str(primaryout(m,15)), ' m3/d']);
disp(['   Influent flow rate to AS = ', num2str(ASinput(m,15)), ' m3/d']);
disp(['   Internal recirculation flow rate = ', num2str(rec(m,15)), ' m3/d']);
disp(['   Settler feed flow rate = ', num2str(feed(m,15)), ' m3/d']);
disp(['   Settler effluent flow rate = ', num2str(settler(m,37)), ' m3/d']);
disp(['   Returned sludge flow rate = ', num2str(settler(m,15)), ' m3/d']);
disp(['   Wastage sludge flow rate = ', num2str(settler(m,22)), ' m3/d']);
disp(['   Thickener feed flow rate = ', num2str(thickenerin(m,15)), ' m3/d']);
disp(['   Thickener sludge flow rate = ', num2str(thickenerout(m,15)), ' m3/d']);
disp(['   Thickener effluent flow rate = ', num2str(thickenerout(m,36)), ' m3/d']);
disp(['   Digester feed flow rate = ', num2str(digesterinpreinterface(m,15)), ' m3/d']);
disp(['   Digester output flow rate = ', num2str(digesteroutpostinterface(m,15)), ' m3/d']);
disp(['   Dewatering feed flow rate = ', num2str(dewateringin(m,15)), ' m3/d']);
disp(['   Dewatering sludge flow rate = ', num2str(dewateringout(m,15)), ' m3/d']);
disp(['   Dewatering effluent flow rate = ', num2str(dewateringout(m,36)), ' m3/d']);
disp(['   Storage tank feed flow rate (not bypass) = ', num2str(storagein(m,15)), ' m3/d']);
disp(['   Storage tank effluent flow rate = ', num2str(storageout(m,15)), ' m3/d']);
disp(['   Storage tank bypass flow rate = ', num2str(storagebypass(m,15)), ' m3/d']);
disp(['   Effluent flow rate = ', num2str(effluent(m,15)), ' m3/d']);
disp(['   Sludge disposal flow rate = ', num2str(sludge(m,15)), ' m3/d']);

disp(' ')
disp('Primary clarifier effluent')
disp('**************************')
disp(['   SI =  ', num2str(primaryout(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(primaryout(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(primaryout(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(primaryout(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(primaryout(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(primaryout(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(primaryout(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(primaryout(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(primaryout(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(primaryout(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(primaryout(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(primaryout(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(primaryout(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(primaryout(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(primaryout(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(primaryout(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(primaryout(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(primaryout(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(primaryout(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(primaryout(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(primaryout(m,21)), ' mg xxx/l']);
end;

disp(' ')
disp('Primary clarifier underflow')
disp('***************************')
disp(['   SI =  ', num2str(primaryout(m,22)), ' mg COD/l']);
disp(['   SS =  ', num2str(primaryout(m,23)), ' mg COD/l']);
disp(['   XI =  ', num2str(primaryout(m,24)), ' mg COD/l']);
disp(['   XS =  ', num2str(primaryout(m,25)), ' mg COD/l']);
disp(['   XBH = ', num2str(primaryout(m,26)), ' mg COD/l']);
disp(['   XBA = ', num2str(primaryout(m,27)), ' mg COD/l']);
disp(['   XP =  ', num2str(primaryout(m,28)), ' mg COD/l']);
disp(['   SO =  ', num2str(primaryout(m,29)), ' mg -COD/l']);
disp(['   SNO = ', num2str(primaryout(m,30)), ' mg N/l']);
disp(['   SNH = ', num2str(primaryout(m,31)), ' mg N/l']);
disp(['   SND = ', num2str(primaryout(m,32)), ' mg N/l']);
disp(['   XND = ', num2str(primaryout(m,33)), ' mg N/l']);
disp(['   SALK = ', num2str(primaryout(m,34)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(primaryout(m,35)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(primaryout(m,36)), ' m3/d']);
disp(['   Temperature = ', num2str(primaryout(m,37)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(primaryout(m,38)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(primaryout(m,39)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(primaryout(m,40)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(primaryout(m,41)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(primaryout(m,42)), ' mg xxx/l']);
end;

disp(' ')
disp('Input to AS')
disp('***********')
disp(['   SI =  ', num2str(ASinput(m,1)/ASinput(m,15)), ' mg COD/l']);
disp(['   SS =  ', num2str(ASinput(m,2)/ASinput(m,15)), ' mg COD/l']);
disp(['   XI =  ', num2str(ASinput(m,3)/ASinput(m,15)), ' mg COD/l']);
disp(['   XS =  ', num2str(ASinput(m,4)/ASinput(m,15)), ' mg COD/l']);
disp(['   XBH = ', num2str(ASinput(m,5)/ASinput(m,15)), ' mg COD/l']);
disp(['   XBA = ', num2str(ASinput(m,6)/ASinput(m,15)), ' mg COD/l']);
disp(['   XP =  ', num2str(ASinput(m,7)/ASinput(m,15)), ' mg COD/l']);
disp(['   SO =  ', num2str(ASinput(m,8)/ASinput(m,15)), ' mg -COD/l']);
disp(['   SNO = ', num2str(ASinput(m,9)/ASinput(m,15)), ' mg N/l']);
disp(['   SNH = ', num2str(ASinput(m,10)/ASinput(m,15)), ' mg N/l']);
disp(['   SND = ', num2str(ASinput(m,11)/ASinput(m,15)), ' mg N/l']);
disp(['   XND = ', num2str(ASinput(m,12)/ASinput(m,15)), ' mg N/l']);
disp(['   SALK = ', num2str(ASinput(m,13)/ASinput(m,15)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(ASinput(m,14)/ASinput(m,15)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(ASinput(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(ASinput(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(ASinput(m,17)/ASinput(m,15)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(ASinput(m,18)/ASinput(m,15)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(ASinput(m,19)/ASinput(m,15)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(ASinput(m,20)/ASinput(m,15)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(ASinput(m,21)/ASinput(m,15)), ' mg xxx/l']);
end;


disp(' ')
disp('AS reactor 1')
disp('************')
disp(['   SI =  ', num2str(reac1(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(reac1(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(reac1(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(reac1(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(reac1(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(reac1(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(reac1(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(reac1(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(reac1(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(reac1(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(reac1(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(reac1(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(reac1(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(reac1(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(reac1(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(reac1(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(reac1(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(reac1(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(reac1(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(reac1(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(reac1(m,21)), ' mg xxx/l']);
end;

disp(' ')
disp('AS reactor 2')
disp('************')
disp(['   SI =  ', num2str(reac2(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(reac2(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(reac2(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(reac2(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(reac2(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(reac2(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(reac2(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(reac2(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(reac2(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(reac2(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(reac2(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(reac2(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(reac2(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(reac2(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(reac2(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(reac2(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(reac2(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(reac2(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(reac2(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(reac2(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(reac2(m,21)), ' mg xxx/l']);
end;

disp(' ')
disp('AS reactor 3')
disp('************')
disp(['   SI =  ', num2str(reac3(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(reac3(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(reac3(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(reac3(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(reac3(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(reac3(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(reac3(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(reac3(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(reac3(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(reac3(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(reac3(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(reac3(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(reac3(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(reac3(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(reac3(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(reac3(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(reac3(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(reac3(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(reac3(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(reac3(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(reac3(m,21)), ' mg xxx/l']);
end;

disp(' ')
disp('AS reactor 4')
disp('************')
disp(['   SI =  ', num2str(reac4(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(reac4(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(reac4(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(reac4(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(reac4(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(reac4(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(reac4(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(reac4(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(reac4(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(reac4(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(reac4(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(reac4(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(reac4(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(reac4(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(reac4(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(reac4(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(reac4(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(reac4(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(reac4(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(reac4(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(reac4(m,21)), ' mg xxx/l']);
end;

disp(' ')
disp('AS reactor 5')
disp('************')
disp(['   SI =  ', num2str(reac5(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(reac5(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(reac5(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(reac5(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(reac5(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(reac5(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(reac5(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(reac5(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(reac5(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(reac5(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(reac5(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(reac5(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(reac5(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(reac5(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(reac5(m,15)), ' m3/d']);
disp(['   temperature = ', num2str(reac5(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(reac5(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(reac5(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(reac5(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(reac5(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(reac5(m,21)), ' mg xxx/l']);
end;

disp(' ')
disp('Settler underflow')
disp('*****************')
disp(['   SI =  ', num2str(settler(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(settler(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(settler(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(settler(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(settler(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(settler(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(settler(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(settler(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(settler(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(settler(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(settler(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(settler(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(settler(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(settler(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(settler(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(settler(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(settler(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(settler(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(settler(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(settler(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(settler(m,21)), ' mg xxx/l']);
end;

disp(' ')
disp('Settler effluent')
disp('****************')
disp(['   SI =  ', num2str(settler(m,23)), ' mg COD/l']);
disp(['   SS =  ', num2str(settler(m,24)), ' mg COD/l']);
disp(['   XI =  ', num2str(settler(m,25)), ' mg COD/l']);
disp(['   XS =  ', num2str(settler(m,26)), ' mg COD/l']);
disp(['   XBH = ', num2str(settler(m,27)), ' mg COD/l']);
disp(['   XBA = ', num2str(settler(m,28)), ' mg COD/l']);
disp(['   XP =  ', num2str(settler(m,29)), ' mg COD/l']);
disp(['   SO =  ', num2str(settler(m,30)), ' mg -COD/l']);
disp(['   SNO = ', num2str(settler(m,31)), ' mg N/l']);
disp(['   SNH = ', num2str(settler(m,32)), ' mg N/l']);
disp(['   SND = ', num2str(settler(m,33)), ' mg N/l']);
disp(['   XND = ', num2str(settler(m,34)), ' mg N/l']);
disp(['   SALK = ', num2str(settler(m,35)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(settler(m,36)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(settler(m,37)), ' m3/d']);
disp(['   Temperature = ', num2str(settler(m,38)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(settler(m,39)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(settler(m,40)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(settler(m,41)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(settler(m,42)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(settler(m,43)), ' mg xxx/l']);
end;

disp(' ')
disp('Settler internal (1 is top layer)')
disp('*********************************')
disp(['   TSS1 = ', num2str(settler(m,44)), ' mg SS/l']);
disp(['   TSS2 = ', num2str(settler(m,45)), ' mg SS/l']);
disp(['   TSS3 = ', num2str(settler(m,46)), ' mg SS/l']);
disp(['   TSS4 = ', num2str(settler(m,47)), ' mg SS/l']);
disp(['   TSS5 = ', num2str(settler(m,48)), ' mg SS/l']);
disp(['   TSS6 = ', num2str(settler(m,49)), ' mg SS/l']);
disp(['   TSS7 = ', num2str(settler(m,50)), ' mg SS/l']);
disp(['   TSS8 = ', num2str(settler(m,51)), ' mg SS/l']);
disp(['   TSS9 = ', num2str(settler(m,52)), ' mg SS/l']);
disp(['   TSS10 = ', num2str(settler(m,53)), ' mg SS/l']);

disp(' ')
disp('Thickener effluent')
disp('******************')
disp(['   SI =  ', num2str(thickenerout(m,22)), ' mg COD/l']);
disp(['   SS =  ', num2str(thickenerout(m,23)), ' mg COD/l']);
disp(['   XI =  ', num2str(thickenerout(m,24)), ' mg COD/l']);
disp(['   XS =  ', num2str(thickenerout(m,25)), ' mg COD/l']);
disp(['   XBH = ', num2str(thickenerout(m,26)), ' mg COD/l']);
disp(['   XBA = ', num2str(thickenerout(m,27)), ' mg COD/l']);
disp(['   XP =  ', num2str(thickenerout(m,28)), ' mg COD/l']);
disp(['   SO =  ', num2str(thickenerout(m,29)), ' mg -COD/l']);
disp(['   SNO = ', num2str(thickenerout(m,30)), ' mg N/l']);
disp(['   SNH = ', num2str(thickenerout(m,31)), ' mg N/l']);
disp(['   SND = ', num2str(thickenerout(m,32)), ' mg N/l']);
disp(['   XND = ', num2str(thickenerout(m,33)), ' mg N/l']);
disp(['   SALK = ', num2str(thickenerout(m,34)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(thickenerout(m,35)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(thickenerout(m,36)), ' m3/d']);
disp(['   Temperature = ', num2str(thickenerout(m,37)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(thickenerout(m,38)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(thickenerout(m,39)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(thickenerout(m,40)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(thickenerout(m,41)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(thickenerout(m,42)), ' mg xxx/l']);
end;

disp(' ')
disp('Thickener underflow')
disp('*******************')
disp(['   SI =  ', num2str(thickenerout(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(thickenerout(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(thickenerout(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(thickenerout(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(thickenerout(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(thickenerout(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(thickenerout(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(thickenerout(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(thickenerout(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(thickenerout(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(thickenerout(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(thickenerout(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(thickenerout(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(thickenerout(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(thickenerout(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(thickenerout(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(thickenerout(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(thickenerout(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(thickenerout(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(thickenerout(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(thickenerout(m,21)), ' mg xxx/l']);
end;

disp(' ')
disp('Anaerobic digester influent (pre ASM2ADM interface)')
disp('***************************************************')
disp(['   SI =  ', num2str(digesterinpreinterface(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(digesterinpreinterface(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(digesterinpreinterface(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(digesterinpreinterface(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(digesterinpreinterface(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(digesterinpreinterface(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(digesterinpreinterface(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(digesterinpreinterface(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(digesterinpreinterface(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(digesterinpreinterface(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(digesterinpreinterface(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(digesterinpreinterface(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(digesterinpreinterface(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(digesterinpreinterface(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(digesterinpreinterface(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(digesterinpreinterface(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(digesterinpreinterface(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(digesterinpreinterface(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(digesterinpreinterface(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(digesterinpreinterface(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(digesterinpreinterface(m,21)), ' mg xxx/l']);
end;

printadmresults(digesterin, digesterout, DIGESTERPAR, ACTIVATE);

disp('Anaerobic digester output (post ADM2ASM interface)')
disp('**************************************************')
disp(['   SI =  ', num2str(digesteroutpostinterface(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(digesteroutpostinterface(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(digesteroutpostinterface(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(digesteroutpostinterface(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(digesteroutpostinterface(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(digesteroutpostinterface(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(digesteroutpostinterface(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(digesteroutpostinterface(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(digesteroutpostinterface(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(digesteroutpostinterface(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(digesteroutpostinterface(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(digesteroutpostinterface(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(digesteroutpostinterface(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(digesteroutpostinterface(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(digesteroutpostinterface(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(digesteroutpostinterface(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(digesteroutpostinterface(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(digesteroutpostinterface(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(digesteroutpostinterface(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(digesteroutpostinterface(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(digesteroutpostinterface(m,21)), ' mg xxx/l']);
end;

disp(' ')
disp('Dewatering effluent (reject water)')
disp('**********************************')
disp(['   SI =  ', num2str(dewateringout(m,22)), ' mg COD/l']);
disp(['   SS =  ', num2str(dewateringout(m,23)), ' mg COD/l']);
disp(['   XI =  ', num2str(dewateringout(m,24)), ' mg COD/l']);
disp(['   XS =  ', num2str(dewateringout(m,25)), ' mg COD/l']);
disp(['   XBH = ', num2str(dewateringout(m,26)), ' mg COD/l']);
disp(['   XBA = ', num2str(dewateringout(m,27)), ' mg COD/l']);
disp(['   XP =  ', num2str(dewateringout(m,28)), ' mg COD/l']);
disp(['   SO =  ', num2str(dewateringout(m,29)), ' mg -COD/l']);
disp(['   SNO = ', num2str(dewateringout(m,30)), ' mg N/l']);
disp(['   SNH = ', num2str(dewateringout(m,31)), ' mg N/l']);
disp(['   SND = ', num2str(dewateringout(m,32)), ' mg N/l']);
disp(['   XND = ', num2str(dewateringout(m,33)), ' mg N/l']);
disp(['   SALK = ', num2str(dewateringout(m,34)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(dewateringout(m,35)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(dewateringout(m,36)), ' m3/d']);
disp(['   Temperature = ', num2str(dewateringout(m,37)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(dewateringout(m,38)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(dewateringout(m,39)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(dewateringout(m,40)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(dewateringout(m,41)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(dewateringout(m,42)), ' mg xxx/l']);
end;

disp(' ')
disp('Storage tank effluent')
disp('*********************')
disp(['   SI =  ', num2str(storageout(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(storageout(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(storageout(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(storageout(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(storageout(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(storageout(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(storageout(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(storageout(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(storageout(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(storageout(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(storageout(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(storageout(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(storageout(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(storageout(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(storageout(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(storageout(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(storageout(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(storageout(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(storageout(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(storageout(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(storageout(m,21)), ' mg xxx/l']);
end;
disp(['   Storage tank water volume = ', num2str(storageout(m,22)), ' m3']);


disp(' ')
disp('Storage tank bypass')
disp('*******************')
disp(['   SI =  ', num2str(storagebypass(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(storagebypass(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(storagebypass(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(storagebypass(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(storagebypass(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(storagebypass(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(storagebypass(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(storagebypass(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(storagebypass(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(storagebypass(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(storagebypass(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(storagebypass(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(storagebypass(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(storagebypass(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(storagebypass(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(storagebypass(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(storagebypass(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(storagebypass(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(storagebypass(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(storagebypass(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(storagebypass(m,21)), ' mg xxx/l']);
end;


disp(' ')
disp('Storage tank output + bypass')
disp('****************************')
disp(['   SI =  ', num2str(storagetotout(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(storagetotout(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(storagetotout(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(storagetotout(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(storagetotout(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(storagetotout(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(storagetotout(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(storagetotout(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(storagetotout(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(storagetotout(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(storagetotout(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(storagetotout(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(storagetotout(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(storagetotout(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(storagetotout(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(storagetotout(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(storagetotout(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(storagetotout(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(storagetotout(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(storagetotout(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(storagetotout(m,21)), ' mg xxx/l']);
end;

disp(' ')
disp('WWTP effluent')
disp('*************')
disp(['   SI =  ', num2str(effluent(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(effluent(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(effluent(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(effluent(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(effluent(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(effluent(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(effluent(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(effluent(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(effluent(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(effluent(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(effluent(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(effluent(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(effluent(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(effluent(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(effluent(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(effluent(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(effluent(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(effluent(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(effluent(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(effluent(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(effluent(m,21)), ' mg xxx/l']);
end;

disp(' ')
disp('WWTP sludge disposal')
disp('********************')
disp(['   SI =  ', num2str(sludge(m,1)), ' mg COD/l']);
disp(['   SS =  ', num2str(sludge(m,2)), ' mg COD/l']);
disp(['   XI =  ', num2str(sludge(m,3)), ' mg COD/l']);
disp(['   XS =  ', num2str(sludge(m,4)), ' mg COD/l']);
disp(['   XBH = ', num2str(sludge(m,5)), ' mg COD/l']);
disp(['   XBA = ', num2str(sludge(m,6)), ' mg COD/l']);
disp(['   XP =  ', num2str(sludge(m,7)), ' mg COD/l']);
disp(['   SO =  ', num2str(sludge(m,8)), ' mg -COD/l']);
disp(['   SNO = ', num2str(sludge(m,9)), ' mg N/l']);
disp(['   SNH = ', num2str(sludge(m,10)), ' mg N/l']);
disp(['   SND = ', num2str(sludge(m,11)), ' mg N/l']);
disp(['   XND = ', num2str(sludge(m,12)), ' mg N/l']);
disp(['   SALK = ', num2str(sludge(m,13)), ' mol HCO3/m3']);
disp(['   TSS = ', num2str(sludge(m,14)), ' mg SS/l']);
disp(['   Flow rate = ', num2str(sludge(m,15)), ' m3/d']);
disp(['   Temperature = ', num2str(sludge(m,16)), ' degC']);
if (ACTIVATE > 0.5)
    disp(['   S_DUMMY_1 =  ', num2str(sludge(m,17)), ' mg xxx/l']);
    disp(['   S_DUMMY_2 =  ', num2str(sludge(m,18)), ' mg xxx/l']);
    disp(['   S_DUMMY_3 =  ', num2str(sludge(m,19)), ' mg xxx/l']);
    disp(['   X_DUMMY_1 =  ', num2str(sludge(m,20)), ' mg xxx/l']);
    disp(['   X_DUMMY_2 =  ', num2str(sludge(m,21)), ' mg xxx/l']);
end;

reactorvol = VOL1+VOL2+VOL3+VOL4+VOL5;
settlervol = DIM(1)*DIM(2);
totalvol = reactorvol+settlervol;

sludge1 = reac1(m,14)*VOL1 + reac2(m,14)*VOL2 + reac3(m,14)*VOL3 + reac4(m,14)*VOL4 + reac5(m,14)*VOL5;
sludge2 = (reac1(m,5)+reac1(m,6))*VOL1 + (reac2(m,5)+reac2(m,6))*VOL2 + (reac3(m,5)+reac3(m,6))*VOL3 + (reac4(m,5)+reac4(m,6))*VOL4 + (reac5(m,5)+reac5(m,6))*VOL5 + sum(((reac5(m,5)+reac5(m,6)).*(settler(m,44:53)./reac5(m,14))).*(settlervol/nooflayers));
sludge3 = sludge1 + sum(settler(m,44:53).*(settlervol/nooflayers));
waste1 = settler(m,14)*settler(m,22) + settler(m,44)*settler(m,37);
waste2 = sum(settler(m,5:6).*settler(m,22)) + sum(settler(m,27:28).*settler(m,37));
sludge_age1 = sludge1/waste1;
sludge_age2 = sludge2/waste2;
sludge_age3 = sludge3/waste1;

disp(' ')
disp('Other variables')
disp('***************')
disp(['   Trad. sludge age (XS + XP + XI + XBH + XBA in reactors) = ', num2str(sludge_age1), ' days']);
disp(['   Spec. sludge age (XBH + XBA in reactors and settler) = ', num2str(sludge_age2), ' days']);
disp(['   Spec. sludge age 2 (XS + XP + XI + XBH + XBA in reactors and settler) = ', num2str(sludge_age3), ' days']);
disp(['   Primary clarifier hydraulic retention time = ', num2str(VOL_P/(in(m,15))*24), ' hours']);
disp(['   Hydraulic retention time in AS + settler = ', num2str(totalvol/(in(m,15))*24), ' hours']);
disp(['   AS reactors hydraulic retention time = ', num2str(reactorvol/(in(m,15))*24), ' hours']);
disp(['   Anaerobic digester retention time = ', num2str(DIM_D(1)/(digesterinpreinterface(m,15))), ' days']);
disp(['   Thickening factor at bottom of settler (TSSu/TSSfeed) = ', num2str(settler(m,54))]);
disp(['   Thinning factor at top of settler (TSSeff/TSSfeed) = ', num2str(settler(m,55))]);
disp(['   Sludge blanket level (secondary clarifier) = ', num2str(settler(m,166)), ' meters']);

disp(' ')
disp('Dimensions')
disp('**********')
disp(['   Primary clarifier volume = ', num2str(VOL_P), ' m3']);
if reac1(m,8)<0.2
disp('   Reactor 1 is anoxic')
else
disp('   Reactor 1 is aerobic')
end
disp(['   Volume reactor 1 = ', num2str(VOL1), ' m3']);
if reac2(m,8)<0.2
disp('   Reactor 2 is anoxic')
else
disp('   Reactor 2 is aerobic')
end
disp(['   Volume reactor 2 = ', num2str(VOL2), ' m3']);
if reac3(m,8)<0.2
disp('   Reactor 3 is anoxic')
else
disp('   Reactor 3 is aerobic')
end
disp(['   Volume reactor 3 = ', num2str(VOL3), ' m3']);
if reac4(m,8)<0.2
disp('   Reactor 4 is anoxic')
else
disp('   Reactor 4 is aerobic')
end
disp(['   Volume reactor 4 = ', num2str(VOL4), ' m3']);
if reac5(m,8)<0.2
disp('   Reactor 5 is anoxic')
else
disp('   Reactor 5 is aerobic')
end
disp(['   Volume reactor 5 = ', num2str(VOL5), ' m3']);
disp(['   Settler height = ', num2str(DIM(2)), ' m']);
disp(['   Settler area = ', num2str(DIM(1)), ' m2']);
disp(['   Settler volume = ', num2str(settlervol), ' m3']);
disp(['   Anaerobic digester volume (liquid) = ', num2str(DIM_D(1)), ' m3']);
disp(['   Anaerobic digester head space volume = ', num2str(DIM_D(2)), ' m3']);
disp(['   Total storage tank volume (10-90% used) = ', num2str(VOL_S), ' m3']);

