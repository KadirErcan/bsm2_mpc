% Controller performance module for BSM2 (default closed loop)
%
% cut away first and last samples, i.e. t=smaller than starttime and 
% t = larger than stoptime
% 2008-07-30 file updated, UJ
% NOTE: this file must obviously be updated for every new control strategy tested
%
% Copyright: Ulf Jeppsson, IEA, Lund University, Lund, Sweden

start=clock; 
disp(' ')
disp('***** Evaluation of BSM2 controllers initiated *****')
disp(['Start time for BSM2 controller evaluation (hour:min:sec) = ', num2str(round(start(4:6)))]); %Display start time of evaluation
disp(' ')

starttime = 245;
stoptime = 609;
startindex=max(find(t <= starttime));
stopindex=min(find(t >= stoptime));

time_eval=t(startindex:stopindex);

sampletime = time_eval(2)-time_eval(1);
totalt=time_eval(end)-time_eval(1);

n=length(time_eval);

%cut out the parts of the files to be used
reac3part=reac3(startindex:(stopindex-1),:);
reac4part=reac4(startindex:(stopindex-1),:);
reac5part=reac5(startindex:(stopindex-1),:);
settlerpart=settler(startindex:(stopindex-1),:);
Qwregpart=Qwreg(startindex:(stopindex),:);
SO4regpart=SO4reg(startindex:(stopindex),:);
SO4sensorpart=SO4sensor(startindex:(stopindex),:);
kla3inpart=kla3in(startindex:(stopindex),:);
kla4inpart=kla4in(startindex:(stopindex),:);
kla5inpart=kla5in(startindex:(stopindex),:);

SO4refvec=ones(n-1,1)*SO4ref;
Qwrefvecfull=[ones(96*182,1)*Qw_low; ones(96*182,1)*Qw_high; ones(96*182,1)*Qw_low; ones(96*63,1)*Qw_high];
Qwrefvec = Qwrefvecfull(startindex:(stopindex-1),:);

% Controller performance
SO4error=SO4refvec-reac4part(:,8);
SO4errorsquare=SO4error.*SO4error;
deltauvec_SO4=abs(SO4regpart(2:end,3)-SO4regpart(1:(end-1),3));
Qwerror=Qwrefvec-settlerpart(:,22);
Qwerrorsquare=Qwerror.*Qwerror;
deltauvec_Qw=abs(Qwregpart(2:end,1)-Qwregpart(1:(end-1),1));

timevector = time_eval(2:end)-time_eval(1:(end-1));

SO4mean = mean(abs(SO4error));
SO4vec = SO4error.*timevector;
IAE_SO4vec=abs(SO4error).*timevector;
ISE_SO4vec=SO4errorsquare.*timevector;
Qwmean = mean(abs(Qwerror));
Qwvec = Qwerror.*timevector;
IAE_Qwvec=abs(Qwerror).*timevector;
ISE_Qwvec=Qwerrorsquare.*timevector;

deltauvecsquare_SO4=deltauvec_SO4.*deltauvec_SO4;
deltauvecsquare_Qw=deltauvec_Qw.*deltauvec_Qw;

deltauvec2_SO4=deltauvec_SO4.*timevector;
ISE_SO4deltauvec=deltauvecsquare_SO4.*timevector;
deltauvec2_Qw=deltauvec_Qw.*timevector;
ISE_Qwdeltauvec=deltauvecsquare_Qw.*timevector;

IAE_SO4 = sum(IAE_SO4vec);
ISE_SO4 = sum(ISE_SO4vec);
maxDEV_SO4 = max(abs(SO4error));
maxDEV_SO4u = max(SO4regpart(:,3))-min(SO4regpart(:,3));
IAE_Qw = sum(IAE_Qwvec);
ISE_Qw = sum(ISE_Qwvec);
maxDEV_Qw = max(abs(Qwerror));
maxDEV_Qwu = max(Qwregpart(:,1))-min(Qwregpart(:,1));

SO4errormean = sum(SO4vec)/totalt;
SO4errormeansquare = ISE_SO4/totalt;
SO4_deltau_errormean = sum(deltauvec2_SO4)/totalt;
SO4_deltau_errormeansquare = sum(ISE_SO4deltauvec)/totalt;
Qwerrormean = sum(Qwvec)/totalt;
Qwerrormeansquare = ISE_Qw/totalt;
Qw_deltau_errormean = sum(deltauvec2_Qw)/totalt;
Qw_deltau_errormeansquare = sum(ISE_Qwdeltauvec)/totalt;

errorvar_SO4 = SO4errormeansquare - (SO4errormean.*SO4errormean);
errorvar_deltau_SO4 = SO4_deltau_errormeansquare - (SO4_deltau_errormean.*SO4_deltau_errormean);
errorvar_Qw = Qwerrormeansquare - (Qwerrormean.*Qwerrormean);
errorvar_deltau_Qw = Qw_deltau_errormeansquare - (Qw_deltau_errormean.*Qw_deltau_errormean);

disp(' ')
disp(['Performance of active controllers during time ',num2str(time_eval(1)),' to ',num2str(time_eval(end)),' days'])
disp('*************************************************************')
disp(' ')
disp('Sludge wastage flow rate controller (timer based)')
disp('=================================================')
disp(' ')
disp(['Time 0 to 182 days: Qw = ',num2str(Qw_low),' m3/d'])
disp(['Time 183 to 364 days: Qw = ',num2str(Qw_high),' m3/d'])
disp(['Time 365 to 546 days: Qw = ',num2str(Qw_low),' m3/d'])
disp(['Time 547 to 609 days: Qw = ',num2str(Qw_high),' m3/d'])
disp(' ')
disp(['Controlled variable = Qw, setpoint = see above'])
disp('-------------------------------------------------------')
disp(['Average value of error (mean(e)) = ',num2str(mean(Qwerror)),' (m3/d)'])
disp(['Average value of absolute error (mean(|e|)) = ',num2str(Qwmean),' (m3/d)'])
disp(['Integral of absolute error (IAE) = ',num2str(IAE_Qw),' m3'])
disp(['Integral of square error (ISE) = ',num2str(ISE_Qw),' (m3^2)/d'])
disp(['Maximum absolute deviation from timer based setpoint (max(e)) = ',num2str(maxDEV_Qw),' m3/d'])
disp(['Standard deviation of error (std(e)) = ',num2str(sqrt(errorvar_Qw)),' m3/d'])
disp(['Variance of error (var(e)) = ',num2str(errorvar_Qw),' (m3/d)^2'])
disp(' ')
disp('Manipulated variable (MV), Qw')
disp('--------------------------------')
disp(['Maximum absolute variation of MV (max-min) = ',num2str(maxDEV_Qwu),' m3/d'])
disp(['Maximum absolute variation of MV in one sample (max delta) = ',num2str(max(diff(Qwregpart(:,1)))),' m3/d'])
disp(['Average value of MV (mean(Qw)) = ',num2str(mean(Qwregpart(:,1))),' m3/d'])
disp(['Standard deviation of MV (std(delta(Qw))) = ',num2str(sqrt(errorvar_deltau_Qw)),' m3/d'])
disp(['Variance of MV (var(delta(Qw))) = ',num2str(errorvar_deltau_Qw),' (m3/d)^2'])
disp(' ')
disp(' ')
disp('Oxygen controller for second aerobic reactor (reactor 4)')
disp('========================================================')
disp(' ')
disp(['PI controller with anti-windup: K = ',num2str(KSO4),' 1/d/(g (-COD)/m3)'])
disp(['                                Ti = ',num2str(TiSO4),' days'])
disp(['                                Tt = ',num2str(TtSO4),' days'])
disp('NOTE: same KLa value used for reactor 3 DO control and KLa4/2 used for reactor 5 DO control')
disp(' ')
disp(['Controlled variable = SO (reactor 4), setpoint = ',num2str(SO4ref),' mg (-COD)/l'])
disp('--------------------------------------------------------------')
disp(' ')
disp(['Average value of error (mean(e)) = ',num2str(mean(SO4error)),' (mg (-COD)/l)'])
disp(['Average value of absolute error (mean(|e|)) = ',num2str(SO4mean),' (mg (-COD)/l)'])
disp(['Integral of absolute error (IAE) = ',num2str(IAE_SO4),' (mg (-COD)/l)*d'])
disp(['Integral of square error (ISE) = ',num2str(ISE_SO4),' (mg (-COD)/l)^2*d'])
disp(['Maximum absolute deviation from oxygen setpoint (max(e)) = ',num2str(maxDEV_SO4),' mg (-COD)/l'])
disp(['Standard deviation of error (std(e)) = ',num2str(sqrt(errorvar_SO4)),' mg (-COD)/l'])
disp(['Variance of error (var(e)) = ',num2str(errorvar_SO4),' (mg (-COD)/l)^2'])
disp(' ')
disp('Manipulated variable (MV), KLa (reactor 4)')
disp('------------------------------------------')
disp('Based on KLa controller output prior to KLa actuator model')
disp(['Maximum absolute variation of MV (max-min) = ',num2str(maxDEV_SO4u),' 1/d'])
disp(['Maximum absolute variation of MV in one sample (max delta) = ',num2str(max(diff(SO4regpart(:,3)))),' 1/d'])
disp(['Average value of MV (mean(KLa4)) = ',num2str(mean(SO4regpart(:,3))),' 1/d'])
disp(['Standard deviation of MV (std(delta(KLa4))) = ',num2str(sqrt(errorvar_deltau_SO4)),' 1/d'])
disp(['Variance of MV (var(delta(KLa4))) = ',num2str(errorvar_deltau_SO4),' (1/d)^2'])
disp(' ')

disp(' ')
disp('Plotting of BSM2 controller evaluation results has been initiated')
disp('*****************************************************************')
disp(' ')
movingaveragewindow = 96; % even number
timeshift = movingaveragewindow/2;
b = ones(1,movingaveragewindow)./movingaveragewindow;
    
figure(23)
plot(time_eval(1:(end-1)),Qwerror)
hold on
xlabel('time (days)')
ylabel('error in Qw control (m^3/d)')
title('Qwref - Qwtruevalue')
hold off
xlim([starttime stoptime])

figure(24)
plot(time_eval(1:(end-1)),settlerpart(:,22))
hold on
plot(time_eval(1:(end-1)),Qwrefvec(:),'r--')
xlabel('time (days)')
ylabel('Qw control (m^3/d)')
title('Qwref (red), Qwtruevalue (blue)')
hold off
xlim([starttime stoptime])

figure(25)
plot(time_eval,Qwregpart(:,1), 'r')
hold on 
plot(time_eval,Qwregpart(:,2), 'b')
xlabel('time (days)')
ylabel('Wastage flow rate output (m^3/d)')
title('Qwcontrolleroutput (red), Qwactuatoroutput (blue)')
hold off
xlim([starttime stoptime])

figure(26)
plot(time_eval(1:(end-1)),SO4error)
hold on
filteredout=filter(b,1,SO4error);
filteredout=filteredout(movingaveragewindow:end);
plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'g','LineWidth',1)
xlabel('time (days)')
ylabel('error in SO4 control (mg (-COD)/l)')
title('SO4ref - SO4truevalue (raw and filtered)')
hold off
xlim([starttime stoptime])

figure(27)
plot(time_eval,SO4sensorpart(:,6),'g')
hold on
plot(time_eval(1:(end-1)),reac4part(:,8))
filteredout=filter(b,1,SO4sensorpart(:,6));
filteredout=filteredout(movingaveragewindow:end);
plot(time_eval(timeshift:(end-timeshift)),filteredout,'g--','LineWidth',1)
filteredout=filter(b,1,reac4part(:,8));
filteredout=filteredout(movingaveragewindow:end);
plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'y','LineWidth',1)
plot(time_eval(1:(end-1)),SO4refvec(:),'r--')
xlabel('time (days)')
ylabel('SO4 control (mg (-COD)/l)')
title('SO4ref (red), SO4truevalue (blue(raw)/yellow(filt)), SO4sensorvalue (green)')
hold off
xlim([starttime stoptime])

figure(28)
plot(time_eval,SO4regpart(:,3), 'r')
hold on 
plot(time_eval,kla4inpart, 'b')
filteredout=filter(b,1,SO4regpart(:,3));
filteredout=filteredout(movingaveragewindow:end);
plot(time_eval(timeshift:(end-timeshift)),filteredout,'g','LineWidth',1)
filteredout=filter(b,1,kla4inpart);
filteredout=filteredout(movingaveragewindow:end);
plot(time_eval(timeshift:(end-timeshift)),filteredout,'y','LineWidth',1)
xlabel('time (days)')
ylabel('KLa4 output (1/d)')
title('KLa4controlleroutput (red/green), KLa4actuatoroutput (blue/yellow)')
hold off
xlim([starttime stoptime])

figure(29)
plot(time_eval(1:(end-1)),reac3part(:,8))
hold on
filteredout=filter(b,1,reac3part(:,8));
filteredout=filteredout(movingaveragewindow:end);
plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'r','LineWidth',1)
xlabel('time (days)')
ylabel('SO3 (mg (-COD)/l)')
title('Resulting SO concentration in reactor 3 (input = KLa4), (raw and filtered)')
hold off
xlim([starttime stoptime])

figure(30)
plot(time_eval(1:(end-1)),reac5part(:,8))
hold on
filteredout=filter(b,1,reac5part(:,8));
filteredout=filteredout(movingaveragewindow:end);
plot(time_eval(timeshift:(end-timeshift-1)),filteredout,'r','LineWidth',1)
xlabel('time (days)')
ylabel('SO5 (mg (-COD)/l)')
title('Resulting SO concentration in reactor 5 (input = KLa4/2), (raw and filtered)')
xlim([starttime stoptime])

disp('Plotting of BSM2 controller evaluation results has been completed')
disp('*****************************************************************')
disp(' ')

stop=clock;
disp('***** Evaluation of BSM2 controllers successfully finished *****')
disp(['End time (hour:min:sec) = ', num2str(round(stop(4:6)))]); %Display simulation stop time
disp(' ')
