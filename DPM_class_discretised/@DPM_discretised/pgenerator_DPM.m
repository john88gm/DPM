% Poisson process generator
function [n,t,safe_command] = pgenerator_DPM(O,number0,time0,Lin,N_of_bins)

% Define functions
V = O.volume;
J = O.nucleationrate;
T = O.Tprofile;
saturation = O.solubility;
T0 = O.Tprofile(0);

% Supersaturation with n-dependent and t-dependent concentration
S = @(T,t) instant_conc(time0,Lin,t,O,N_of_bins,number0)./saturation(T);

% n_lim = number0+1;
% i = 2;
% k = 0;
Time = O.sol_time(end);
%%
n = 0;
%%
l = @(t,V) V.*J(S(T(t),t),T(t));
L = @(t) integral(@(x) l(x,V),time0,t);

% Initialise
time = time0;
number = number0;
safe_command = 0;

t_temp = time0;

%%


r = rand;

% % for Aspirin
told = time+1000.*(time==0);
options = optimoptions('fsolve','Display','off','TolFun',1e-15);
[t_temporary,fval,exitflag] = fsolve(@(x) log(1-r) + L(x),told,options);
%   t_temp      disp('Evaluating trial...')
if(exitflag>0)
    t_temp = t_temporary;
else
    t_temp = time0;
end

% mini safe command
if(time>=Time && number == 0)
    disp(['Nothing happened'])
    n = number;
    t = time0;
    safe_command = 1;
elseif(time>=Time && number ~= 0)
    disp(['The process could not develop a new nucleus'])
    n = number;
    t = time0;
    safe_command = 1;
elseif(exitflag<=0)
    disp(['The process could not develop a new nucleus'])
    n = number;
    t = time0;
    safe_command = 1;
elseif(t_temp>time && t_temp<=Time)
    t = t_temp;
    n = number +1;
end % end if

end % function end