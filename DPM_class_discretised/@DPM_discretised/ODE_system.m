function [fun] = ODE_system(O)

%% Initial values
% Local time (initial time)
t = O.sol_time(1);

% Initial time
TIME0 = t; % is a vector

% Initial concentration
c = O.in_conc;

% Initial temperature
T = O.Tprofile(TIME0);

% Initial supersaturation
cs = evalanonfunc(O.solubility,T);

% Initial supersaturation
S = c/cs;

% Initial population
N = O.init_pop;
N_vector = N;
% Final population
nmax = O.max_pop;

% Crystals per bin
Nc = O.bin;
% Initial length
L = O.init_length(1);
% Dummy variable for crystal length
Ldummy = [];

% Initialise nucleation times
% O.nucleation_time = 1e50.*ones(O.max_pop+1,1);
O.nucleation_time(1,1) = 0;

% dummy vector of population
ndummy = [];

% Output matrix
Y = [L S t N c];

% Initialise the number of bins
N_of_bins = 1.*(N==Nc) + 1.*(mod(N,Nc) == 0).*(N>Nc);
safe_command = 0;

%% Integration
while(c(end)>= O.solubility(O.Tprofile(1e8)))
    
    %%
    Nbin = N;
    TIME_update = TIME0;
    
    if(safe_command == 0)
        
        while(Nbin-Nc*N_of_bins~=Nc && safe_command == 0)
            % Update the number of crystals
            [Ntemp,TIME,safe_command] = pgenerator_DPM(O,Nbin,TIME_update,L(end,:),N_of_bins);
            Nbin = Ntemp;
            O.nucleation_time(1+Nbin,1) = TIME;
            TIME_update = TIME;            
        end
        
    end
    
    %% Define the time span for the integrator
    if(TIME<O.sol_time(end) && safe_command == 0 && numel(O.nucleation_time)-1>1)
        % The system could generate the required amount of nuclei in
        % the desired time span
        interruption = 0;
        
        N_of_bins = N_of_bins + 1.*(mod(Nbin,Nc) == 0);
        N = numel(O.nucleation_time)-1;
        % Display the number of bins
        if(mod(N_of_bins,5) == 0)
        disp(['Bin n.: ' num2str(N_of_bins) ' completed'])
        end
                
        if(N_of_bins==1)
            t_i = O.nucleation_time(2);
        else            
            t_i = TIME0;
        end
        t_f = TIME;
        TIME0 = TIME;
        
        % Assure that the desired solution points are computed in the solution
        additional_times = O.sol_time(O.sol_time>=t_i & O.sol_time<=t_f);
        
        if(N==nmax)
            additional_times = [additional_times O.sol_time(O.sol_time>=t_f & O.sol_time<=max(O.sol_time))];
        end % if end
        
        additional_times;               
        
        time_span = sort(unique([additional_times O.nucleation_time((N_of_bins-1)*Nc+1:N_of_bins*Nc)' linspace(t_i,t_f,ceil(10+t_f-t_i)/3)]));
        time_span(time_span==0) = []; % Remove the zero, otherwise the code starts solving the ODEs of growth from this timepoint
        IntegrationODE = 1;
        
    elseif(TIME>O.sol_time(end) || safe_command == 1 && c(end) >= O.solubility(O.Tprofile(1e6)))
        % The system could not generate the required amount of nuclei in
        % the desired time span
        interruption = 1;
        
        N_of_bins = N_of_bins;
        N = numel(O.nucleation_time)-1;              
        
        t_i = TIME;
        
        t_f = max(O.sol_time);
        TIME0 = TIME;
        
        % Assure that the desired solution points are computed in the solution
        additional_times = O.sol_time(O.sol_time>=t_i & O.sol_time<=t_f);
        
        if(N==nmax)
            additional_times = [additional_times O.sol_time(O.sol_time>=t_f & O.sol_time<=max(O.sol_time))];
        end % if end
        
        time_span = sort(unique([additional_times linspace(t_i,t_f,ceil(10+t_f-t_i)/4)]));
        
        IntegrationODE = 1;
        
    elseif(numel(O.nucleation_time)-1==1 && safe_command == 0 && TIME<O.sol_time(end))
        % The system has generated only the first nucleus
        N_of_bins = N_of_bins;
        N = 1;
        
        interruption = 1;
        
        t_i = TIME; % Nucleation time 
        N_of_bins = 1; % Number of bins
        
        TIME0 = TIME;
        time_span = 0; % No time span
                
        IntegrationODE = 0; % Do not integrate any ODE
        

    end % if end

    %%
    switch(IntegrationODE)
        case 1 % Integrate the ODE system
            %% Formulate the new initial conditions
            if isempty(Ldummy)
                Ldummy0(1,1) = 0;
            elseif(interruption==0 && ~isempty(Ldummy))
                Ldummy0 = Ldummy(end,1:end);
                Ldummy0(1,end+1) = 0;
            elseif(interruption==1  && ~isempty(Ldummy))
                Ldummy0 = Ldummy(end,1:end);
            end
            %%
            % Clear the dummy variables
            clear cdummy Ldummy tdummy N_vectordummy additional_times
            
            % Compute the solution of the current ODEs system
            options = odeset('RelTol',1e-8,'BDF','on');
            [tdummy,Ldummy] = ode15s(@(tau,L) crystal_sizes(tau,L,O,N-1,N_of_bins),time_span,Ldummy0,options);
            % Concentration
            cdummy = O.in_conc - sum((Ldummy.^3).*O.bin,2).*O.kv.*O.rhoc./O.volume./O.rhos(O.Tprofile(1));
            % 3rd moment
            Moment3 = sum((Ldummy(end,:).^3).*O.bin,2).*O.kv./O.volume;            
            % Number of particles
            N_vectordummy = (N-1).*ones(numel(tdummy),1);
            %% Updating 
            if(N_of_bins-1>1 && interruption == 0)
                L(:,end+1)  = 0;
            else
                L  = L;
            end% if end
            
            %% Update the solution vectors
            N_vector = vertcat(N_vector(1:end-1),N_vectordummy);
            t = vertcat(t(1:end-1),tdummy);
            L = vertcat(L(1:end-1,:),Ldummy);
            c = vertcat(c(1:end-1),cdummy);
            
            %% Exit conditions
            if(Moment3>=1e-4)                
                index_a = find(sum((Ldummy.^3).*O.bin,2).*O.kv./O.volume>=1e-4,1,'First');
                index_b = find(sum((Ldummy.^3).*O.bin,2).*O.kv./O.volume<1e-4,1,'Last');
                t_a = tdummy(index_a);
                t_b = tdummy(index_b);
                disp('Crystals formed and detected')
                break
            elseif(time_span(end)>=O.sol_time(end) || TIME>O.sol_time(end) || abs(c(end) - O.solubility(O.Tprofile(1e6)))./O.solubility(O.Tprofile(1e6))<1e-3 || N>nmax)
                disp('The population required could not be achieved with the selected discretisation')
                break
            end % if end
            %%
        case 0 % No particles present, the system is not integrated
            
            if(N<Nc && N_of_bins ==0)
                disp('The population required could not be achieved with the selected discretisation')
                break
            elseif(time_span(end)>=O.sol_time(end) || TIME>O.sol_time(end) || abs(c(end) - O.solubility(O.Tprofile(1e6)))./O.solubility(O.Tprofile(1e6))<1e-3 || N>nmax)
                disp('The population required could not be achieved with the selected discretisation')
                break
            end % if end
            
    end
    
end % while end

% if(N<nmax)
%     disp('The population required could not be achieved in the selected time span')
% end % if end

%% Include the nucleation times and the other extra time in the solution times vector
try
    O.sol_time = unique(sort(vertcat(O.sol_time',O.nucleation_time,t_a,t_b)));
catch
    O.sol_time = unique(sort(vertcat(O.sol_time',O.nucleation_time)));
end
% Compute the supersaturation
S = c./O.solubility(O.Tprofile(t));
% Update the solution Matrix
Y = [L S t N_vector c];
%% Define outputs
O.calc_lengths = Y(:,1:size(L,2));
O.calc_conc = Y(:,end);
O.crystal_number = Y(:,end-1);
O.calc_time = Y(:,end-2);
O.calc_sup = Y(:,end-3);

end % function