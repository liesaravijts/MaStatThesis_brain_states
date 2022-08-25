function state_dynamics=get_state_dynamics(T,spath,Gamma)
    % get_state_dynamics Calculates state summary metrics.
    %   state_dynamics=get_state_dynamics(T,spath,Gamma)
    %   
    %   FO: fractional occupancy [N x K]
    %   maxFO: maximum fractional occupancy [N x 1]
    %   SwitchingRate: rate of switching between states [N x 1]
    %   StateVisitsStartFinish: start and finish time points of state visits [N x K]
    %   StateLifeTimes: duration (number of time points) of state visits [1 x K]
    %   StateIntervals: intervals (number of time points) between state visits [1 x K]
    %   
    %   Makes use of the <a href="matlab:web('https://github.com/OHBA-analysis/HMM-MAR')">HMM-MAR Toolbox</a> (Vidaurre et al., 2016, 2017).
    %   
    %   See also run_hmm, run_gmm, run_swkmeans.
    
    state_dynamics = struct();
    
    if exist('Gamma','var')
        state_dynamics.FO = getFractionalOccupancy(Gamma,T);
        state_dynamics.maxFO = getMaxFractionalOccupancy(Gamma,T);
        state_dynamics.SwitchingRate = getSwitchingRate(Gamma,T); 
    else
        state_dynamics.FO = getFractionalOccupancy(spath,T);
        state_dynamics.maxFO = getMaxFractionalOccupancy(spath,T);
        state_dynamics.SwitchingRate = getSwitchingRate(spath,T);
    end
    
    state_dynamics.StateLifeTimes = getStateLifeTimes(spath,T);
    state_dynamics.StateIntervals = getStateIntervalTimes(spath,T);
    state_dynamics.StateVisitsStartFinish = getStartFinish_StateVisits(spath,T);
    
end
