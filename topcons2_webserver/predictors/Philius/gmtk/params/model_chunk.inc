%##############################################################################
%################## Copyright 2008 University of Washington ###################
%##############################################################################


    % -------------------------------------------------------------------------
    % changeState : in all frames beyond the 0th one, this variable is set
    % based on either stateCountDown(-1) or, if stateCountDown(0)=0 (which is
    % used to indicate that we are in an _exponential_ state), then changeState
    % will be chosen randomly based on a CPT that is a function of the current
    % state ...
    %
    % the switching parent which is in [0,K], is mapped down to just three
    % values using the mapping function "counter_map" as follows:
    %		if parent=0 --> 0
    %		if parent=1 --> 1
    %		if parent>1 --> 2
    % and these values in turn control which of the three conditionalparents
    % 'clauses' will be switched "on" ...

    variable: changeState {
	type: discrete hidden cardinality BINARY_CARD;
	weight: scale WEIGHT_changeState;
	switchingparents: stateCountDown(-1) using mapping ("counter_map");
	conditionalparents: 
		lengthClass(-1) using DenseCPT("leaveExpState_CPT") % 0
	      |	nil using DeterministicCPT("setBinaryTRUE")         % 1
	      | nil using DeterministicCPT("setBinaryFALSE");       % > 1
    }

    % -------------------------------------------------------------------------
    % topoChange : in all frames beyond the 0th one, this variable is set
    % by checking if the current topoClass is different from the previous
    % topoClass ...

    variable: topoChange {
	type: discrete hidden cardinality BINARY_CARD;
        switchingparents: topoClass(-1), topoClass(0) using mapping ("areDiff_map");
	conditionalparents: 
                  nil using DeterministicCPT("setBinaryFALSE")  % not diff
                | nil using DeterministicCPT("setBinaryTRUE");  % yes diff
    }

    variable: topo10 {
	type: discrete hidden cardinality TOPO_10_CARD;
	switchingparents: nil;
	conditionalparents: 
		topoChange(0), topoClass(0) using DeterministicCPT("two2oneMap");
    }

    % -------------------------------------------------------------------------
    % state : this variable indicates what state we are currently in ... 
    % in all frames after the 0th one, the current state will either be the
    % same as the previous state, OR if changeState=TRUE, it will be chosen
    % based on what is allowed to follow state(-1) 

    variable: state {
	type: discrete hidden cardinality STATE_CARD;
	weight: scale WEIGHT_state;
	switchingparents: changeState(0) using mapping ("oneD_map");
	conditionalparents: 
		state(-1) using DeterministicCPT("copyState")	% 0
	      | state(-1) using DenseCPT("stateTrans_CPT");     % 1
    }

    % -------------------------------------------------------------------------

    variable: spBit {
        type: discrete hidden cardinality BINARY_CARD;
        switchingparents: spBit(-1) using mapping ("oneD_map");
        conditionalparents:
                  topoClass(0) using DeterministicCPT("isSP")   % spBit(-1) = 0
                | nil using DeterministicCPT("setBinaryTRUE");  % spBit(-1) = 1
    }

    variable: iNumTM {
        type: discrete hidden cardinality NUM_TM_CARD;
        switchingparents: state(-1), state(0) using mapping ("areDiff_map");
        conditionalparents:
                  iNumTM(-1) using DeterministicCPT("copyNumTM")  % not diff
                | state(-1), iNumTM(-1) using DeterministicCPT("incNumTM");  % yes diff
    }

    % -------------------------------------------------------------------------
    % stateCountDown : this variable is used to count *down* the duration
    % of a particular state ... so the parents are :
    %	- from current frame : changeState (as a switching parent)
    %		if changeState=FALSE, decrement stateCountDown(-1)
    %		if changeState=TRUE,  choose duration for new state
    % note that the "decLength" decrements any incoming value that is greater
    % than 0 -- if the incoming (parent) value is 0, then the new value is
    % also 0 ...

    variable: stateCountDown {
	type: discrete hidden cardinality MAX_STATE_LENGTH_CARD;
	weight: scale WEIGHT_stateCountDown;
	switchingparents: changeState(0) using mapping ("oneD_map");
	conditionalparents: 
	    stateCountDown(-1) using DeterministicCPT("decLength")	% 0 
	  | lengthClass(0) using DenseCPT("length_CPT");		% 1 
    }

    % -------------------------------------------------------------------------
